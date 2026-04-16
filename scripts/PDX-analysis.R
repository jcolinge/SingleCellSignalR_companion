##################################################################
# Example code to process PDX samples submitted to bulk RNA-seq  #
# and perform subsequent LRI inference with SinfleCellSignalR v2 #
##################################################################

library(BulkSignalR)
library(SingleCellSignalR)
library(data.table)
library(ComplexHeatmap)
library(edgeR)
library(eulerr)


LRdb <- getResource("LRdb")
reactome <- getResource("Reactome")
gobp <- getResource("GO-BP")

if (.Platform$OS.type=="windows")
  windowsFonts("Arial" = windowsFont("Arial"))

res.folder <- "res1/"

# RNA-seq ======================================================================

counts.1 <- fread("data/hsa-counts.txt",sep="\t",data.table=F)
counts.2 <- fread("data/mmu-counts.txt",sep="\t",data.table=F)
sample.names <- sub("_align","",colnames(counts.1)[-(1:5)])


# protein coding

counts.1c <- counts.1[counts.1$type=="protein_coding",]
counts.1c <- counts.1c[!duplicated(counts.1c$symbol),]
write.table(counts.1c,file="protein-coding-counts-hsa.txt",quote=F,sep="\t",row.names=F)
counts.2c <- counts.2[counts.2$type=="protein_coding",]
counts.2c <- counts.2c[!duplicated(counts.2c$symbol),]
write.table(counts.2c,file="protein-coding-counts-mmu.txt",quote=F,sep="\t",row.names=F)

counts.h <- counts.1c[,-(1:5)]
rownames(counts.h) <- counts.1c[[2]]
colnames(counts.h) <- sample.names
good.h <- rowSums(counts.h>1)>=2
counts.h <- counts.h[good.h,]

counts.m <- counts.2c[,-(1:5)]
rownames(counts.m) <- counts.2c[[2]]
colnames(counts.m) <- paste0("mmu_",sample.names)
good.m <- rowSums(counts.m>1)>=2
counts.m <- counts.m[good.m,]


# normalization on the sum of human and murine alignments because of PDX

lib.size <- apply(counts.h,2,quantile,prob=0.75)+apply(counts.m,2,quantile,prob=0.75)
lib.size[1:2] <- sum(lib.size[1:2])
mls <- median(lib.size)
fact <- lib.size/mls
ncounts.h <- sweep(counts.h,2,fact,"/")
ncounts.m <- sweep(counts.m,2,fact,"/")
write.table(ncounts.h,file="normalized-counts-hsa.txt",quote=F,sep="\t")
write.table(ncounts.m,file="normalized-counts-mmu.txt",quote=F,sep="\t")


# Identify stromal / non-stromal parts =========================================

# Get the ortholog mapping using BSR functionality

ortholog.dict <- findOrthoGenes(
  from_organism = "mmusculus",
  from_values = rownames(ncounts.m)
)
ncounts.m2h <- convertToHuman(
  counts = ncounts.m,
  dictionary = ortholog.dict
)


# intersect on mapped orthologs and compute ratios

common <- intersect(rownames(ncounts.h),rownames(ncounts.m2h))
ncounts.hm <- cbind(ncounts.h[common,],ncounts.m2h[common,])
write.table(ncounts.hm,file="normalized-counts-hsa-mmu.txt",quote=F,sep="\t")
lratios <- (log1p(ncounts.h[common,])-log1p(ncounts.m2h[common,]))/log(10)


# ratios for the expression in human control sample, should display no mouse read counts

NL <- grep("Ctrl_h_",colnames(lratios))
med.ch <- lratios[,NL]
signal <- ncounts.h[common,NL]
ts <- 0
med.ch <- med.ch[signal>ts]
hist(med.ch,breaks=50)
thres.h.01 <- quantile(med.ch,prob=0.01)
abline(v=thres.h.01,col="red")
sum(med.ch<0)/length(med.ch) # 0.005822981
pdf(paste0(res.folder,"histogram-ratios-human-control.pdf"),width=4,height=2,useDingbats=F,pointsize=8)
par(mgp=c(1.75,0.7,0),mar=c(3,3,3,2))
hist(med.ch,breaks=25,main="human / mouse read counts",xlab="Log10[ (hsa+1) / (mmu+1) ]",xlim=c(-7,7))
abline(v=thres.h.01,col="red")
dev.off()


# ratios for the expression in mouse control sample, should display no human read counts

NC <- grep("Ctrl_m_",colnames(lratios))
med.cm <- lratios[,NC]
signal <- ncounts.m2h[common,NC]
ts <- 0
med.cm <- med.cm[signal>ts]
hist(med.cm,breaks=50)
thres.m.99 <- quantile(med.cm,prob=0.99)
sum(med.cm>0)/length(med.cm) # 0.006057745
abline(v=thres.m.99,col="blue")
pdf(paste0(res.folder,"histogram-ratios-mouse-control.pdf"),width=4,height=2,useDingbats=F,pointsize=8)
par(mgp=c(1.75,0.7,0),mar=c(3,3,3,2))
hist(med.cm,breaks=30,main="human / mouse read counts",xlab="Log10[ (hsa+1) / (mmu+1) ]",xlim=c(-7,7))
abline(v=thres.m.99,col="blue")
dev.off()


# ratios for the expression in PDX samples, should display both positive and negative ratios

PDX <- grep("^P",colnames(lratios))
med <- rowMeans(lratios[,PDX])
signal <- rowMeans(ncounts.h[common,PDX]) + rowMeans(ncounts.m2h[common,PDX])
ts <- 0
med <- med[signal>ts]
hist(med,breaks=50)
abline(v=0,col="orange")
pdf(paste0(res.folder,"histogram-ratios-pdx.pdf"),width=4,height=2,useDingbats=F,pointsize=8)
par(mgp=c(1.75,0.7,0),mar=c(3,3,3,2))
hist(med,breaks=30,main="human / mouse read counts",xlab="Log10[ (hsa+1) / (mmu+1) ]",xlim=c(-7,7))
abline(v=0,col="orange")
dev.off()


# Identify stromal & cancer cell genes ===========================================

# select stromal genes

signal <- rowMeans(ncounts.m2h[common,PDX])
pdx.ts <- 0
stromal <- rowMeans(lratios[,PDX]) < 0
is.stromal <- stromal & signal > pdx.ts
sum(is.stromal)
stromal.genes <- common[is.stromal]
write.table(stromal.genes,file=paste0(res.folder,"stromal-genes.txt"),row.names=F,quote=F,col.names=F)


# select cancer cell genes

signal <- rowMeans(ncounts.h[common,PDX])
pdx.ts <- 0
cancer <- rowMeans(lratios[,PDX]) > 0
is.cancer <- cancer & signal > pdx.ts
sum(is.cancer)
cancer.genes <- common[is.cancer]
write.table(cancer.genes,file=paste0(res.folder,"cancer-genes.txt"),row.names=F,quote=F,col.names=F)

intersect(stromal.genes,cancer.genes)


# Control heatmaps =============================================================

# cancer cell genes only

c.mat <- t(scale(t(data.matrix(ncounts.h[cancer.genes,]))))
dim(c.mat)
pdf(paste0(res.folder,"hm-cancer-genes.pdf"),width=3,height=5,pointsize=8)
simpleHeatmap(c.mat,row.names=F)
dev.off()


# stromal genes only

s.mat <- t(scale(t(data.matrix(ncounts.m2h[stromal.genes,]))))
dim(s.mat)
pdf(paste0(res.folder,"hm-stromal-genes.pdf"),width=3,height=5,pointsize=8)
simpleHeatmap(s.mat,row.names=F)
dev.off()


# both

cutExtremeValues <- function(m, p=0.01) {
  thres.lo <- stats::quantile(m, prob = p)
  thres.hi <- stats::quantile(m, prob = 1 - p)
  m[m > thres.hi] <- thres.hi
  m[m < thres.lo] <- thres.lo
  m
}

b.mat <- rbind(c.mat,s.mat)
di.spl <- dist(t(b.mat))
hc.spl <- hclust(di.spl, method = "ward.D")
dend.spl <- as.dendrogram(hc.spl)
di.gene <- dist(b.mat)
hc.gene <- hclust(di.gene, method = "ward.D")
dend.row <- as.dendrogram(hc.gene)
is.stromal <- rownames(b.mat)%in%stromal.genes
origin <- rep("cancer",length(is.stromal))
origin[is.stromal] <- "stroma"
col.origin <- setNames(c("red","green"),c("cancer","stroma"))
annot.genes <- rowAnnotation(origin=origin,
                       col=list(origin=col.origin)
)
b.mat.cut <- cutExtremeValues(b.mat)
cols <- circlize::colorRamp2(
  breaks = c(min(b.mat.cut), 0, max(b.mat.cut)),
  colors = c("royalblue3", "white", "orange")
)
clusters <- cutree(hc.spl,3)
clusters[clusters==1] <- "Ctrl"
clusters[clusters==2] <- "Cold"
clusters[clusters==3] <- "Hot"
clust.col <- setNames(c("lightgray","cyan","pink"),c("Ctrl","Cold","Hot"))
annot.samples <- columnAnnotation(group=clusters,
                       col=list(group=clust.col[clusters])
)
pdf(paste0(res.folder,"hm-cancer-and-stromal-genes.pdf"),width=4,height=7,pointsize=8)
Heatmap(b.mat.cut,
        cluster_rows = dend.row,
        cluster_columns = dend.spl, col = cols,
        show_row_names = FALSE,
        use_raster = TRUE, raster_device = "png",
        raster_quality = 8, raster_by_magick = FALSE,
        show_row_dend = TRUE,
        left_annotation = annot.genes,
        top_annotation = annot.samples
)
dev.off()


# ==============================================================================
# Analysis through the SCSR framework 
# ==============================================================================

# define "populations" that are here sample types actually

counts <- data.matrix(cbind(ncounts.h[common,],ncounts.m2h[common,]))
colnames(counts)
cond <- rep("PDX",ncol(ncounts.h))
cond[grep("^Ctrl_",names(ncounts.h))] <- "CTRL"
pop <- c(cond,paste0("mmu_",cond))
pop[colnames(counts)=="Ctrl_m_BAT"] <- "dummy"
pop[colnames(counts)=="mmu_Ctrl_h_hPanc"] <- "dummy"
pop


# SCSR object

scsr <- SCSRNet(counts=counts,
                populations=pop,
                normalize=FALSE, # already with special 2-species scheme
                method="UQ hsa+mmu",
                log.transformed=FALSE
)
scsr


# Define a PDX-specific pairwise comparison function adapted to bulk data based on edgeR

getDiffTable <- function(obj,pop){
  
  if (pop!="PDX" && pop!="mmu_PDX")
    stop("populations must be PDX (to be compared with CTRL in same species)")
  
  # define the comparison
  conditions <- populations(obj)
  if (pop == "PDX")
    ctrl <- "CTRL"
  else
    ctrl <- "mmu_CTRL"
  cl <- factor(conditions)
  design <- model.matrix(~0+cl)
  colnames(design) <- gsub("^cl","",colnames(design))
  comp <- paste(pop,"-",ctrl)
  cm <- makeContrasts(contrasts=comp,levels=design)
  
  # apply edgeR
  dge <- DGEList(ncounts(bsrdmComp(obj)),genes=rownames(ncounts(bsrdmComp(obj))))
  rownames(dge$counts) <- rownames(ncounts(bsrdmComp(obj)))
  # dge <- calcNormFactors(dge)  no edgeR normalization here !
  dge <- estimateDisp(dge, design, robust=T)
  fit.dge <- glmFit(dge, design)
  lrt <- glmLRT(fit.dge, contrast=cm[,comp])
  sel.r <- topTags(lrt, n=nrow(dge$counts))
  tab <- sel.r$table[,c("PValue","logFC","FDR")]
  tab$expr <- rowMeans(dge$counts[, which(conditions==pop)])
  colnames(tab)[1] <- "pval"
  rownames(tab) <- rownames(dge$counts)
  
  # kill the ligand or receptor that are not expressed in a given species
  if (pop == "PDX"){
    bad <- !(rownames(tab) %in% cancer.genes)
    tab$pval[bad] <- 1
    tab$logFC[bad] <- 0
  }
  else{
    bad <- !(rownames(tab) %in% stromal.genes)
    tab$pval[bad] <- 1
    tab$logFC[bad] <- 0
  }
  
  # return the BSRClusterComp object
  BSRClusterComp(bsrdmComp(obj),
                 which(conditions==pop),
                 which(conditions==ctrl),
                 tab
  )
  
} # getDiffTable


# Generate the needed comparisons (each species PDX versus corresponding CTRL)

scsr <- performInferences(scsr,
                          verbose=TRUE,
                          rank.p=0.75,
                          selected.populations=c("PDX",
                                                 "mmu_PDX"
                          ),
                          funDiffExpr=getDiffTable,
                          max.pval=0.05,
                          min.logFC=0.5
)
scsr
save(scsr,file="scsr.rdta")
load("scsr.rdta")


# stroma to cancer cell paracrine interactions

bsrinfc.stroma.cancer <- getParacrines(scsr,"mmu_PDX","PDX")
lri.stroma.cancer <- LRinter(bsrinfc.stroma.cancer)
write.table(lri.stroma.cancer,file=paste0(res.folder,"LRI-stroma-cancer.csv"),sep="\t",quote=F,row.names=F)
bsrinfc.stroma.cancer.red <- reduceToBestPathway(bsrinfc.stroma.cancer)
lri.stroma.cancer.red <- LRinter(bsrinfc.stroma.cancer.red)
bsrinfc.stroma.cancer.redP <- reduceToPathway(bsrinfc.stroma.cancer)
lri.stroma.cancer.redP <- LRinter(bsrinfc.stroma.cancer.redP)


# cancer cell to stroma paracrine interactions

bsrinfc.cancer.stroma <- getParacrines(scsr,"PDX","mmu_PDX")
lri.cancer.stroma <- LRinter(bsrinfc.cancer.stroma)
write.table(lri.cancer.stroma,file=paste0(res.folder,"LRI-cancer-stroma.csv"),sep="\t",quote=F,row.names=F)
bsrinfc.cancer.stroma.red <- reduceToBestPathway(bsrinfc.cancer.stroma)
lri.cancer.stroma.red <- LRinter(bsrinfc.cancer.stroma.red)
bsrinfc.cancer.stroma.redP <- reduceToPathway(bsrinfc.cancer.stroma)
lri.cancer.stroma.redP <- LRinter(bsrinfc.cancer.stroma.redP)


# stroma autocrine interactions

bsrinfc.a.stroma <- getAutocrines(scsr,"mmu_PDX")
lri.a.stroma <- LRinter(bsrinfc.a.stroma)
write.table(lri.a.stroma,file=paste0(res.folder,"LRI-auto-stroma.csv"),sep="\t",quote=F,row.names=F)
bsrinfc.a.stroma.red <- reduceToBestPathway(bsrinfc.a.stroma)
lri.a.stroma.red <- LRinter(bsrinfc.a.stroma.red)
bsrinfc.a.stroma.redP <- reduceToPathway(bsrinfc.a.stroma)
lri.a.stroma.redP <- LRinter(bsrinfc.a.stroma.redP)


# cancer cell to stroma paracrine interactions

bsrinfc.a.cancer <- getAutocrines(scsr,"PDX")
lri.a.cancer <- LRinter(bsrinfc.a.cancer)
write.table(lri.a.cancer,file=paste0(res.folder,"LRI-auto-cancer.csv"),sep="\t",quote=F,row.names=F)
bsrinfc.a.cancer.red <- reduceToBestPathway(bsrinfc.a.cancer)
lri.a.cancer.red <- LRinter(bsrinfc.a.cancer.red)
bsrinfc.a.cancer.redP <- reduceToPathway(bsrinfc.a.cancer)
lri.a.cancer.redP <- LRinter(bsrinfc.a.cancer.redP)


# graphic overview ===================================================


# ligands ---------------------------------

cancer.lig <- union(lri.cancer.stroma$L,lri.a.cancer$L)
stroma.lig <- union(lri.stroma.cancer$L,lri.a.stroma$L)

plot(venn(list(para=unique(lri.cancer.stroma$L),auto=unique(lri.a.cancer$L))))
plot(venn(list(para=unique(lri.stroma.cancer$L),auto=unique(lri.a.stroma$L))))


# global bubble plots ----------------------------

cellNetBubblePlot(scsr)
cellNetHeatmap(scsr, selected.populations=c("PDX","mmu_PDX"))

genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-190236"]) # Signaling by FGFR
cellNetBubblePlot(scsr, genes.to.count=genes)
pdf(paste0(res.folder,"bubble-FGFR-signal.pdf"),width=2.5,height=1.5,pointsize=7)
cellNetBubblePlot(scsr, genes.to.count=genes, use.proportions=TRUE)
dev.off()
cellNetHeatmap(scsr, genes.to.count=genes, use.proportions=TRUE, selected.populations=c("PDX","mmu_PDX"))

genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-194138"]) # Signaling by VEGF
cellNetBubblePlot(scsr, genes.to.count=genes)
pdf(paste0(res.folder,"bubble-VEGF-signal.pdf"),width=2.5,height=1.5,pointsize=7)
cellNetBubblePlot(scsr, genes.to.count=genes, use.proportions=TRUE)
dev.off()
cellNetHeatmap(scsr, genes.to.count=genes, use.proportions=TRUE, selected.populations=c("PDX","mmu_PDX"))



# pairwise bubble plots

lri <- LRinter(reduceToBestPathway(bsrinfc.stroma.cancer.redP))
pw.names <- unique(lri$pw.name[lri$len>=30])
bubblePlotPathwaysLR(bsrinfc.stroma.cancer,pw.names)

lri <- LRinter(reduceToBestPathway(bsrinfc.cancer.stroma.redP))
pw.names <- unique(lri$pw.name[lri$len>=10])
bubblePlotPathwaysLR(bsrinfc.cancer.stroma,pw.names)

lri <- LRinter(reduceToBestPathway(bsrinfc.a.stroma.redP))
pw.names <- unique(lri$pw.name[lri$len>=10])
bubblePlotPathwaysLR(bsrinfc.a.stroma,pw.names)


# signature heatmaps

sign.cancer.stroma.red <- BSRSignatureComp(bsrinfc.cancer.stroma.red,qval.thres=1e-6)
score.cancer.stroma.red <- scoreLRGeneSignatures(bsrdmComp(scsr),sign.cancer.stroma.red)[,-12]
annot.2 <- columnAnnotation(group=clusters[colnames(score.cancer.stroma.red)],
                            col=list(group=clust.col[clusters[colnames(score.cancer.stroma.red)]])
)
pdf(paste0(res.folder,"hm-cancer-stroma-red.pdf"),width=4,height=5,pointsize=8)
simpleHeatmap(score.cancer.stroma.red,bottom.annotation=annot.2)
dev.off()

sign.stroma.cancer.red <- BSRSignatureComp(bsrinfc.stroma.cancer.red,qval.thres=1e-6)
score.stroma.cancer.red <- scoreLRGeneSignatures(bsrdmComp(scsr),sign.stroma.cancer.red)[,-12]
annot.2 <- columnAnnotation(group=clusters[colnames(score.stroma.cancer.red)],
                            col=list(group=clust.col[clusters[colnames(score.stroma.cancer.red)]])
)
pdf(paste0(res.folder,"hm-stroma-cancer-red.pdf"),width=4,height=5,pointsize=8)
simpleHeatmap(score.stroma.cancer.red,bottom.annotation=annot.2)
dev.off()

sign.a.stroma.red <- BSRSignatureComp(bsrinfc.a.stroma.red,qval.thres=1e-6)
score.a.stroma.red <- scoreLRGeneSignatures(bsrdmComp(scsr),sign.a.stroma.red)[,-12]
annot.2 <- columnAnnotation(group=clusters[colnames(score.cancer.stroma.red)],
                            col=list(group=clust.col[clusters[colnames(score.cancer.stroma.red)]])
)
pdf(paste0(res.folder,"hm-auto-stroma-red.pdf"),width=4,height=5,pointsize=8)
simpleHeatmap(score.a.stroma.red,bottom.annotation=annot.2)
dev.off()

sign.a.cancer.red <- BSRSignatureComp(bsrinfc.a.cancer.red,qval.thres=1e-6)
score.a.cancer.red <- scoreLRGeneSignatures(bsrdmComp(scsr),sign.a.cancer.red)[,-12]
annot.2 <- columnAnnotation(group=clusters[colnames(score.cancer.stroma.red)],
                            col=list(group=clust.col[clusters[colnames(score.cancer.stroma.red)]])
)
pdf(paste0(res.folder,"hm-auto-cancer-red.pdf"),width=4,height=5,pointsize=8)
simpleHeatmap(score.a.cancer.red,bottom.annotation=annot.2)
dev.off()


# Checks =======================================================================

tab.stroma <- differentialStats(comparison(bsrdmComp(scsr))[["mmu_PDX_vs_others"]])
tab.cancer <- differentialStats(comparison(bsrdmComp(scsr))[["PDX_vs_others"]])
tab.stroma[cancer.genes,"pval"]
tab.cancer[stromal.genes,"pval"]

auto.c <- unique(paste(lri.a.cancer$L,lri.a.cancer$R,sep="¦"))
auto.s <- unique(paste(lri.a.stroma$L,lri.a.stroma$R,sep="¦"))
intersect(auto.c,auto.s)

para.cs <- unique(paste(lri.cancer.stroma$L,lri.cancer.stroma$R,sep="¦"))
para.sc <- unique(paste(lri.stroma.cancer$L,lri.stroma.cancer$R,sep="¦"))
intersect(para.cs,para.sc)

intersect(auto.c,para.sc)
intersect(auto.s,para.sc)
intersect(auto.c,para.cs)
intersect(auto.s,para.cs)

intersect(stroma.lig,cancer.lig)

