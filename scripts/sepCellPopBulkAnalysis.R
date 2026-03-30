###################################################################
# Example code to process FACS-sorted populations un bulk RNA-seq #
# and perform subsequent LRI inference with SingleCellSignalR v2  #
###################################################################

library(BulkSignalR)
library(SingleCellSignalR)
library(data.table)
library(matrixStats)
library(edgeR)
library(eulerr)
library(igraph)

if (.Platform$OS.type=="windows")
  windowsFonts("Arial" = windowsFont("Arial"))


LRdb <- getResource("LRdb")
LRdb.clean <- LRdb
LRdb.clean <- LRdb.clean[!(LRdb.clean$ligand %in% c("ARRB2","JAK2","NCK1","PSG1",
                                                    "PTPN11","PTPN6","SRC","STAT3",
                                                    "SYK","CD14","HRAS","HSP90AA1",
                                                    "HSP90B1")),]
resetLRdb(LRdb.clean,switch=TRUE)
LRdb <- getResource("LRdb")
reactome <- getResource("Reactome")
gobp <- getResource("GO-BP")


# sample description

sample.descr <- read.csv("sample-descr.txt",sep="\t")
cell.type <- vapply(strsplit(unlist(sample.descr["Sample_characteristics_ch1",])," "),
                      function(x) x[4], "text"
)
names(cell.type) <- names(sample.descr)
sample.type <- gsub("[0-9]+$","",names(sample.descr))
names(sample.type) <- names(sample.descr)
SRR <- fread("map-SRR.txt",data.table=FALSE,header=FALSE)
srr2gsm <- setNames(SRR[[2]],SRR[[1]])

# transcriptomics

mat <- fread("sepppop-counts.txt",data.table=FALSE)
mat <- mat[mat$type=="protein_coding",]
good <- rowSums2(mat[,-(1:5)]>=5)>=3
sum(good)
mat <- mat[good,]
counts <- data.matrix(mat[,-(1:5)])
rownames(counts) <- mat$symbol
colnames(counts) <- srr2gsm[gsub("_align$","",colnames(counts))]
dup <- rownames(counts)[duplicated(rownames(counts))]
counts[dup,]
counts[duplicated(rownames(counts)),]
counts <- counts[!duplicated(rownames(counts)),]

# ==============================================================================
# Use SCSR to perform analyses similar to the original paper by Choi, et al.
# ==============================================================================


# human ortholog genes using BSR mechanism

ortholog.dict <- findOrthoGenes(
  from_organism = "mmusculus",
  from_values = rownames(counts)
)
counts.h <- convertToHuman(
  counts = counts,
  dictionary = ortholog.dict
)


# Create the SCSRNet object

pop <- paste(cell.type,sample.type,sep="_")
scsr <- SCSRNet(counts=counts.h,
                populations=pop,
                normalize=FALSE, # because of edgeR diff analysis
                method="nonorm",
                log.transformed=FALSE
)
scsr


# Define a problem-specific pairwise comparison function adapted to bulk data based on edgeR

getDiffTable <- function(obj,pop){
  
  if (length(grep("_Tum$",pop)) == 0)
    stop("populations must be tumor populations (to be compared with WT counterparts)")

  # define the comparison
  conditions <- populations(obj)
  wt.pop <- gsub("_Tum$","_WT",pop)
  cl <- factor(conditions)
  design <- model.matrix(~0+cl)
  colnames(design) <- gsub("^cl","",colnames(design))
  comp <- paste(pop,"-",wt.pop)
  cm <- makeContrasts(contrasts=comp,levels=design)

  # apply edgeR
  dge <- DGEList(ncounts(bsrdmComp(obj)),genes=rownames(ncounts(bsrdmComp(obj))))
  rownames(dge$counts) <- rownames(ncounts(bsrdmComp(obj)))
  dge <- calcNormFactors(dge)  
  dge <- estimateDisp(dge, design, robust=T)
  fit.dge <- glmFit(dge, design)
  lrt <- glmLRT(fit.dge, contrast=cm[,comp])
  sel.r <- topTags(lrt, n=nrow(dge$counts))
  tab <- sel.r$table[,c("PValue","logFC")]
  tab$expr <- rowMeans(dge$counts[, which(conditions==pop)])
  colnames(tab)[1] <- "pval"
  rownames(tab) <- rownames(dge$counts)
  
  # return the BSRClusterComp object
  BSRClusterComp(bsrdmComp(obj),
                 which(conditions==pop),
                 which(conditions==wt.pop),
                 tab
  )
  
} # getDiffTable


# Generate the needed comparisons (each immune pop tumor versus WT)

scsr <- performInferences(scsr,
                          verbose=TRUE,
                          rank.p=0.75,
                          selected.populations=c("macrophages_Tum",
                                                 "neutrophils_Tum",
                                                 "monocytes_Tum",
                                                 "epithelial_Tum"
                          ),
                          funDiffExpr=getDiffTable,
                          max.pval=0.05,
                          min.logFC=0.5
)
scsr
save(scsr,file="scsr.rdta")


# macrophage to cancer cell paracrine interactions =================================

# bsrinfc.macro.epi <- resetToInitialOrganism(getParacrines(scsr,"macrophages_Tum","epithelial_Tum"),
#                                             conversion.dict=ortholog.dict
# )
bsrinfc.macro.epi <- getParacrines(scsr,"macrophages_Tum","epithelial_Tum")
lri.macro.epi <- LRinter(bsrinfc.macro.epi)
write.table(lri.macro.epi[,1:8],file="LRI-macro-epi.txt",sep="\t",quote=F,row.names=F)
bsrinfc.macro.epi.red <- reduceToBestPathway(bsrinfc.macro.epi)
lri.macro.epi.red <- LRinter(bsrinfc.macro.epi.red)
write.table(lri.macro.epi.red[,1:8],file="LRI-macro-epi-best-pw.txt",sep="\t",quote=F,row.names=F)
bsrinfc.macro.epi.redP <- reduceToPathway(bsrinfc.macro.epi)
lri.macro.epi.redP <- LRinter(bsrinfc.macro.epi.redP)


# monocyte to cancer cell paracrine interactions =================================

bsrinfc.mono.epi <- getParacrines(scsr,"monocytes_Tum","epithelial_Tum")
lri.mono.epi <- LRinter(bsrinfc.mono.epi)
write.table(lri.mono.epi[,1:8],file="LRI-mono-epi.txt",sep="\t",quote=F,row.names=F)
bsrinfc.mono.epi.red <- reduceToBestPathway(bsrinfc.mono.epi)
lri.mono.epi.red <- LRinter(bsrinfc.mono.epi.red)
write.table(lri.mono.epi.red[,1:8],file="LRI-mono-epi-best-pw.txt",sep="\t",quote=F,row.names=F)
bsrinfc.mono.epi.redP <- reduceToPathway(bsrinfc.mono.epi)
lri.mono.epi.redP <- LRinter(bsrinfc.mono.epi.redP)

tab.mono <- differentialStats(comparison(bsrdmComp(scsr))[["monocytes_Tum_vs_others"]])
sort(intersect(LRdb$ligand,rownames(tab.mono)[tab.mono$pval<0.001]))


# neutrophil to cancer cell paracrine interactions =================================

bsrinfc.neutro.epi <- getParacrines(scsr,"neutrophils_Tum","epithelial_Tum")
lri.neutro.epi <- LRinter(bsrinfc.neutro.epi)
write.table(lri.neutro.epi[,1:8],file="LRI-neutro-epi.txt",sep="\t",quote=F,row.names=F)
bsrinfc.neutro.epi.red <- reduceToBestPathway(bsrinfc.neutro.epi)
lri.neutro.epi.red <- LRinter(bsrinfc.neutro.epi.red)
write.table(lri.neutro.epi[,1:8],file="LRI-neutro-epi-best-pw.txt",sep="\t",quote=F,row.names=F)
bsrinfc.neutro.epi.redP <- reduceToPathway(bsrinfc.neutro.epi)
lri.neutro.epi.redP <- LRinter(bsrinfc.neutro.epi.redP)

tab.neutro <- differentialStats(comparison(bsrdmComp(scsr))[["neutrophils_Tum_vs_others"]])
sort(intersect(LRdb$ligand,rownames(tab.neutro)[tab.neutro$pval<0.001]))


# cancer cell autocrine interactions =================================

bsrinfc.epi.auto <- getAutocrines(scsr,"epithelial_Tum")
lri.epi.auto <- LRinter(bsrinfc.epi.auto)
write.table(lri.epi.auto[,1:8],file="LRI-epi-autocrine.txt",sep="\t",quote=F,row.names=F)
bsrinfc.epi.auto.red <- reduceToBestPathway(bsrinfc.epi.auto)
lri.epi.auto.red <- LRinter(bsrinfc.epi.auto.red)
write.table(lri.epi.auto.red[,1:8],file="LRI-epi-autocrine-best-pw.txt",sep="\t",quote=F,row.names=F)
bsrinfc.epi.auto.redP <- reduceToPathway(bsrinfc.epi.auto)
lri.epi.auto.redP <- LRinter(bsrinfc.epi.auto.redP)


# graphic overview ===================================================


# ligands ---------------------------------

lig <- list(macro=unique(lri.macro.epi$L),
            mono=unique(lri.mono.epi$L),
            neutro=unique(lri.neutro.epi$L),
            epi=unique(lri.epi.auto$L)
)

fills = list(fill = rainbow(4)[c(4:1)], alpha = 0.3)
# fills = list(fill = hcl.colors(4, palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE),alpha=0.4)
plot(venn(lig),fills=fills)

pdf("euler-ligands.pdf",width=3.4,height=2.3,pointsize=8,useDingbats=FALSE)  
plot(venn(lig),fills=fills)
dev.off()

nonauto.lig <- setdiff(unlist(lig[1:3]),lig$epi)
macro.lig <- setdiff(lig$macro,unlist(lig[2:4]))
mono.lig <- setdiff(lig$mono,unlist(lig[c(1,3,4)]))
neutro.lig <- setdiff(lig$neutro,unlist(lig[c(1,2,4)]))
epi.lig <- setdiff(lig$epi,unlist(lig[1:3]))


# bubble plots ----------------------------

pdf("global-LRIs.pdf",width=3.2,height=2.3,pointsize=8,useDingbats=FALSE)  
cellNetBubblePlot(scsr)
dev.off()

genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-449147"]) # Signaling by Interleukins
pdf("interleukin-LRIs.pdf",width=3.4,height=2.3,pointsize=8,useDingbats=FALSE)  
cellNetBubblePlot(scsr, genes.to.count=genes, use.proportions=TRUE)
dev.off()

genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-9006936"]) # Signaling by TGFB
pdf("TGFb-LRIs.pdf",width=3.4,height=2.3,pointsize=8,useDingbats=FALSE)  
cellNetBubblePlot(scsr, genes.to.count=genes, use.proportions=TRUE)
dev.off()

genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-195721"]) # Signaling by WNT
pdf("WNT-LRIs.pdf",width=3.4,height=2.3,pointsize=8,useDingbats=FALSE)  
cellNetBubblePlot(scsr, genes.to.count=genes, use.proportions=TRUE)
dev.off()

genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-1474244"]) # ECM reorg
pdf("ECM-LRIs.pdf",width=3.4,height=2.3,pointsize=8,useDingbats=FALSE)  
cellNetBubblePlot(scsr, genes.to.count=genes, use.proportions=TRUE)
dev.off()


# interleukin network -------------------------------------------------------

# IL at macrophages
index <- grep("^IL",lri.macro.epi.red$L)
IL.macro <- setdiff(lri.macro.epi.red$L[index],"IL6ST")

# ILR at epithelial cells
index <- grepl("^IL",lri.macro.epi.red$R)
ILR.epi <- setdiff(lri.macro.epi.red$R[lri.macro.epi.red$L %in% IL.macro & index],"IL6ST")

# IL-ILR from the macro -> epithelial LRIs
good <- lri.macro.epi.red$L %in% IL.macro & lri.macro.epi.red$R %in% ILR.epi
links.para <- lri.macro.epi.red[good,1:2]

# finds IL & ILR targets of the above LRIs, excluding ILR.epi
red.targets <- tgGenes(bsrinfc.macro.epi.red)[good]
red.t.pvals <- tgPval(bsrinfc.macro.epi.red)[good]
red.t.logFC <- tgLogFC(bsrinfc.macro.epi.red)[good]
rec.para <- lri.macro.epi.red$R[good]
targets.para <- NULL
for (i in seq_len(length(red.targets))){
  signif <- red.targets[[i]][red.t.pvals[[i]]<0.05 & red.t.logFC[[i]]>0.5]
  select <- setdiff(signif[grep("^IL",signif)],c(ILR.epi,"IL6ST"))
  if (length(select)>0)
    targets.para <- rbind(targets.para,data.frame(R=rep(rec.para[i],length(select)),
                                                  T=select))
}

# IL & ILR up that would be triggered by other receptors
red.targets <- tgGenes(bsrinfc.macro.epi.red)[!good & lri.macro.epi.red$R!="IL6ST"]
red.t.pvals <- tgPval(bsrinfc.macro.epi.red)[!good & lri.macro.epi.red$R!="IL6ST"]
red.t.logFC <- tgLogFC(bsrinfc.macro.epi.red)[!good & lri.macro.epi.red$R!="IL6ST"]
all.targets <- NULL
for (i in seq_len(length(red.targets)))
  all.targets <- c(all.targets,red.targets[[i]][red.t.pvals[[i]]<0.05 & red.t.logFC[[i]]>0.5])
all.targets <- unique(all.targets)
ILR.others <- setdiff(all.targets[grep("^IL",all.targets)],c(ILR.epi,"IL6ST",targets.para$T))
ILR.others # ==> none !

# autocrines
links.auto <- lri.epi.auto.red[grepl("^IL",lri.epi.auto.red$L) & grepl("^IL",lri.epi.auto.red$R),1:2]
links.auto <- links.auto[links.auto$L!="IL6ST" & links.auto$R!="IL6ST",]

# build the graph

# paracrine LRIs
links <- data.frame(from=paste0(links.para$L,"_mac"),
                    to=paste0(links.para$R,"_epi1"),
                    type=rep("macro_epi",nrow(links.para))
)

# epithelial IL receptors towards their IL & ILR targets
links <- rbind(links,data.frame(from=paste0(targets.para$R,"_epi1"),
                                to=paste0(targets.para$T,"_epi2"),
                                type=rep("intra_epi",nrow(targets.para))
))

# autocrine at epithelial
links <- rbind(links,data.frame(from=paste0(links.auto$L,"_epi2"),
                                to=paste0(links.auto$R,"_epi3"),
                                type=rep("epi_auto",nrow(links.auto))
))

# the network itself
g.interleukin <- graph_from_data_frame(unique(links))
plot(g.interleukin)
write_graph(g.interleukin,file="graph-macro-cancer-interleukin.graphml",format="graphml")
molecules <- union(links$from,links$to)
write.table(data.frame(node=molecules,label=gsub("_.+$","",molecules)),
            file="dic-macro-cancer-interleukin.txt",sep="\t",quote=F,row.names=F)


# WNT network -------------------------------------------------------

# macro
macro.wnt <- unique(lri.macro.epi[lri.macro.epi$pw.id=="R-HSA-195721",1:2])
good <- lri.macro.epi$pw.id=="R-HSA-195721"
red.targets <- tgGenes(bsrinfc.macro.epi)[good]
red.t.pvals <- tgPval(bsrinfc.macro.epi)[good]
red.t.logFC <- tgLogFC(bsrinfc.macro.epi)[good]
targets.macro <- NULL
for (i in seq_len(length(red.targets))){
  signif <- red.targets[[i]][red.t.pvals[[i]]<0.05 & red.t.logFC[[i]]>0.5]
  select <- signif[grep("^(WNT|FZD)",signif)]
  if (length(select)>0)
    targets.macro <- c(targets.macro,select)
}
unique(targets.macro)

# mono
mono.wnt <- unique(lri.mono.epi[lri.mono.epi$pw.id=="R-HSA-195721",1:2])
mono.wnt <- mono.wnt[!mono.wnt$L %in% c("APOE","CDH1"),]
good <- lri.mono.epi$pw.id=="R-HSA-195721"
red.targets <- tgGenes(bsrinfc.mono.epi)[good]
red.t.pvals <- tgPval(bsrinfc.mono.epi)[good]
red.t.logFC <- tgLogFC(bsrinfc.mono.epi)[good]
targets.mono <- NULL
for (i in seq_len(length(red.targets))){
  signif <- red.targets[[i]][red.t.pvals[[i]]<0.05 & red.t.logFC[[i]]>0.5]
  select <- signif[grep("^(WNT|FZD)",signif)]
  if (length(select)>0)
    targets.mono <- c(targets.mono,select)
}
unique(targets.mono)

# neutro
neutro.wnt <- unique(lri.neutro.epi[lri.neutro.epi$pw.id=="R-HSA-195721",1:2])
neutro.wnt <- neutro.wnt[!neutro.wnt$L %in% c("APOE","CDH1"),]
good <- lri.neutro.epi$pw.id=="R-HSA-195721"
red.targets <- tgGenes(bsrinfc.neutro.epi)[good]
red.t.pvals <- tgPval(bsrinfc.neutro.epi)[good]
red.t.logFC <- tgLogFC(bsrinfc.neutro.epi)[good]
targets.neutro <- NULL
for (i in seq_len(length(red.targets))){
  signif <- red.targets[[i]][red.t.pvals[[i]]<0.05 & red.t.logFC[[i]]>0.5]
  select <- signif[grep("^(WNT|FZD)",signif)]
  if (length(select)>0)
    targets.neutro <- c(targets.neutro,select)
}
unique(targets.neutro)


# ==============================================================================
# No network version, good old LR-score only
# ==============================================================================

pop <- paste(cell.type,sample.type,sep="_")
scsr.no <- SCSRNoNet(counts=counts.h,
                populations=pop,
                normalize=FALSE, # because of edgeR diff analysis
                method="nonorm",
                log.transformed=FALSE
)
scsr.no <- performInferences(scsr.no,
                             verbose=TRUE,
                             selected.populations=c("macrophages_Tum",
                                                    "neutrophils_Tum",
                                                    "monocytes_Tum",
                                                    "epithelial_Tum"
                             ),
                             funDiffExpr=getDiffTable,
                             max.pval=0.05,
                             min.logFC=0.5
)
scsr.no

lri.no.macro.epi <- getParacrines(scsr.no,"macrophages_Tum","epithelial_Tum")
lri.no.epi.macro <- getParacrines(scsr.no,"epithelial_Tum","macrophages_Tum")

mac.epi <- unique(paste(lri.no.macro.epi$L,lri.no.macro.epi$R,sep="¦"))
epi.mac <- unique(paste(lri.no.epi.macro$L,lri.no.epi.macro$R,sep="¦"))
intersect(mac.epi,epi.mac)

lri.no.mono.epi <- getParacrines(scsr.no,"monocytes_Tum","epithelial_Tum")
lri.no.epi.mono <- getParacrines(scsr.no,"epithelial_Tum","monocytes_Tum")
lri.no.neutro.epi <- getParacrines(scsr.no,"neutrophils_Tum","epithelial_Tum")
lri.no.epi.neutro <- getParacrines(scsr.no,"epithelial_Tum","neutrophils_Tum")
