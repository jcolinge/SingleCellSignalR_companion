######################################################################
# Example script to process scProt-MS data and perform LRI inference #
# with SingleCellSignalR v2                                          #
######################################################################

library(data.table)
library(BulkSignalR)
library(SingleCellSignalR)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(eulerr)
library(matrixTests)
library(foreach)

# reset to 2019 databases
lrdb2019 <- fread("data/LRdb_122019.txt",data.table=FALSE)
resetLRdb(lrdb2019,switch=TRUE)
gobp2019 <- fread("data/gobp_2019.txt",data.table=FALSE)
resetPathways(gobp2019,resourceName="GO-BP")
react2019 <- fread("data/reactome_2019.txt",data.table=FALSE)
resetPathways(react2019,resourceName="Reactome")
net2019 <- fread("data/Network_2019.sif",data.table=FALSE)
resetNetwork(net2019)

# uniprot
sprot <- fread("data/sp-name-human-2023-03-30.txt",data.table=F)
genes <- sapply(strsplit(sprot$`Gene Names`," "),function(x) x[1])
sp2gene <- setNames(genes,sprot$Entry)

# Reactome for annotating network nodes
reactome <- fread("data/reactome-human-2023-03-30.txt",data.table=F)

# read scProt-MS data
# downloaded from https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing provided by Slavov lab (https://scp.slavovlab.net/Specht_et_al_2019)
specht.norm <- fread("specht-data/Proteins-processed.csv",data.table=FALSE)
sp.norm <- data.matrix(specht.norm[,-1])
rownames(sp.norm) <- sp2gene[specht.norm[[1]]]
sp.norm <- sp.norm[!is.na(rownames(sp.norm)),]
sp.norm<-sp.norm[,-ncol(sp.norm)]

#read cell annotation
# downloaded from https://drive.google.com/file/d/16vf6rjIsk-oK9naAH6BQnCFrlWnYtJsS/view?usp=sharing provided by Slavov lab (https://scp.slavovlab.net/Specht_et_al_2019)
annot <- fread("specht-data/Cells.csv",data.table=F)
cell.type <- setNames(unlist(annot[2,-1]),unlist(annot[1,-1]))

length(intersect(lrdb2019$ligand,rownames(sp.norm)))
length(intersect(lrdb2019$receptor,rownames(sp.norm)))


# Control heatmaps
# -------------------------------------------------------------------------

spn <- sp.norm
thres <- 2
spn[spn>thres] <- thres
spn[spn< -thres] <- -thres
hist(spn,freq=F,breaks=30)

d.genes <- dist(spn)
h.genes <- hclust(d.genes,method="ward.D")
d.cells <- dist(t(spn))
h.cells <- hclust(d.cells,method="ward.D")

color.scale <- colorRamp2(breaks=c(-2,-1,0,1,2), colors=c("royalblue3", "royalblue3", "white", "orange", "orange"))
ha <- HeatmapAnnotation(
  type=cell.type[colnames(spn)],
  col=list(
    type=setNames(c("darkolivegreen2","coral1"),c("sc_m0","sc_u"))
  )
)
pdf("hm-all-cells-specht.pdf",width=3,height=1.5,pointsize = 8,useDingbats = F)
Heatmap(spn,col=color.scale,use_raster=T,raster_device="png",raster_quality=8,
        cluster_rows=h.genes,cluster_columns=h.cells,#column_split=3,column_gap=unit(1,"mm"),
        show_row_dend=F,
        show_row_names=F,show_column_names=F,bottom_annotation=ha)
dev.off()

m.genes <- c("AIMP2","NCL","BAZ1B","VIM","ANXA2","SOGA1","ITGB2","ELANE")
selected <- spn[m.genes,]
pdf("hm-selected-genes-specht.pdf",width=3,height=1.5,pointsize=8,useDingbats=F)
Heatmap(selected,col=color.scale,use_raster=T,raster_device="png",raster_quality=8,
        cluster_columns=h.cells,#column_split=3,column_gap=unit(1,"mm"),
        show_row_dend=F,row_names_gp=gpar(fontsize=7),
        show_row_names=T,show_column_names=F,bottom_annotation=ha)
dev.off()

ind <- order(sp.norm["VIM",])
sel.vim <- selected[,ind]
ha.vim <- HeatmapAnnotation(
  type=cell.type[colnames(sel.vim)],
  col=list(
    type=setNames(c("darkolivegreen2","coral1"),c("sc_m0","sc_u"))
  )
)
pdf("hm-selected-genes-VIM-specht.pdf",width=3,height=1.1,pointsize=8,useDingbats=F)
Heatmap(sel.vim,col=color.scale,use_raster=T,raster_device="png",raster_quality=8,
        cluster_columns=F,#column_split=3,column_gap=unit(1,"mm"),
        show_row_dend=F,row_names_gp=gpar(fontsize=7),
        show_row_names=T,show_column_names=F,bottom_annotation=ha.vim)
dev.off()


# macrophage autocrine LRIs
# -------------------------------------------------------------------------

# VIM sorting to define mono / macro phenotypes
ind <- order(sp.norm["VIM",])
bottom <- ind[1:trunc(ncol(sp.norm)*0.2)] # mono
top <- ind[-(1:trunc(ncol(sp.norm)*0.8))] # macro

# tailored differential analysis function
getDiffMacro <- function(obj,pop){
  
  if (pop != "macro")
    stop("macro versus mono is the only implemented comparison here")
  
  # define the comparison
  p <- populations(obj)
  macro <- which(p=="macro")
  mono <- which(p=="mono")

  # comparison
  mat <- ncounts(bsrdmComp(obj))
  diff.macro <- row_wilcoxon_twosample(mat[,macro],mat[,mono])$pvalue
  FC <- apply(mat,1,function(x) median(x[macro])-median(x[mono]))
  expr <- apply(mat,1,function(x) median(x[macro]))
  macro.table <- data.frame(pval=diff.macro,logFC=FC,expr=expr)
  rownames(macro.table) <- rownames(mat)
  
  # return the BSRClusterComp object
  BSRClusterComp(bsrdmComp(obj),macro,mono,macro.table)
  
} # getDiffMacro


# build SCSR object
populations <- rep("cells",ncol(sp.norm))
populations[top] <- "macro"
populations[bottom] <- "mono"
minFC <- 0.01
scsr <- SCSRNet(counts=sp.norm, populations=populations,
                normalize=FALSE, method="SCoPE2",
                min.count=0, min.LR.found=10, log.transformed=FALSE
)
scsr <- performInferences(scsr,
                          verbose=TRUE,
                          selected.populations="macro",
                          funDiffExpr=getDiffMacro,
                          max.pval=0.01,
                          min.logFC=minFC,
                          max.pw.size=1200,
                          min.pw.size=4,
                          min.positive=4,
                          min.t.logFC=minFC
)
scsr

# LRI tables
bsrinf.comp.macro <- getAutocrines(scsr,"macro")
macro.inter <- LRinter(bsrinf.comp.macro)
write.table(macro.inter,file="paper/macro-autocrine-proteomics.txt",sep="\t",quote=F,row.names=F)
macro.inter.full <- macro.inter
pvals <- tgPval(bsrinf.comp.macro)
good <- lapply(pvals, function(x) x<0.01)
for (i in 1:length(pvals))
  pvals[[i]] <- pvals[[i]][good[[i]]]
targets <- tgGenes(bsrinf.comp.macro)
for (i in 1:length(pvals))
  targets[[i]] <- targets[[i]][good[[i]]]
corr <- tgCorr(bsrinf.comp.macro)
for (i in 1:length(pvals))
  corr[[i]] <- corr[[i]][good[[i]]]
logFC <- tgLogFC(bsrinf.comp.macro)
for (i in 1:length(pvals))
  logFC[[i]] <- logFC[[i]][good[[i]]]
macro.inter.full$targets <- sapply(targets,function(x) paste(x,collapse=","))
macro.inter.full$target.pvals <- sapply(pvals,function(x) paste(x,collapse=","))
macro.inter.full$target.logFC <- sapply(logFC,function(x) paste(x,collapse=","))
macro.inter.full$target.corr <- sapply(corr,function(x) paste(x,collapse=","))
write.table(macro.inter.full,file="paper/macro-autocrine-proteomics-full.txt",sep="\t",quote=F,row.names=F)

# L & R basic statistics
length(unique(macro.inter$L)) # 14
length(unique(macro.inter$R)) # 8
dim(unique(macro.inter[,1:2])) # 22

# graph
g.macro <- getLRIntracellNetwork(bsrinf.comp.macro,qval.thres=0.05,min.cor=0,max.pval=0.001,min.logFC=log2(1.05))
g.macro
plot(g.macro)
write_graph(g.macro,file="paper/macro-intracellnet.graphml",format="graphml")

# references for network
vertices <- V(g.macro)$name
all.targets <- unique(unlist(targets))
immune.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Adaptive Immune System","Neutrophil degranulation")])
hemo.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`=="Hemostasis"])
rtk.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`=="Signaling by Receptor Tyrosine Kinases"])
vascu.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`=="Cell surface interactions at the vascular wall"])
ecm.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("ECM proteoglycans","Non-integrin membrane-ECM interactions",
                                                                "Extracellular matrix organization","Degradation of the extracellular matrix")])
muscle.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Muscle contraction","Smooth Muscle Contraction")])
gpcr.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Signaling by GPCR")])
neuro.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Neuronal System")])
axong.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Axon guidance")])

plot(euler(list(Im = immune.genes, He = hemo.genes, RTK = rtk.genes, Va = vascu.genes,
                ECM = ecm.genes, Mu = muscle.genes, Neu = neuro.genes, Axg = axong.genes)))
pdf(file="paper/euler-macro-all-proteomics.pdf",useDingbats=FALSE,width=3,height=3,pointsize=7)
e <- euler(list(Im=immune.genes,He=hemo.genes,RTK=rtk.genes,Va=vascu.genes,ECM=ecm.genes,
                Mu=muscle.genes,Neu=neuro.genes,Axg=axong.genes))
plot(e,fills=c("khaki1","palegreen","lightskyblue","limegreen","tan1","plum1","wheat3","lightcyan2"),fontsize=7,quantities=list(fontsize=7))
dev.off()
print(e)
lt = list(Immune = immune.genes, Hemost=hemo.genes,RTK = rtk.genes, Vascul=vascu.genes,ECM = ecm.genes, Muscle = muscle.genes,
          Neuro = neuro.genes, Axon.g = axong.genes)
m1 = make_comb_mat(lt)
m1
pdf(file="paper/upset-macro-all-proteomics.pdf",useDingbats=FALSE,width=8,height=3,pointsize=7)
UpSet(m1)
dev.off()

plot(euler(list(Im = immune.genes, RTK = rtk.genes, ECM = ecm.genes, Mu = muscle.genes, Neu = neuro.genes, Axg = axong.genes)))
pdf(file="paper/euler-macro-proteomics.pdf",useDingbats=FALSE,width=1.5,height=1.5,pointsize=7)
e <- euler(list(Im = immune.genes, RTK = rtk.genes, ECM = ecm.genes, Mu = muscle.genes, Neu = neuro.genes, Axg = axong.genes))
plot(e,fills=c("khaki1","lightskyblue","tan1","plum1","wheat3","lightcyan2"),fontsize=7,quantities=list(fontsize=7))
dev.off()
print(e)
lt = list(Immune = immune.genes, RTK = rtk.genes, ECM = ecm.genes, Muscle = muscle.genes, Neuro = neuro.genes, Axon.g = axong.genes)
m1 = make_comb_mat(lt)
m1
pdf(file="paper/upset-macro-proteomics.pdf",useDingbats=FALSE,width=6,height=3,pointsize=7)
UpSet(m1)
dev.off()

# annotations for Cytoscape rendering
pathway <- setNames(rep("other",length(vertices)),vertices)
pathway[axong.genes] <- "axon.g"
pathway[neuro.genes] <- "neuro"
pathway[muscle.genes] <- "muscle"
pathway[ecm.genes] <- "ecm"
pathway[immune.genes] <- "immun"
pathway[setdiff(rtk.genes,c(immune.genes,ecm.genes))] <- "rtk"
pathway[macro.inter$L] <- "ligand"
pathway[macro.inter$R] <- "receptor"
macro.annotations <- data.frame(name=vertices,
                          detected=vertices%in%rownames(sp.norm),
                          macro.target=vertices%in%all.targets,
                          pathway=pathway)
rownames(macro.annotations) <- macro.annotations$name
table(macro.annotations$pathway)
write.table(macro.annotations,file="paper/proteomics-macro-annotations.txt",sep="\t",quote=F,row.names=F)


# =======================================================================================
# scRNA-seq -----------------------------------------------------------------------------
# =======================================================================================

# data must be first downloaded from GEO website: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142392
# there are two samples (replicates) to download: GSM4226877 & GSM4226878

# read data
rna.1 <- fread("specht-data/GSM4226877_rna_data_Bio_Replicate_1.csv",data.table=F)
rna.2 <- fread("specht-data/GSM4226878_rna_data_Bio_Replicate_2.csv",data.table=F)
rna <- cbind(data.matrix(rna.1[,-1]),data.matrix(rna.2[,-1]))
origin <- c(rep(1,ncol(rna.1)-1),rep(2,ncol(rna.2)-1))
rownames(rna) <- rna.1[[1]]

# basic filtering of low-quality cells and non-expressed genes
ngenes <- colSums(rna>0)
pdf("num-genes-RNA-specht.pdf",width=4,height=3,pointsize=8,useDingbats=F)
plot(sort(ngenes,decreasing=T),type="l",main="",ylab="Number of genes",xlab="Cells")
abline(h=6000,col="red")
abline(h=3500,col="red")
dev.off()
rna <- rna[,ngenes>3500 & ngenes < 6000]
origin <- origin[ngenes>3500 & ngenes < 6000]
dim(rna)
ncells <- rowSums(rna>=2)
sum(ncells>0.15*ncol(rna))
rna <- rna[ncells>0.15*ncol(rna),]
dim(rna)
# [1] 4296 7177
# save(rna,file="specht-combined-scRNA-seq-data.rdta")
# load("specht-combined-scRNA-seq-data.rdta")

# TC normalization
tot <- colSums(rna)
rna.TC <- sweep(rna,2,tot/median(tot),"/")
rna.median <- apply(rna.TC,1,median)
rna.norm <- rna.TC-rna.median
rna.mad <- apply(rna.TC,1,mad)
rna.norm <- rna.norm/rna.mad

# control heatmaps
CV <- rna.mad/rna.median
top.1000 <- head(order(CV,decreasing=T),1000)
rna.spn <- rna.norm[top.1000,]
rna.spn[rna.spn>8] <- 8
rna.spn[rna.spn< -8] <- -8
hist(rna.spn,freq=F,breaks=30)
rna.spn[rna.spn>4] <- 4
rna.spn[rna.spn< -1] <- -1
rna.spn.d <- rna.norm[top.1000,]
d.genes <- dist(rna.spn.d)
h.genes <- hclust(d.genes,method="ward.D")
d.cells <- dist(t(rna.spn.d))
h.cells <- hclust(d.cells,method="ward.D")
ha <- HeatmapAnnotation(
  replicate=origin,
  col=list(
    replicate=setNames(c("lightblue","red"),c("1","2"))
  )
)
color.scale <- colorRamp2(breaks=c(-1,0,4),colors=c("royalblue3","white","orange"))
pdf("scRNA-seq-hm-all-cells-specht.pdf",width=3,height=1.5,pointsize = 8,useDingbats = F)
Heatmap(rna.spn,col=color.scale,use_raster=T,raster_device="png",raster_quality=8,
        cluster_rows=h.genes,cluster_columns=h.cells,#column_split=3,column_gap=unit(1,"mm"),
        show_row_dend=F,show_column_dend=T,
        show_row_names=F,show_column_names=F,bottom_annotation=ha)
dev.off()

m.genes <- c("AIMP2","NCL","BAZ1B","VIM","ANXA2","SOGA1","ITGB2","ELANE")
rna.selected <- rna.norm[intersect(m.genes,rownames(rna.norm)),]
pdf("scRNA-seq-hm-selected-genes-specht.pdf",width=4,height=1.7,pointsize=8,useDingbats=F)
Heatmap(rna.selected,col=color.scale,use_raster=T,raster_device="png",raster_quality=8,
        cluster_columns=h.cells,#column_split=3,column_gap=unit(1,"mm"),
        show_row_dend=F,row_names_gp=gpar(fontsize=7),
        show_row_names=T,show_column_names=F,bottom_annotation=ha)
dev.off()

# macrophage sorting
rna.ind <- order(rna.norm["VIM",]+rna.norm["ITGB2",]-rna.norm["ELANE",]-rna.norm["NCL",])
pdf("control-markers-macro-mono.pdf",width=5,height=5,pointsize=8,useDingbats=F)
par(mfrow=c(2,2))
plot(rna.TC["VIM",rna.ind],pch=".",main="VIM",xlab="Monocytes to Macrophages",ylab="UMI count")
plot(rna.TC["ITGB2",rna.ind],pch=".",main="ITGB2",xlab="Monocytes to Macrophages",ylab="UMI count")
plot(rna.TC["NCL",rna.ind],pch=".",main="NCL",xlab="Monocytes to Macrophages",ylab="UMI count")
plot(rna.TC["ELANE",rna.ind],pch=".",main="ELANE",xlab="Monocytes to Macrophages",ylab="UMI count")
dev.off()
rna.bottom <- rna.ind[1:trunc(ncol(rna.norm)*0.2)] # mono
rna.top <- rna.ind[-(1:trunc(ncol(rna.norm)*0.8))] # macro


# macrophage autocrine LRIs ----------------------------------------------

# quick check for potential LR pairs in scRNA-seq
rna.diff.macro <- apply(rna.TC,1,function(x) wilcox.test(x[rna.top],x[rna.bottom])$p.value)
rna.FC.macro <- apply(rna.TC,1,function(x) log2((median(x[rna.top])+0.5)/(median(x[rna.bottom])+0.5)))
rna.macro.table <- data.frame(pval=rna.diff.macro,logFC=rna.FC.macro)
rownames(rna.macro.table) <- rownames(rna.TC)
minFC <- log2(1.1)
found <- lrdb2019$ligand%in%rownames(rna.macro.table) & lrdb2019$receptor%in%rownames(rna.macro.table)
sum(found)
LRfound <- lrdb2019[found,1:2]
LRstat <- cbind(LRfound,rna.macro.table[LRfound$ligand,],rna.macro.table[LRfound$receptor,])
macro.candidates <- LRstat[LRstat[[3]]<=0.01 & LRstat[[5]]<=0.01 & LRstat[[4]]>=minFC & LRstat[[6]]>=minFC,]

length(intersect(lrdb2019$ligand,rownames(rna.TC)))
length(intersect(lrdb2019$receptor,rownames(rna.TC)))


# analysis with SCSR (similar to proteomics) --------------------

# tailored differential analysis function, same as proteomics but for the logFC
rnaGetDiffMacro <- function(obj,pop,subsample.size=length(top),n.resample = 10){
  
  if (pop != "macro")
    stop("macro versus mono is the only implemented comparison here")
  
  # define the comparison
  p <- populations(obj)
  macro <- which(p=="macro")
  mono <- which(p=="mono")
  
  # comparison
  mat <- ncounts(bsrdmComp(obj))
  # diff.macro <- row_wilcoxon_twosample(mat[,macro],mat[,mono])$pvalue
  d <- foreach(k=seq_len(n.resample),.combine=cbind) %do% {
                          Ap <- sample(macro, subsample.size, replace = TRUE)
                          Bp <- sample(mono, subsample.size, replace = TRUE)
                          row_wilcoxon_twosample(mat[,Ap],mat[,Bp])$pvalue
                        }
  diff <- apply(d, 1, stats::median, na.rm = TRUE)
  # typically caused by all the values being equals
  diff[is.na(diff)] <- 1 
  expr <- apply(mat,1,function(x) median(x[macro]))
  logFC <- (log1p(expr) - log1p( apply(mat,1,function(x) median(x[mono])) )) / log(2)
  macro.table <- data.frame(pval=diff,logFC=logFC,expr=expr)
  rownames(macro.table) <- rownames(mat)
  
  # return the BSRClusterComp object
  BSRClusterComp(bsrdmComp(obj),macro,mono,macro.table)
  
} # rnaGetDiffMacro



rna.populations <- rep("cells",ncol(rna.TC))
rna.populations[rna.top] <- "macro"
rna.populations[rna.bottom] <- "mono"
rna.scsr <- SCSRNet(counts=rna.TC, populations=rna.populations,
                    normalize=FALSE, method="TC",
                    min.count=0, min.LR.found=10, log.transformed=FALSE
)
rna.scsr <- performInferences(rna.scsr,
                              verbose=TRUE,
                              selected.populations="macro",
                              funDiffExpr=rnaGetDiffMacro,
                              max.pval=0.01,
                              min.logFC=minFC,
                              max.pw.size=1200,
                              min.pw.size=4,
                              min.positive=4,
                              min.t.logFC=minFC
)
rna.scsr

# LRI tables
rna.bsrinf.comp.macro <- getAutocrines(rna.scsr,"macro")
rna.macro.inter <- LRinter(rna.bsrinf.comp.macro)
rna.t <- differentialStats(comparison(bsrdmComp(rna.scsr))[["macro_vs_others"]])
write.table(rna.macro.inter,file="paper/macro-autocrine-scRNA-seq.txt",sep="\t",quote=F,row.names=F)
rna.macro.inter.full <- rna.macro.inter
rna.pvals <- tgPval(rna.bsrinf.comp.macro)
good <- lapply(rna.pvals, function(x) x<0.01)
for (i in 1:length(rna.pvals))
  rna.pvals[[i]] <- rna.pvals[[i]][good[[i]]]
rna.targets <- tgGenes(rna.bsrinf.comp.macro)
for (i in 1:length(rna.pvals))
  rna.targets[[i]] <- rna.targets[[i]][good[[i]]]
rna.corr <- tgCorr(rna.bsrinf.comp.macro)
for (i in 1:length(rna.pvals))
  rna.corr[[i]] <- rna.corr[[i]][good[[i]]]
rna.logFC <- tgLogFC(rna.bsrinf.comp.macro)
for (i in 1:length(rna.pvals))
  rna.logFC[[i]] <- rna.logFC[[i]][good[[i]]]
rna.macro.inter.full$targets <- sapply(rna.targets,function(x) paste(x,collapse=","))
rna.macro.inter.full$target.pvals <- sapply(rna.pvals,function(x) paste(x,collapse=","))
rna.macro.inter.full$target.logFC <- sapply(rna.logFC,function(x) paste(x,collapse=","))
rna.macro.inter.full$target.corr <- sapply(rna.corr,function(x) paste(x,collapse=","))
write.table(rna.macro.inter.full,file="paper/macro-autocrine-scRNA-seq-full.txt",sep="\t",quote=F,row.names=F)

length(unique(rna.macro.inter$L)) # 29
length(unique(rna.macro.inter$R)) # 17
dim(unique(rna.macro.inter[,1:2])) # 42

# graph construction
g.rna.macro <- getLRIntracellNetwork(rna.bsrinf.comp.macro,qval.thres=0.05,min.cor=0,max.pval=0.001,min.logFC=minFC)
g.rna.macro
plot(g.rna.macro)
write_graph(g.rna.macro,file="paper/scRNA-seq-macro-intracellnet.graphml",format="graphml")

rna.vertices <- V(g.rna.macro)$name
rna.all.targets <- unique(unlist(rna.targets))
rna.hemo.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`=="Hemostasis"])
rna.vascu.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`=="Cell surface interactions at the vascular wall"])
rna.ecm.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("ECM proteoglycans","Non-integrin membrane-ECM interactions",
                                                                                   "Extracellular matrix organization","Degradation of the extracellular matrix")])
rna.integr.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Integrin cell surface interactions")])
rna.il.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Signaling by Interleukins","Interleukin-4 and Interleukin-13 signaling",
                                                                                  "Interleukin-2 family signaling","Interleukin-3, Interleukin-5 and GM-CSF signaling")])
rna.tgfb.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("TGF-beta receptor signaling activates SMADs",
                                                                                    "TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)",
                                                                                    "Signaling by TGF-beta family members","Signaling by TGF-beta Receptor Complex")])
rna.gpcr.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("G alpha (i) signalling events")])
rna.immune.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Adaptive Immune System","Neutrophil degranulation")])
rna.nr.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Signaling by Nuclear Receptors")])
rna.mapk.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Oncogenic MAPK signaling")])
rna.muscle.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Muscle contraction","Smooth Muscle Contraction")])
rna.neuro.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Neuronal System")])
rna.axong.genes <- intersect(rna.vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Axon guidance")])

lt = list(Hemo=rna.hemo.genes,Vasc=rna.vascu.genes,ECM=rna.ecm.genes,Itgr=rna.integr.genes,IL=rna.il.genes,
          Imm=immune.genes,TGFb=rna.tgfb.genes,GPCR=rna.gpcr.genes,NR=rna.nr.genes,MAPK=rna.mapk.genes)
m1 = make_comb_mat(lt)
pdf(file="paper/upset-macro-all-transcriptomics.pdf",useDingbats=FALSE,width=12,height=4.5,pointsize=7)
UpSet(m1)
dev.off()

lt = list(ECM=rna.ecm.genes,IL=rna.il.genes,Imm=immune.genes,MAPK=rna.mapk.genes,
          GPCR=rna.gpcr.genes,NR=rna.nr.genes)
m1 = make_comb_mat(lt)
pdf(file="paper/upset-macro-transcriptomics.pdf",useDingbats=FALSE,width=7,height=4,pointsize=7)
UpSet(m1)
dev.off()

# annotations for Cytoscape rendering
pathway <- setNames(rep("other",length(rna.vertices)),rna.vertices)
pathway[rna.mapk.genes] <- "mapk"
pathway[rna.immune.genes] <- "immun"
pathway[rna.nr.genes] <- "nr"
pathway[rna.gpcr.genes] <- "gpcr"
pathway[rna.il.genes] <- "il"
pathway[rna.ecm.genes] <- "em"
pathway[rna.macro.inter$L] <- "ligand"
pathway[rna.macro.inter$R] <- "receptor"
macro.annotations <- data.frame(name=rna.vertices,
                                detected=rna.vertices%in%rownames(rna.TC),
                                macro.target=rna.vertices%in%rna.all.targets,
                                pathway=pathway)
rownames(macro.annotations) <- macro.annotations$name
table(macro.annotations$pathway)
write.table(macro.annotations,file="paper/transcriptomics-macro-annotations.txt",sep="\t",quote=F,row.names=F)


# merged networks ------------------------------------------

el  = unique(rbind(as_data_frame(g.macro)[,1:2],as_data_frame(g.rna.macro)[,1:2]))
keys <- apply(el,1,function(x) paste(sort(x),collapse="||"))
dup <- duplicated(keys)
sum(dup)
GU = graph_from_data_frame(el[!dup,],directed=F)
GU
write_graph(GU,file="paper/union-macro-intracellnet.graphml",format="graphml")
union.vertices <- V(GU)$name
status <- setNames(rep("union",length(union.vertices)),union.vertices)
status[setdiff(vertices,rna.vertices)] <- "prot"
status[setdiff(rna.vertices,vertices)] <- "rna"
union.annotations <- data.frame(name=union.vertices,
                                status=status)
write.table(union.annotations,file="paper/union-macro-annotations.txt",sep="\t",quote=F,row.names=F)

# betweenness study
bet <- betweenness(GU,directed=FALSE,normalized=TRUE)

pdf(file="paper/betweenness.pdf",useDingbats=FALSE,width=1.5,height=1.5,pointsize=7)
boxplot(list(prot=bet[union.annotations$status=="prot"],
             rna=bet[union.annotations$status=="rna"],
             both=bet[union.annotations$status=="union"]),pch=".",lwd=1)
dev.off()
wilcox.test(bet[union.annotations$status=="prot"],bet[union.annotations$status=="union"]) # < 2.2E-16
wilcox.test(bet[union.annotations$status=="rna"],bet[union.annotations$status=="union"]) # < 2.2E-16
wilcox.test(bet[union.annotations$status=="prot"],bet[union.annotations$status=="rna"]) # 0.061


# Comparisons proteomics / transcriptomics ----------------------------------------

# global

proteo.LR <- unique(paste(macro.inter$L,macro.inter$R,sep=" / "))
rna.LR <- unique(paste(rna.macro.inter$L,rna.macro.inter$R,sep=" / "))
pdf(file="paper/euler-LR.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(scRNA=rna.LR,scProt=proteo.LR))
plot(e,fills=c("#B2DFFFFF","#FFA0A0FF"),fontsize=7,quantities=list(fontsize=7))
dev.off()
rna.rownames <- rownames(rna.TC)
prot.rownames <- rownames(sp.norm)
pot.LR <- rbind(lrdb2019[lrdb2019$ligand%in%rna.rownames & lrdb2019$receptor%in%rna.rownames,1:2],
                lrdb2019[lrdb2019$ligand%in%prot.rownames & lrdb2019$receptor%in%prot.rownames,1:2])
pot.LR <- unique(pot.LR)
N <- nrow(pot.LR)
K <- 6+14 # proteomics
k <- 6 # intersection
n <- 36+6 # transcriptomics
phyper(q=k-1,m=K,n=N-K,k=n,lower.tail=F) # 0.34 non-significant

pdf(file="paper/euler-L.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(scRNA=unique(rna.macro.inter$L),scProt=unique(macro.inter$L)))
plot(e,fills=c("#B2DFFFFF","#FFA0A0FF"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-R.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(scRNA=unique(rna.macro.inter$R),scProt=unique(macro.inter$R)))
plot(e,fills=c("#B2DFFFFF","#FFA0A0FF"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-network-nodes.pdf",useDingbats=FALSE,width=1.5,height=1.5,pointsize=7)
e <- euler(list(scRNA=V(g.rna.macro)$name,scProt=V(g.macro)$name))
plot(e,fills=c("#B2DFFFFF","#FFA0A0FF"),fontsize=7,quantities=list(fontsize=7))
dev.off()
N <- length(union(rna.rownames,prot.rownames))
K <- 115+136 # proteomics
k <- 115 # intersection
n <- 242+115 # transcriptomics
phyper(q=k-1,m=K,n=N-K,k=n,lower.tail=F) # ==> 1.7E-79 very signif large

pdf(file="paper/euler-network-data.pdf",useDingbats=FALSE,width=1.5,height=1.5,pointsize=7)
e <- euler(list(scRNA=rna.rownames,scProt=prot.rownames))
plot(e,fills=c("#B2DFFFFF","#FFA0A0FF"),fontsize=7,quantities=list(fontsize=7))
dev.off()
N <- length(union(rna.rownames,prot.rownames))
K <- length(prot.rownames) # proteomics
k <- 1237 # intersection
n <- length(rna.rownames) # transcriptomics
phyper(q=k-1,m=K,n=N-K,k=n,lower.tail=F)
phyper(q=k,m=K,n=N-K,k=n) # ==> very significant small

# by pathway

pdf(file="paper/euler-macro-compare-immune.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(Im.r=rna.immune.genes,Imm.p=immune.genes))
plot(e,fills=c("white","khaki1"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-macro-compare-ecm.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(ecm.r=rna.ecm.genes,ecm.p=ecm.genes))
plot(e,fills=c("white","tan1"),fontsize=7,quantities=list(fontsize=7))
dev.off()

il.genes <- intersect(vertices,reactome$`Gene name`[reactome$`Reactome name`%in%c("Signaling by Interleukins","Interleukin-4 and Interleukin-13 signaling",
                                                                                  "Interleukin-2 family signaling","Interleukin-3, Interleukin-5 and GM-CSF signaling")])
pdf(file="paper/euler-macro-compare-IL.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(IL.r=rna.il.genes,IL.p=il.genes))
plot(e,fills=c("white","lightpink"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-macro-compare-gpcr.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(GPCR.r=rna.gpcr.genes,GPCR.p=gpcr.genes))
plot(e,fills=c("white","cyan"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-macro-compare-muscle.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(lip.r=rna.muscle.genes,lip.p=muscle.genes))
plot(e,fills=c("white","plum1"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-macro-compare-axong.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(NR.r=rna.axong.genes,NR.p=axong.genes))
plot(e,fills=c("#FFFFFFFF","lightcyan2"),fontsize=7,quantities=list(fontsize=7))
dev.off()

pdf(file="paper/euler-macro-compare-neuro.pdf",useDingbats=FALSE,width=1,height=1,pointsize=7)
e <- euler(list(pyru.r=rna.neuro.genes,pyru.p=neuro.genes))
plot(e,fills=c("#FFFFFFFF","wheat3"),fontsize=7,quantities=list(fontsize=7))
dev.off()
