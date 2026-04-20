library(SingleCellSignalR)
library(BulkSignalR)

# get the data from the library
data(example_dataset, package='SingleCellSignalR')

# trivial normalization (log2 transform)
mat <- log1p(data.matrix(example_dataset[,-1]))/log(2)
rownames(mat) <- example_dataset[[1]]

# define cell populations
rme <- rowMeans(mat)
mmat <- mat[rme>0.05,]
d <- dist(t(mmat))
h <- hclust(d, method="ward.D")
pop <- paste0("pop_", cutree(h, 5))
table(pop)


# ==============================================================================
# analysis without pathway below the receptors is very
# quick but simplistic and it does not return BSRInferenceComp
# objects, i.e., follow up operations are limited

scsr.no <- SCSRNoNet(counts=mat,
                          populations=pop,
                          normalize=FALSE,
                          method="log2-transform",
                          log.transformed=TRUE)
scsr.no

# infer L-R interactions
scsr.no <- performInferences(scsr.no,
                             verbose=TRUE,
                             min.logFC=0.001,
                             max.pval=0.05)
scsr.no

# retrieves different sets of interactions

# Note 1: L-R interactions are directly returned as a data frame
# in the case of SCSRNoNet objects

# Note 2: SingleCellSignalR version 1 original LR-score is computed
# as well and returned in a column named 'LR.score' (in addition to
# the independent P-value obtained by BulkSignalR statistical model).
# Check performInferences parameters to see how thresholds can be
# imposed on the various score and filtering done using both or just
# one score. Below, no filtering is applied to LR-score (min.LR.score=0)

inter.a.pop.1 <- getAutocrines(scsr.no, "pop_1")
inter.a.pop.5 <- getAutocrines(scsr.no, "pop_5")
inter.p.pop.1.2 <- getParacrines(scsr.no, "pop_1", "pop_2")
inter.p.pop.1.3 <- getParacrines(scsr.no, "pop_1", "pop_3")

cellNetBubblePlot(scsr.no)
cellNetHeatmap(scsr.no, bar.num=TRUE, bar.plot.height=unit(2, "cm"))

# version without resampling (not recommended as results depend on population sizes)
scsr.no.no <- SCSRNoNet(counts=mat,
                        populations=pop,
                        normalize=FALSE,
                        method="log2-transform",
                        log.transformed=TRUE)
scsr.no.no

# infer L-R interactions
scsr.no.no <- performInferences(scsr.no.no,
                                subsample.size=-1,
                                verbose=TRUE,
                                min.logFC=log2(1.1),
                                max.pval=0.05)
scsr.no.no
cellNetHeatmap(scsr.no.no, bar.num=TRUE, bar.plot.height=unit(2, "cm"))



# ==============================================================================
# Deeper analysis including pathways below the receptors and getting
# results as BSRInferenceComp objects. Takes more time though.

# We also examplify the use of a parallel computing backend. Here,
# we use doParallel to be compatible with Windows, which cannot fork.
library(doParallel)
registerDoParallel(min(c(8, detectCores()-1)))

scsr.with <- SCSRNet(counts=mat,
                     populations=pop,
                     normalize=FALSE,
                     method="log2-transform",
                     log.transformed=TRUE)
scsr.with

# infer L-R interactions, takes longer
scsr.with <- performInferences(scsr.with,
                               verbose=TRUE,
                               min.logFC=0.001,
                               max.pval=0.05,
                               max.pw.size=200,
                               min.t.logFC=0.005)
scsr.with

# retrieves different sets of interactions

# Note 1: L-R interactions are returned as BSRInferenceComp objects,
# i.e., LRinter must be used if we want to get them as a data.frame

# Note 2: SingleCellSignalR version 1 original LR-score is also computed
# here

# Note 3: since results are BRSInferenceComp objects, all the follow up
# tools to work with such objects in the BulkSignalR library are available

winter.a.pop.1 <- LRinter(getAutocrines(scsr.with, "pop_1"))

# P-values tend to get very small, impose support by target gene regulation
sum(winter.a.pop.1$qval<1e-3 & winter.a.pop.1$rank.pval<0.05)

# we can also apply BulkSignalR reduction operations
winter.a.pop.1.red <- LRinter(reduceToBestPathway(getAutocrines(scsr.with, "pop_1")))
sum(winter.a.pop.1.red$qval<1e-3 & winter.a.pop.1.red$rank.pval<0.05)

winter.a.pop.5 <- LRinter(getAutocrines(scsr.with, "pop_5"))
sum(winter.a.pop.5$qval<1e-3 & winter.a.pop.5$rank.pval<0.05)

winter.p.pop.1.2 <- LRinter(getParacrines(scsr.with, "pop_1", "pop_2"))
winter.p.pop.1.2.red <- LRinter(reduceToBestPathway(getParacrines(scsr.with, "pop_1", "pop_2")))

# plots with counts of LR interactions
cellNetBubblePlot(scsr.with)
cellNetHeatmap(scsr.no, bar.num=TRUE, bar.plot.height=unit(2, "cm"))
cellNetHeatmap(scsr.with, bar.num=TRUE, bar.plot.height=unit(2, "cm"))

# plots focusing on a given pathway
LRdb <- getResource("LRdb")
reactome <- getResource("Reactome")
genes <- intersect(LRdb$receptor, reactome$`Gene name`[reactome$`Reactome ID`=="R-HSA-449147"]) # Signaling by Interleukins
cellNetBubblePlot(scsr.with, genes.to.count=genes, use.proportions=TRUE)
cellNetBubblePlot(scsr.with, genes.to.count=genes, use.proportions=TRUE,
                  low.color="gray", high.color="royalblue3")
cellNetHeatmap(scsr.with, genes.to.count=genes, use.proportions=TRUE)
cellNetHeatmap(scsr.with, genes.to.count=genes, high.color="tomato")
