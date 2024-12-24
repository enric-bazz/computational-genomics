#Set working direceptortory to filepath, data files are in the same direceptortory
dir <- rstudioapi::getActiveDocumentContext()
setwd(dirname(dir$path))


## Load data
library(SingleCellExperiment)
load("sce_count.RData")
load('ref.RData')
n_HVG <- 1000 #number of top expressed genes to consider
percVar <- 70 #percentage of explained variance for PCA

# Filter zeros and normalization
allzero_genes <- rowSums(counts(sce))==0 
sce <- sce[which(!allzero_genes),] 
sce <- scuttle::logNormCounts(sce) #perform normalization + store normalized counts in logcounts assay

# Top highly variable genes
var_model <- scran::modelGeneVar(sce)
top_HVGs <- scran::getTopHVGs(var_model, n=n_HVG)
metadata(sce)$HVGs <- top_HVGs

# PCA
n_PCs <- 93
sce <- scater::runPCA(sce, subset_row=top_HVGs, ncomponents=n_PCs, name='PCA') #PCA of high variance genes, fixed #PCs


# Clustering
SNN_graph <- scran::buildSNNGraph(sce, k=10, use.dimred='PCA') #10 neighbors
clust <- igraph::cluster_louvain(SNN_graph)
sce$Louvain <- factor(igraph::membership(clust))

# SingleR annotation
logcounts(ref) <- counts(ref) #copy assay
pred <- SingleR::SingleR(test=sce, ref=ref, labels=ref$cell_ontology_class, 
                                   clusters = sce$Louvain, #annotation on clusters 
                                   de.method = "wilcox",assay.type.ref = "counts", assay.type.test = "counts")
sce$cell_annotation <- factor(pred$labels[sce$Louvain], levels = unique(pred$labels))


# tSNE tuning
set.seed(123)
sce <- scater::runTSNE(sce, subset_row = top_HVGs, dimred="PCA",
                       assay.type='logcounts', perplexity=50, name="tSNE_50")
set.seed(123)
sce <- scater::runTSNE(sce, subset_row = top_HVGs, dimred="PCA",
                       assay.type='logcounts', perplexity=100, name="tSNE_100")
set.seed(123)
sce <- scater::runTSNE(sce, subset_row = top_HVGs, dimred="PCA",
                       assay.type='logcounts', perplexity=400, name="tSNE_400")


png(file='tSNE.png', width = 1800, height = 600, res=120)
par(mfrow = c(1, 3))
# create t-SNE plots 
tsne_50 <- as.data.frame(reducedDim(sce, "tSNE_50"))
tsne_50$cell_annotation <- factor(sce$cell_annotation)
plot(tsne_50$TSNE1, tsne_50$TSNE2, col = tsne_50$cell_annotation, 
     main = "Perplexity = 50", xlab = "tSNE1", ylab = "tSNE2")
legend("bottomright",legend=levels(tsne_50$cell_annotation), pch=1, col=unique(tsne_50$cell_annotation),bty='n')

tsne_100 <- as.data.frame(reducedDim(sce, "tSNE_100"))
tsne_100$cell_annotation <- factor(sce$cell_annotation)
plot(tsne_100$TSNE1, tsne_100$TSNE2, col = tsne_100$cell_annotation, 
     main = "Perplexity = 100", xlab = "tSNE1", ylab = "tSNE2")
legend("topleft",legend=levels(tsne_100$cell_annotation), pch=1, col=unique(tsne_100$cell_annotation),bty='n')


tsne_400 <- as.data.frame(reducedDim(sce, "tSNE_400"))
tsne_400$cell_annotation <- factor(sce$cell_annotation)
plot(tsne_400$TSNE1, tsne_400$TSNE2, col = tsne_400$cell_annotation, 
     main = "Perplexity = 400", xlab = "tSNE1", ylab = "tSNE2")
legend("topleft",legend=levels(tsne_400$cell_annotation), pch=1, col=unique(tsne_400$cell_annotation),bty='n')
dev.off()