
### FIRST FUNCTION

## INPUT
# sce: A SingleCellExperimentobject that is supposed to have the same format of the sce_count data provided in stem
# n_HVG: a numeric value specifying the number of top highly variable genes to select
# percVar: a numeric value indicating the percentage of variance to be explained by the principal components PCs


process_sce <- function(sce, n_HVG = 1000, percVar = 70) {
  
  library(SingleCellExperiment)
  library(scuttle)
  library(scran)
  library(scater)
  
  # filter out low-expression genes (genes with zero counts across all cells)
  allzero <- rowMeans(counts(sce) == 0) == 1  # Identifies genes with zero counts across all cells
  # table(allzero)  
  sce <- sce[which(!allzero), ]  # Removes genes with zero expression in all cells
  
  # log-normalize the raw counts using scuttle
  sce <- logNormCounts(sce) # Stores normalized counts in 'logcounts' assay
  
  # identify the top n_HVG highly variable genes (HVG) using variance modelling
  var_model <- modelGeneVar(sce) 
  
  top_HVGs <- getTopHVGs(var_model, n = n_HVG, var.field = "bio") # top HVGs based on biological variance
  metadata(sce)$HVGs <- top_HVGs 
  
  # principal components using the HVG subset
  sce <- runPCA(sce, subset_row = top_HVGs, ncomponents = n_HVG/10, name = "PCA") 
  
  # calculate cumulative variance explained by PCs
  explainedVar_PCA <- cumsum(attr(reducedDim(sce, "PCA"), which = "percentVar"))
  
  # find the number of PCs required to explain at least percVar% of variance
  n_PCs <- which(explainedVar_PCA >= percVar)[1]  # Find the first PC index where the variance meets the threshold
  
  cat("The number of Principal Components necessary to explain at least", percVar, "% of variance is:", n_PCs, "\n")
  
  # PCA again with the adjusted number of PCs
  sce <- runPCA(sce, subset_row = top_HVGs, ncomponents = n_PCs, name = "PCA")
  
  # t-SNE on top HVGs
  sce <- runTSNE(sce, subset_row = top_HVGs, assay.type = "logcounts", perplexity = 100, name = "TSNE")
  
  return(sce)
}


## OUTPUT
# Gives as output the update and processed sce-object, which includes the following
# logcountsassay with log-normalized gene expression data
# Results of PCA and t-SNE in their respective slots (i.e. «PCA» and «TSNE»)
# metadata field «HVGs» listing the top HVGs


load('sce_count.RData')
processed_sce <- process_sce(sce, n_HVG = 1000, percVar = 70)

# the results
print(processed_sce)
plotPCA(processed_sce)
plotTSNE(processed_sce)



### SECOND FUNCTION

## INPUT
# sce: a SingleCellExperimentobject having the same format of the output of process_sce() function
# ref:a SingleCellExperimentobject to use as reference dataset for annotation (cell_ontology class in, colData) 
# pdffilename: the name of the PDF file where the tSNE plots will be saved

annotation <- function(sce, ref, pdffilename) {
  library(SingleCellExperiment)
  library(SingleCellExperiment)
  library(SingleR)
  library(igraph)
  
  # cell clustering with Louvain algorithm with k=10 nearest neighborhood
  SNN_graph <- buildSNNGraph(sce, k = 10, use.dimred = "PCA")
  # add Louvain clustering to colData
  clust <- cluster_louvain(SNN_graph)
  sce$Louvain <- factor(membership(clust))
  
  # cell type annotation with SingleR
  ref <- logNormCounts(ref)
  pred <- SingleR(test = sce, ref = ref, labels = ref$cell_ontology_class, 
                  clusters = sce$Louvain, de.method = "wilcox", 
                  assay.type.ref = "counts", assay.type.test = "counts")
  sce$cell_annotation <- factor(pred$labels[sce$Louvain], levels = unique(pred$labels))
  
  # create t-SNE plots 
  tsne_data <- reducedDim(sce, "TSNE")
  tsne_df <- as.data.frame(tsne_data)
  
  # add metadata
  tsne_df$Louvain <- factor(sce$Louvain)
  tsne_df$cell_annotation <- factor(sce$cell_annotation)
  
  pdf(pdffilename)
  
  # Louvain Clusters
  plot(tsne_df$TSNE1, tsne_df$TSNE2, col = tsne_df$Louvain, 
       main = "t-SNE: Louvain Clusters", xlab = "TSNE1", ylab = "TSNE2")
  
  
  # cell Type Annotations
  plot(tsne_df$TSNE1, tsne_df$TSNE2, col = tsne_df$cell_annotation, 
       main = "t-SNE: Cell Type Annotations", xlab = "TSNE1", ylab = "TSNE2")
  
  dev.off()
  cat("PDF file created successfully\n")

  return(sce)
}


## OUTPUT 
# the update and processed sceobject, which includes “cell_annotation”column in colDataslot

load("ref.RData")  # Load reference data outside the function
processed_sce <- annotation(processed_sce, ref, "tsne_plots.pdf")



### THIRD FUNCTION

## INPUT
# sce: a SingleCellExperiment object having the same format of the output of annotation() function
# LRdb: a dataframe with the known LR pairs (two columns “ligand” and “receptor” with respectively the ligand and receptor gene names)

CCCanalysis <- function(sce, LRdb) {
  library(SingleCellExperiment)
  
  # init an empty data frame to store the results
  results <- data.frame(ligand = character(),
                        receptor = character(),
                        sender = character(),
                        receiver = character(),
                        LRscore = numeric(),
                        stringsAsFactors = FALSE)
  
  # cell type annotations
  cell_types <- unique(sce$cell_annotation)
  
  # iterate over each unique cell type pair (sender-receiver)
  for (sender in cell_types) {
    for (receiver in cell_types) {
      # subset cells by sender and receiver types
      sender_cells <- sce[, sce$cell_annotation == sender]
      receiver_cells <- sce[, sce$cell_annotation == receiver]
      
      # iterate over each ligand-receptor pair
      for (i in seq_len(nrow(LRdb))) {
        ligand <- LRdb$ligand[i]
        receptor <- LRdb$receptor[i]
        
        # calculation of mean log-normalized expression of the ligand gene across all sender cells
        if (ligand %in% rownames(sender_cells)) {
          ligand_expr <- rowMeans(logcounts(sender_cells)[ligand, , drop = FALSE])
        } else {
          ligand_expr <- 0
        }
        
        # of the receptor gene across all receiver cells
        if (receptor %in% rownames(receiver_cells)) {
          receptor_expr <- rowMeans(logcounts(receiver_cells)[receptor, , drop = FALSE])
        } else {
          receptor_expr <- 0
        }
        
        # LRscore as the product of the ligand mean expression and receptor mean expression
        LRscore <- ligand_expr * receptor_expr
        
      
        results <- rbind(results, data.frame(ligand = ligand,
                                             receptor = receptor,
                                             sender = sender,
                                             receiver = receiver,
                                             LRscore = LRscore,
                                             stringsAsFactors = FALSE))
      }
    }
  }
  
  # sorting in descending order by LRscore
  results <- results[order(results$LRscore, decreasing = TRUE), ]
  
  return(results)
}

## OUTPUT
# a dataframe, where each row corresponds to a ligand-receptor pair for a specific sender-receiver cell type pair.
# The columns include c(“ligand”, “receptor”, “sender”, “receiver”, “LRscore”), sorted in descending order by LRscore to prioritize the strongest potential interactions


load("sce_count.RData")
load("ref.RData")
load("LRdb.RData")

processed_sce <- process_sce(sce, n_HVG = 1000, percVar = 70)
annotated_sce <- annotation(processed_sce, ref, "tsne_plots.pdf")
interaction_results <- CCCanalysis(annotated_sce, LRdb)
print(interaction_results)
