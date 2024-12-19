#sce_count.RData -> raw count data from a single-cell RNA experiment of mouse arteries, stored as SingleCellExperiment object (SCE)
#ref.RData -> contains an SCE object w\ reference dataset to be used for cell type annotation
#LRdb.RData -> contains data.frame w\ ligand-receptor pairs

## FIRST FUNCTION ##

#input:
 #sce: SCE object supposed to have same format of sce_count 
 #n_HVG: # of top highly variable genes to select
 #percVar: % of variance to be explained by the principal components PCs

#output:
 #updated and processed sce object, which includes the following:
  #logcounts: assay w\ log-normalized gene expression data
  #results of PCA and t-SNE in their respective slots 
  #metadata field HVGs listing the top HVGs

#steps to be performed: 
 #filter out low-expression genes ( = genes w\ expression value of 0 in all cells)
 #log-normalize the raw count data using scran package
 #identify the top n_HVG highly variable genes using variance modelling
 #determine and print the number of PCs required to explain at least a percVar % of the total variance; adjust PCA to include only the required number of PCs (you may want to adjust the ncomponents argument in runPCA() )
 # perform t-SNE using top HVG subset

# Required libraries 
library(SingleCellExperiment)
library(scran)
library(scater)
library(Rtsne)
library(scuttle)
library(SingleR)
library(igraph)

process_SCE <- function(sce, n_HVG = 1000, percVar = 70){
  # Filter out low-expression genes ( = genes with expression value of 0 in all cells)
  sce <- sce[rowSums(counts(sce) > 0) > 0, ]
  
  # Log-normalize raw count data using scran
  # Compute size factors
  sce <- computeLibraryFactors(sce) 
  # Log-normalize
  sce <- logNormCounts(sce)         
  
  # Identify top n_HVG highly variable genes using variance modelling
  var_out <- modelGeneVar(sce)
  hvg <- getTopHVGs(var_out, n = n_HVG, var.field='bio')
  # Subset the SCE object to the top HVGs
  sce <- sce[hvg, ]    
  # Store the HVGs in the metadata
  metadata(sce)$HVGs <- hvg       
  
  # Perform PCA and determine number of PCs explaining at least percVar % of the variance
  # Perform PCA with default 10 components
  sce <- runPCA(sce, subset_row = hvg, ncomponents = 10) 
  # Get variance explained by PC
  pca_var <- attr(reducedDim(sce, "PCA"), "percentVar") 
  # Cumulative variance
  cumulative_var <- cumsum(pca_var)   
  # Find the number of PCs needed
  n_PCs <- which(cumulative_var >= percVar)[1]         
  cat("Number of PCs required to explain at least", percVar, "% variance:", n_PCs, "\n")
  # Re-run PCA with the required PCs
  sce <- runPCA(sce, subset_row = hvg, ncomponents = n_PCs) 
  
  # Perform t-SNE using top HVG subset
  # Run t-SNE
  tsne_out <- Rtsne(reducedDim(sce, "PCA")[, 1:n_PCs], perplexity = 30) 
  # Store t-SNE results
  reducedDim(sce, "TSNE") <- tsne_out$Y                               
  
  # Return the processed SCE object
  return(sce) 
  
}

#to test#
load('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 3/sce_count.RData')
processed_sce <- process_sce(sce, n_HVG = #value, percVar = #value) #un-comment and insert values

# the results
print(processed_sce)
plotPCA(processed_sce)
plotTSNE(processed_sce)


## SECOND FUNCTION##

#input:
 #sce
 #ref: an SCE object to use as reference dataset for annotation ( cell ontology class in colData)
 #pdffilename: name of the PDF file where the tSNE plots will be saved
 
#output:
 #updated and processed sce object, which includes 'cell_annotation' column in colData slot

#steps to be performed:
 #compute cell clustering with Louvain algorithm with k=10 nearest neighborhood
 #perform cell type annotation using SingleR at cluster level (resultinf from the cell clustering), using ref as reference dataset ( use assay.type.ref = 'counts')
 #generate and save 2 t-SNE plots to the specified pdffilename, coloring the plots based on clusters (plot 1) and cell type annotation (plot 2)

annotation <- function(sce, reference_data, output_pdf) {
  
  # Build a Shared Nearest Neighbor (SNN) graph using Louvain algorithm
  # k = 10 nearest neighbors based on PCA-reduced dimensions
  snn_graph <- buildSNNGraph(sce, k = 10, use.dimred = "PCA")
  
  # Louvain clustering on the SNN graph
  louvain_clusters <- cluster_louvain(snn_graph)
  
  # Add Louvain clustering results to SingleCellExperiment  
  sce$Louvain_Cluster <- factor(membership(louvain_clusters))
  
  # Perform cell type annotation using SingleR
  # Log normalize the reference dataset
  reference_data <- logNormCounts(reference_data)
  
  # SingleR cell type annotation
  cell_type_predictions <- SingleR(test = sce, 
                                   ref = reference_data, 
                                   labels = reference_data$cell_ontology_class, 
                                   clusters = sce$Louvain_Cluster, 
                                   de.method = "wilcox", 
                                   assay.type.ref = "counts", 
                                   assay.type.test = "counts")
  
  # Add the predicted cell annotations to the SingleCellExperiment 
  sce$Cell_Type_Annotation <- factor(cell_type_predictions$labels[sce$Louvain_Cluster], 
                                     levels = unique(cell_type_predictions$labels))
  
  # Generate t-SNE plots for visualization
  # Extract t-SNE coordinates
  tsne_coordinates <- reducedDim(sce, "TSNE")
  tsne_data <- as.data.frame(tsne_coordinates)
  
  # Add Louvain cluster and cell type annotation as metadata for t-SNE data
  tsne_data$Louvain_Cluster <- factor(sce$Louvain_Cluster)
  tsne_data$Cell_Type_Annotation <- factor(sce$Cell_Type_Annotation)
  
  # Create PDF for plots
  pdf(output_pdf)
  
  # Plot t-SNE with Louvain Clusters
  plot(tsne_data$TSNE1, tsne_data$TSNE2, col = tsne_data$Louvain_Cluster, 
       main = "t-SNE: Louvain Clusters", xlab = "TSNE1", ylab = "TSNE2")
  
  # Plot t-SNE with Cell Type Annotations
  plot(tsne_data$TSNE1, tsne_data$TSNE2, col = tsne_data$Cell_Type_Annotation, 
       main = "t-SNE: Cell Type Annotations", xlab = "TSNE1", ylab = "TSNE2")
  
  # Close PDF 
  dev.off()
  
  # Confirmation message
  cat("PDF file created successfully at", output_pdf, "\n")
  
  # Return the updated SingleCellExperiment 
  return(sce)
}


#to test#
load('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 3/ref.RData')  # Load reference data outside the function
processed_sce <- annotation(processed_sce, ref, "tsne_plots.pdf")


## THIRD FUNCTION##

# Perform a simple cell-cell communication analysis ( = analyze potential communication between different
# cell typed in a single-cell-RNA-seq dataset identified throught the previous function)
# Objective: compute an interaction score (LRscore) for each ligand-receptor pair (LR) across all possible combinations of sender (expressing ligand) and receiver (expressing receptor) cell types (cell type pairs)

# The code must follow this structure:
# for each unique cell type pair (e.g 'cell type A' as sender and ' cell type B' as receiver)
 # for each ligand-receptor pair LR (e.d 'gene 1' as ligand, 'gene 2 ' as receptor)
  # calculate the mean log-normalized expression of the ligand gene across all sender cells (e.g. means expression of 'gene 1 ' in 'cell type A')
  # calculate the mean log-normalized expression of the receptor gene across all receiver cells (e.g mean expression of 'gene 2' in 'cell B type')
  # compute LRscore as the product of the ligand mean expression and receptor mean expression

#input: 
 #sce
 #LRdb: dataframe with the known LR pairs ( 2 columns 'ligand' and 'receptor' with respectively the ligand and receptor gene names)

#output: 
 #dataframe where each row corresponds to a ligand-receptor pair for a specific sender-receiver cell type pair: teh columns include
 # c('ligand','receptor','sender','receiver','LRscore') sorted in descending order by LRscore to prioritize the strongest potential interactions

CCCanalysis <- function(sce,LRdb){
  # Extract required data
  # Log-normalized expression
  logcounts_data <- logcounts(sce)         
  # Cell type annotations
  cell_annotations <- colData(sce)$cell_annotation    
  
  # Initialize results dataframe
  results <- data.frame(
    ligand = character(),
    receptor = character(),
    sender = character(),
    receiver = character(),
    LRscore = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Unique cell types
  cell_types <- unique(cell_annotations)
  
  # Loop through all cell type pairs
  for (sender in cell_types) {
    for (receiver in cell_types) {
      # Get sender and receiver cells
      sender_cells <- which(cell_annotations == sender)
      receiver_cells <- which(cell_annotations == receiver)
      
      # Loop through each ligand-receptor pair in LRdb
      for (i in seq_len(nrow(LRdb))) {
        ligand <- LRdb$ligand[i]
        receptor <- LRdb$receptor[i]
        
        # Check if the ligand and receptor are present in the dataset
        if (ligand %in% rownames(logcounts_data) && receptor %in% rownames(logcounts_data)) {
          # Compute mean expression for ligand in sender cells
          ligand_mean <- mean(logcounts_data[ligand, sender_cells], na.rm = TRUE)
          
          # Compute mean expression for receptor in receiver cells
          receptor_mean <- mean(logcounts_data[receptor, receiver_cells], na.rm = TRUE)
          
          # Compute LRscore
          LRscore <- ligand_mean * receptor_mean
          
          # Append result to the dataframe
          results <- rbind(results, data.frame(
            ligand = ligand,
            receptor = receptor,
            sender = sender,
            receiver = receiver,
            LRscore = LRscore,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Sort results by LRscore in descending order
  results <- results[order(-results$LRscore), ]
  
  return(results)
}

# to test#
load("/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 3/sce_count.RData")
load("/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 3/ref.RData")
load("/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 3/LRdb.RData")

processed_sce <- process_sce(sce, n_HVG = 1000, percVar = 70)
annotated_sce <- annotation(processed_sce, ref, "tsne_plots.pdf")
interaction_results <- CCCanalysis(annotated_sce, LRdb)
print(interaction_results)
