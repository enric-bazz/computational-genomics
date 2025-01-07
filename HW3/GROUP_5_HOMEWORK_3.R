### Function 1: process_sce()
## INPUT
# - sce     : a SingleCellExperiment object that is supposed to have the same format of the sce_count data provided in stem

# - n_HVG   : a numeric value specifying the number of top highly variable genes to select

# - percVar : a numeric value indicating the percentage of variance to be explained by the principal components PCs


## OUTPUT
# Gives as output the update and processed sce object, which includes the following
# • A logcounts assay with log-normalized gene expression data
# • Results of PCA and t-SNE in their respective slots (i.e. «PCA» and «TSNE»)
# • A metadata field «HVGs» listing the top HVGs

# CODE
process_sce <- function(sce, n_HVG, percVar) {
  
  ## Load SCE package
  library(SingleCellExperiment)
  
  ## Filter out low expression genes (zero count across all cells)
  allzero_genes <- rowSums(counts(sce))==0 #boolean vector identifying all-zero-counts genes as TRUE
  sce <- sce[which(!allzero_genes),] #excludes low expression genes (conditional indexing)
  
  ## Log-normalization of raw counts
  sce <- scuttle::logNormCounts(sce) #performs normalization and stores normalized counts in the new logcounts assay
  # NOTE: In the Lab, normalization follows the scaling factor computation based on cell pooling
  # cl <- quickCluster(sce)
  # sce <- computeLibraryFactors(sce, clusters=cl) 
  # In the scuttle package version 1.16.0 the argument "clusters=cl" is ignored and scaling factors computed internally with logNormCounts()
  # coincide with the ones coming from computeLibraryFactors(), so these steps are neglected.
  
  ## Highly Variable Genes identification
  # 1. Variance modelling
  var_model <- scran::modelGeneVar(sce) #computes variance of log-normalized expression values: technical + biological
  # 2. Identify top=n_HVG high variance genes
  top_HVGs <- scran::getTopHVGs(var_model, 
                                n = n_HVG, # number of genes to select
                                var.field = "bio") #variance component used
  metadata(sce)$HVGs <- top_HVGs #stores HVGs in metadata field of SCE object
  
  ## PCA on highly variable genes subset
  sce <- scater::runPCA(sce, 
                        subset_row = top_HVGs, #genes subset
                        ncomponents = n_HVG/10, #number of PCs computed
                        name = "PCA")  #name of dimensionality reduction in SCE object
  # NOTE: the number of components is set high enough to capture most variance, the default value of 50 is not sufficient.
  # For values of percVar larger than the trial one (70) an even larger number of components might be necessary.
  
  ## Number of PCs based on cumulative explained variance
  explainedVar_PCA <- cumsum(attr(reducedDim(sce, "PCA"), which = "percentVar")) #computes cumulative variance explained by PCs
  n_PCs <- min(which(explainedVar_PCA>=percVar)) #returns minimum index of logical TRUE = minimum number of PCs to explain at least percVar % of the total
  
  # Display number of selected PCs
  print(paste0("The number of Principal Components necessary to explain at least ", toString(percVar), "% of the total variance is: ",toString(n_PCs)))
  
  ## Adjusted PCA restricted to first n_PCs principal components
  sce <- scater::runPCA(sce, 
                        subset_row = top_HVGs, #genes subset
                        ncomponents = n_PCs, #number of PCs computed
                        name = "PCA") #name of dimensionality reduction in SCE object
  
  ## tSNE on top HVGs (default to 2-dimensional space)
  sce <- scater::runTSNE(sce, 
                         subset_row = top_HVGs, #genes subset
                         assay.type = 'logcounts', #assay used
                         perplexity = 100,  #perplexity hyper-parameter value
                         dimred="PCA", #using PCA improves speed and reduces random noise (source: scater documentation on Bioconductor)
                         name = "TSNE") #name of dimensionality reduction in SCE object
  # NOTE: Perplexity value tuning 
  # The optimal perplexity value is selected via trial and error taking into consideration these two factors:
  # - the intent of capturing the global structure more than the local one with the purpose of visualizing cluster-level annotation
  # - the size of the dataset, especially in terms of the number of cells (approx. 10.000 columns/observations)
  
  ## OUTPUT
  return(sce)
}






### Function 2: annotation()
## INPUT
# - sce         : a SingleCellExperiment object having the same format of the output of process_sce() function

# - ref         : a SingleCellExperiment object to use as reference dataset for annotation (cell_ontology_class in colData)

# - pdffilename :   the name of the PDF file where the tSNE plots will be saved


## OUTPUT
# Gives as output the update and processed sce object, which includes “cell_annotation” column in colData slot

# CODE
annotation <- function(sce, ref, pdffilename) {
  
  ## Load SCE package
  library(SingleCellExperiment)
  
  ## Cell clustering with Louvain algorithm
  # 1. Build Shared-Nearest-Neighbor graph of test data
  SNN_graph <- scran::buildSNNGraph(sce, 
                                    k = 10, #number of neighbors
                                    use.dimred = "PCA" #dimensionally reduced space to build the graph on
  ) 
  # 2. Perform Louvain clustering
  clust <- igraph::cluster_louvain(SNN_graph)
  # 3. Saves clusters field in SCE object
  sce$Louvain <- factor(igraph::membership(clust))
  
  ## Cell type annotation at cluster level: SingleR package
  logcounts(ref) <- counts(ref) #add "proxy" logcounts assay to reference  
  pred <- SingleR::SingleR(test = sce, #test set
                           ref = ref, #reference set
                           labels = ref$cell_ontology_class, #labels in reference for annotation
                           clusters = sce$Louvain, #annotation at clusters level
                           de.method = "wilcox", #differential expression test
                           assay.type.ref = "counts", #reference assay
                           assay.type.test = "counts" #test assay
  )
  # NOTE: If the 'logcounts' assay is not added on the reference data (ref), SingleR() returns an error (assay.type.ref="counts" argument is ignored). 
  # To ensure that the reference assay used in the SingleR annotation contains the non-logarithmic counts, instead of adding the actual logcounts assay to the ref object,
  # the raw counts assay is copied under the "logcounts" name.
  
  ## Add annotation column on SCE object colData slot
  sce$cell_annotation <- factor(pred$labels[sce$Louvain], levels = unique(pred$labels))
  
  ## tSNE plots in PDF file
  pdf(file = pdffilename, onefile = TRUE)
  print(scater::plotTSNE(sce, colour_by = "Louvain")) #Louvain clusters highlighted
  print(scater::plotTSNE(sce, colour_by = "cell_annotation")) #SingleR cell types highlighted
  dev.off()
  # Confirmation message
  cat("PDF file created successfully at", pdffilename, "\n")
  
  ## OUTPUT
  return(sce)
  
} 




### Function 3: CCCanalysis()
## INPUT
# - sce  : a SingleCellExperiment object having the same format of the output of annotation() function

# - LRdb : a dataframe with the known LR pairs (two columns “ligand” and “receptor” with respectively the ligand and receptor gene names)


## OUTPUT
# Gives as output a dataframe, where each row corresponds to a ligand-receptor pair for a specific sender-receiver
# cell type pair. The columns include c(“ligand”, “receptor”, “sender”, “receiver”, “LRscore”),
# sorted in descending order by LRscore to prioritize the strongest potential interactions

# CODE
CCCanalysis <- function(sce, LRdb) {
  
  ## Load SCE package
  library(SingleCellExperiment)
  
  # NOTE: the SCE object in input (sce) is assumed to have the exact format of the annotation() output
  
  ## Retrieve all cell types annotated in the SCE object 
  cell_types <- levels(sce$cell_annotation)
  
  ## Build a matrix of all possible sender-receiver pairs: matrix with all partial permutations with repetitions of the cell types
  cell_type_pairs <- gtools::permutations(n = length(cell_types), #source vector size = number of cell types
                                          r = 2, #target vector size
                                          v = cell_types, #source vector = cell types character vector
                                          repeats.allowed = TRUE #include self-communication pairs
  ) 
  # Assign sender and receiver names to the columns
  colnames(cell_type_pairs) <- c('sender','receiver')
  
  ## Exclude ligand/receptor genes not present in the dataset (sce object)
  ligand_present <- LRdb$ligand %in% rowData(sce)$SYMBOL #boolean of non-missing ligand genes
  receptor_present <- LRdb$receptor %in% rowData(sce)$SYMBOL #boolean of non-missing receptor genes
  LRdb <- LRdb[ligand_present | receptor_present,] #exclude LRdb rows with a missing ligand or receptor gene
  
  ## Initialize output data frame columns as single vectors
  ligand <- receptor <- sender <- receiver <- character() #empty character vectors
  LRscore_vec <- numeric() #empty numeric vector
  
  ## LR score computation cycle
  for (s_r in 1:nrow(cell_type_pairs)) {    #cycles through each sender-receiver cell types pair (s_r index)
    
    for (l_r in 1:nrow(LRdb)) {   #cycles through every ligand-receptor genes pair (l_r index)
      
      # Append current sender/receiver cell types to corresponding vectors
      sender <- c(sender, cell_type_pairs[s_r, 'sender'])
      receiver <- c(receiver, cell_type_pairs[s_r, 'receiver'])
      #values are retrieved from the current row and the two respective columns of the 'cell_type_pairs' matrix
      
      # Append current ligand/receptor genes to corresponding vectors
      ligand <- c(ligand, LRdb$ligand[l_r])
      receptor <- c(receptor, LRdb$receptor[l_r])
      #values are retrieved from the current row and the two respective fields of the 'LRdb' data frame
      
      ## Compute mean log-normalized expressions:
      # Select the rows of SCE object by gene name and the columns by annotated cell type 
      # NOTE: logcounts(sce) is a matrix-like object so conditional indexing can be exploited
      
      # Sender cell - ligand gene
      mean_LNE_sender <- mean(logcounts(sce)[LRdb$ligand[l_r], #select row corresponding to the current ligand gene
                                             sce$cell_annotation==cell_type_pairs[s_r,'sender']]) #select columns corresponding to cells annotated as the sender type
      
      # Receiver cell - receptor gene
      mean_LNE_reciever <- mean(logcounts(sce)[LRdb$receptor[l_r], #select row corresponding to the current receptor gene
                                               sce$cell_annotation==cell_type_pairs[s_r,'receiver']]) #select columns corresponding to cells annotated as the receiver type
      
      ## Compute LRScore as the product of mean log-expressions
      LRScore <- mean_LNE_sender*mean_LNE_reciever
      
      # Append current score to vector
      LRscore_vec <- c(LRscore_vec,LRScore)
      
    } # ligand-receptor for
    
  } # sender-receiver for
  
  ## Build output data frame
  CCC_results <- data.frame(
    ligand = ligand,
    receptor = receptor,
    sender = sender,
    receiver = receiver,
    LRScore = LRscore_vec)
  
  ## Sort data frame based on descending LRScore
  CCC_results <- CCC_results[order(CCC_results$LRScore, decreasing = TRUE),]
  # Re-assign row names: they get "mixed" with the sorting procedure
  rownames(CCC_results) <- 1:nrow(CCC_results)
  
  ## OUTPUT
  return(CCC_results)
}