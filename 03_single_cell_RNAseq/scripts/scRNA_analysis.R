### Function 1: process_sce()
#
# INPUT
# - sce     : SingleCellExperiment object containing raw count data
# - n_HVG   : number of top highly variable genes to select
# - percVar : percentage of variance to be explained by the retained PCs
#
# OUTPUT
# Returns the processed SingleCellExperiment object containing:
# - logcounts assay with log-normalized expression data
# - PCA representation stored in the "PCA" slot
# - t-SNE representation stored in the "TSNE" slot
# - "HVGs" metadata field containing the selected highly variable genes

process_sce <- function(sce, n_HVG = 1000, percVar = 70) {
    
    ## Filter out genes with zero expression across all cells
    
    allzero_genes <- rowSums(counts(sce)) == 0
    sce <- sce[!allzero_genes, ]
    
    ## Log-normalization of raw counts
    
    cl <- scran::quickCluster(sce)
    sce <- scran::computeSumFactors(sce, clusters = cl)
    sce <- scuttle::logNormCounts(sce)
    
    ## Highly Variable Genes identification
    
    var_model <- scran::modelGeneVar(sce)
    
    top_HVGs <- scran::getTopHVGs(
        var_model,
        n = n_HVG,
        var.field = "bio"
    )
    
    metadata(sce)$HVGs <- top_HVGs
    
    ## Initial PCA on highly variable genes
    
    n_PCs_initial <- min(100, nrow(sce), ncol(sce) - 1)
    
    sce <- scater::runPCA(
        sce,
        subset_row = top_HVGs,
        ncomponents = n_PCs_initial,
        name = "PCA"
    )
    
    ## Number of PCs based on cumulative explained variance
    
    explainedVar_PCA <- cumsum(
        attr(reducedDim(sce, "PCA"), which = "percentVar")
    )
    
    PCs_required <- which(explainedVar_PCA >= percVar)
    
    if (length(PCs_required) == 0) {
        stop(
            paste0(
                "The first ", n_PCs_initial,
                " PCs do not explain ", percVar,
                "% of the variance."
            )
        )
    }
    
    n_PCs <- min(PCs_required)
    
    print(
        paste0(
            "The number of Principal Components necessary to explain at least ",
            percVar,
            "% of the total variance is: ",
            n_PCs
        )
    )
    
    ## PCA restricted to the required number of components
    
    sce <- scater::runPCA(
        sce,
        subset_row = top_HVGs,
        ncomponents = n_PCs,
        name = "PCA"
    )
    
    ## t-SNE
    
    sce <- scater::runTSNE(
        sce,
        subset_row = top_HVGs,
        assay.type = "logcounts",
        perplexity = 100,
        dimred = "PCA",
        name = "TSNE"
    )
    
    return(sce)
}




### Function 2: annotation()
#
# INPUT
# - sce         : processed SingleCellExperiment object
# - ref         : reference SingleCellExperiment containing cell-type labels
#                in cell_ontology_class
# - pdffilename : name of the PDF file for the t-SNE plots
#
# OUTPUT
# Returns the annotated SingleCellExperiment object containing:
# - "Louvain" cluster assignments
# - "cell_annotation" cell-type annotations
# - t-SNE plots saved to the specified PDF file

annotation <- function(sce, ref, pdffilename) {
    
    ## Cell clustering with Louvain algorithm
    
    SNN_graph <- scran::buildSNNGraph(
        sce,
        k = 10,
        use.dimred = "PCA"
    )
    
    set.seed(100)
    clust <- igraph::cluster_louvain(SNN_graph)
    
    sce$Louvain <- factor(igraph::membership(clust))
    
    ## Cell-type annotation with SingleR
    
    # Required because the reference object does not contain a logcounts assay; SingleR expects it internally even when assay.type.ref="counts"
    logcounts(ref) <- counts(ref)
    
    pred <- SingleR::SingleR(
        test = sce,
        ref = ref,
        labels = ref$cell_ontology_class,
        clusters = sce$Louvain,
        de.method = "wilcox",
        assay.type.ref = "counts",
        assay.type.test = "counts"
    )
    
    cluster_labels <- setNames(pred$labels, rownames(pred))
    
    sce$cell_annotation <- factor(
        cluster_labels[as.character(sce$Louvain)]
    )
    
    ## t-SNE visualizations
    
    pdf(file = pdffilename, onefile = TRUE)
    
    print(
        scater::plotTSNE(
            sce,
            colour_by = "Louvain"
        )
    )
    
    print(
        scater::plotTSNE(
            sce,
            colour_by = "cell_annotation"
        )
    )
    
    dev.off()
    
    cat("PDF file created successfully at", pdffilename, "\n")
    
    return(sce)
}




### Function 3: CCCanalysis()
#
# INPUT
# - sce  : annotated SingleCellExperiment object from annotation()
# - LRdb : ligand-receptor database
#
# OUTPUT
# Returns a data frame containing:
# - ligand
# - receptor
# - sender
# - receiver
# - LRscore
#
# Rows are sorted by decreasing LRscore.

CCCanalysis <- function(sce, LRdb) {
    
    ## Retrieve annotated cell types
    
    cell_types <- levels(sce$cell_annotation)
    
    ## Generate all sender-receiver cell-type pairs
    
    cell_type_pairs <- gtools::permutations(
        n = length(cell_types),
        r = 2,
        v = cell_types,
        repeats.allowed = TRUE
    )
    
    colnames(cell_type_pairs) <- c("sender", "receiver")
    
    ## Keep ligand-receptor pairs whose genes are present in the dataset
    
    ligand_present <- LRdb$ligand %in% rowData(sce)$SYMBOL
    receptor_present <- LRdb$receptor %in% rowData(sce)$SYMBOL
    
    LRdb <- LRdb[ligand_present & receptor_present, ]
    
    ## Initialize output vectors
    
    ligand <- receptor <- sender <- receiver <- character()
    LRscore_vec <- numeric()
    
    ## Compute interaction scores
    
    for (s_r in 1:nrow(cell_type_pairs)) {
        
        for (l_r in 1:nrow(LRdb)) {
            
            current_sender <- cell_type_pairs[s_r, "sender"]
            current_receiver <- cell_type_pairs[s_r, "receiver"]
            
            ligand <- c(ligand, LRdb$ligand[l_r])
            receptor <- c(receptor, LRdb$receptor[l_r])
            sender <- c(sender, current_sender)
            receiver <- c(receiver, current_receiver)
            
            mean_LNE_sender <- mean(
                logcounts(sce)[
                    LRdb$ligand[l_r],
                    sce$cell_annotation == current_sender
                ]
            )
            
            mean_LNE_receiver <- mean(
                logcounts(sce)[
                    LRdb$receptor[l_r],
                    sce$cell_annotation == current_receiver
                ]
            )
            
            LRscore <- mean_LNE_sender * mean_LNE_receiver
            
            LRscore_vec <- c(LRscore_vec, LRscore)
        }
    }
    
    ## Build and sort output
    
    CCC_results <- data.frame(
        ligand = ligand,
        receptor = receptor,
        sender = sender,
        receiver = receiver,
        LRscore = LRscore_vec
    )
    
    CCC_results <- CCC_results[
        order(CCC_results$LRscore, decreasing = TRUE),
    ]
    
    rownames(CCC_results) <- seq_len(nrow(CCC_results))
    
    return(CCC_results)
}