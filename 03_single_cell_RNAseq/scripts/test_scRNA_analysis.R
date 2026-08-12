### Test script for scRNA-seq analysis functions

main <- function(input_dir, output_dir, script_dir) {

    library(SingleCellExperiment)
    
    ## Create output directory if necessary
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    ## Load functions
    
    source(file.path(script_dir, "scRNA_analysis.R"))
    
    ## Load data
    
    load(file.path(input_dir, "sce_count.RData"))
    load(file.path(input_dir, "ref.RData"))
    load(file.path(input_dir, "LRdb.RData"))
    
    ## Parameters
    
    n_HVG <- 1000
    percVar <- 70
    
    ## Process scRNA-seq data
    
    sce <- process_sce(
        sce,
        n_HVG = n_HVG,
        percVar = percVar
    )
    
    ## Annotate cell populations
    
    sce <- annotation(
        sce,
        ref,
        file.path(output_dir, "tSNE_annotation.pdf")
    )
    
    ## t-SNE comparison
    
    perplexities <- c(50, 100, 400)
    
    pdf(
        file = file.path(output_dir, "tSNE_perplexity_comparison.pdf"),
        width = 15,
        height = 5
    )
    
    par(mfrow = c(1, 3))
    
    for (perplexity in perplexities) {
        
        tsne_name <- paste0("TSNE_", perplexity)
        
        set.seed(123)
        
        sce <- scater::runTSNE(
            sce,
            subset_row = sce$HVGs,
            dimred = "PCA",
            assay.type = "logcounts",
            perplexity = perplexity,
            name = tsne_name
        )
        
        print(
            scater::plotReducedDim(
                sce,
                dimred = tsne_name,
                colour_by = "Louvain"
            ) +
                ggplot2::ggtitle(paste0("Perplexity = ", perplexity))
        )
    }
    
    dev.off()
    
    ## Cell-cell communication analysis
    
    CCC_results <- CCCanalysis(
        sce,
        LRdb
    )
    
    write.table(
        CCC_results,
        file = file.path(output_dir, "CCC_results.txt"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    
    cat("Analysis completed successfully.\n")
}


## Entry point

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop(
        "Usage: Rscript scripts/test_scRNA_analysis.R <input_dir> <output_dir>"
    )
}

script_path <- commandArgs()[grep(
    "--file=",
    commandArgs()
)]

script_dir <- dirname(
    normalizePath(
        sub("--file=", "", script_path)
    )
)

main(
    input_dir = normalizePath(args[1]),
    output_dir = normalizePath(args[2], mustWork = FALSE),
    script_dir = script_dir
)