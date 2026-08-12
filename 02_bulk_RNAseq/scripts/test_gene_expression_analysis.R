### Test script for gene_expression_analysis.R

run_tests <- function(input_dir = "data", output_dir = "output") {
    
    # Load functions
    source(file.path("scripts", "gene_expression_analysis.R"))
    
    # Create output directory if it does not already exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Input files
    rawdata_file <- file.path(input_dir, "raw_count.txt")
    annot_file <- file.path(input_dir, "gene_annot.txt")
    
    # Load input data
    rawdata <- read.table(
        rawdata_file,
        header = TRUE,
        sep = "\t",
        row.names = 1
    )
    
    annot <- read.table(
        annot_file,
        header = TRUE,
        sep = "\t",
        row.names = "Symbol",
        quote = ""
    )
    
    # Sample groups are encoded in the sample names
    groups <- colnames(rawdata)
    
    
    ### 1. MvA plots
    
    MvAplot(
        rawdata,
        file.path(output_dir, "MvAplots.pdf")
    )
    
    cat("MvA plots saved to:", file.path(output_dir, "MvAplots.pdf"), "\n")
    
    
    ### 2. TMM normalization
    
    TMM_results <- TMMnorm(
        rawdata,
        annot
    )
    
    normalized_data <- TMM_results$Normalized_exprData
    scaling_factors <- TMM_results$ScalingFactors
    
    cat("\nTMM normalization completed.\n")
    cat("\nScaling factors:\n")
    print(scaling_factors)
    
    cat("\nNormalized expression matrix dimensions:\n")
    print(dim(normalized_data))
    
    
    ### 3. Differential expression analysis
    
    DE_results <- DEbyEdgeR(
        rawdata,
        groups
    )
    
    DE_output_file <- file.path(output_dir, "DE_results.txt")
    
    sink(DE_output_file)
    
    cat("Differential expression summary\n")
    cat("===============================\n\n")
    print(DE_results$DEtestOutcomes)
    
    cat("\n\nPer-gene differential expression results\n")
    cat("=======================================\n\n")
    print(DE_results$DEtestResults)
    
    sink()
    
    cat("DE results saved to:", DE_output_file, "\n")
}


### Entry point

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1) args[1] else "data"
output_dir <- if (length(args) >= 2) args[2] else "output"

run_tests(input_dir, output_dir)