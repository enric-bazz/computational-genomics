# test_gene_expression_analysis.R
#
# Runs the three functions defined in gene_expression_analysis.R end to end on
# the provided raw count and annotation files, writing MvAplots.pdf,
# DE_summary.txt and DE_results.txt to an output directory.
#
# This script lives in 02_bulk_RNAseq/scripts/. By default it reads its input
# from the sibling 02_bulk_RNAseq/data/ folder and writes results to
# 02_bulk_RNAseq/output/.
#
# Usage:
#   - From the command line:   Rscript test_gene_expression_analysis.R [input_dir] [output_dir]
#     (both arguments are optional; they default to 02_bulk_RNAseq/data and 02_bulk_RNAseq/output)
#   - From an interactive R/RStudio session:
#       source("test_gene_expression_analysis.R")
#       run_gene_expression_analysis()                                  # use all defaults
#       run_gene_expression_analysis(input_dir = "...", output_dir = "...")

# Locate the directory this script lives in, so input/output paths can default
# to it regardless of the current working directory. Works both when run via
# Rscript (--file= command line argument) and when run via source() (the
# 'ofile' variable set by source() on the call stack)
get_script_dir <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg))))
    }

    frames <- sys.frames()
    ofiles <- Filter(Negate(is.null), lapply(frames, function(f) f$ofile))
    if (length(ofiles) > 0) {
        return(dirname(normalizePath(as.character(ofiles[[length(ofiles)]]))))
    }

    # Fallback if the script's own path can't be determined
    getwd()
}

# Capture the script directory while the source()/Rscript frame is still on the
# call stack
THIS_SCRIPT_DIR <- get_script_dir()

DEFAULT_DATA_DIR <- normalizePath(file.path(THIS_SCRIPT_DIR, "..", "data"), mustWork = FALSE)
DEFAULT_OUTPUT_DIR <- normalizePath(file.path(THIS_SCRIPT_DIR, "..", "output"), mustWork = FALSE)


run_gene_expression_analysis <- function(input_dir = DEFAULT_DATA_DIR,
                                         output_dir = DEFAULT_OUTPUT_DIR,
                                         alpha = 0.05) {

    # Load functions
    source(file.path(THIS_SCRIPT_DIR, "gene_expression_analysis.R"))

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Load input data
    rawdata <- read.table(
        file.path(input_dir, "raw_count.txt"),
        header = TRUE,
        sep = "\t",
        row.names = 1
    )

    annot <- read.table(
        file.path(input_dir, "gene_annot.txt"),
        header = TRUE,
        sep = "\t",
        row.names = "Symbol",
        quote = ""
    )

    # Sample groups are encoded in the sample names (Group1_*, Group2_*, CTRL*)
    groups <- colnames(rawdata)


    ### 1. MvA plots

    mva_file <- file.path(output_dir, "MvAplots.pdf")
    MvAplot(rawdata, mva_file)
    cat("MvA plots saved to:", mva_file, "\n")


    ### 2. TMM normalization

    TMM_results <- TMMnorm(rawdata, annot)

    cat("\nTMM normalization completed.\n")
    cat("\nScaling factors:\n")
    print(TMM_results$ScalingFactors)
    cat("\nNormalized expression matrix dimensions:", dim(TMM_results$Normalized_exprData), "\n")


    ### 3. Differential expression analysis

    DE_results <- DEbyEdgeR(rawdata, groups, alpha = alpha)

    summary_file <- file.path(output_dir, "DE_summary.txt")
    write.table(
        data.frame(statistic = names(DE_results$DEtestOutcomes),
                   value = unname(DE_results$DEtestOutcomes)),
        file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE
    )

    results_file <- file.path(output_dir, "DE_results.txt")
    write.table(
        DE_results$DEtestResults,
        file = results_file, sep = "\t", quote = FALSE, col.names = NA
    )

    cat("\nDifferential expression summary (alpha =", alpha, "):\n")
    print(DE_results$DEtestOutcomes)
    cat("\nWrote:\n -", summary_file, "\n -", results_file, "\n")

    invisible(list(TMM = TMM_results, DE = DE_results))
}


# When run as a standalone script (e.g. via Rscript), execute the analysis using
# the command line arguments for input/output directories, if given
if (!interactive() && sys.nframe() == 0) {
    cli_args <- commandArgs(trailingOnly = TRUE)
    input_dir <- if (length(cli_args) >= 1) cli_args[1] else DEFAULT_DATA_DIR
    output_dir <- if (length(cli_args) >= 2) cli_args[2] else DEFAULT_OUTPUT_DIR
    run_gene_expression_analysis(input_dir = input_dir, output_dir = output_dir)
}
