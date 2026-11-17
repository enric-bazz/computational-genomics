# test_scRNA_analysis.R
#
# Runs the three functions defined in scRNA_analysis.R end to end on the
# provided single-cell data, writing tSNE_annotation.pdf,
# tSNE_perplexity_comparison.pdf and CCC_results.txt to an output directory.
#
# This script lives in 03_single_cell_RNAseq/scripts/. By default it reads its
# input from the sibling 03_single_cell_RNAseq/data/ folder and writes results
# to 03_single_cell_RNAseq/output/.
#
# Usage:
#   - From the command line:   Rscript test_scRNA_analysis.R [input_dir] [output_dir]
#     (both arguments are optional; they default to 03_single_cell_RNAseq/data
#     and 03_single_cell_RNAseq/output)
#   - From an interactive R/RStudio session:
#       source("test_scRNA_analysis.R")
#       run_scRNA_analysis()                                  # use all defaults
#       run_scRNA_analysis(input_dir = "...", output_dir = "...")

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


run_scRNA_analysis <- function(input_dir = DEFAULT_DATA_DIR,
                               output_dir = DEFAULT_OUTPUT_DIR,
                               n_HVG = 1000,
                               percVar = 70,
                               perplexities = c(50, 100, 400),
                               seed = 123) {

    library(SingleCellExperiment)

    # Load functions
    source(file.path(THIS_SCRIPT_DIR, "scRNA_analysis.R"))

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    ## Load data (objects: sce, ref, LRdb)

    load(file.path(input_dir, "sce_count.RData"))
    load(file.path(input_dir, "ref.RData"))
    load(file.path(input_dir, "LRdb.RData"))

    ## Process scRNA-seq data
    # quickCluster() and runTSNE() are stochastic: fix the seed so that the
    # processed object and the t-SNE in tSNE_annotation.pdf are reproducible

    set.seed(seed)
    sce <- process_sce(sce, n_HVG = n_HVG, percVar = percVar)

    ## Annotate cell populations

    sce <- annotation(sce, ref, file.path(output_dir, "tSNE_annotation.pdf"))

    ## t-SNE perplexity comparison
    # Same PCA and Louvain clusters, only the perplexity changes. When dimred is
    # given, runTSNE() works on the PCs (computed on the HVGs) and ignores
    # subset_row, which is passed only for consistency with process_sce()

    tsne_plots <- lapply(perplexities, function(perplexity) {

        tsne_name <- paste0("TSNE_", perplexity)

        set.seed(seed)
        sce_tsne <- scater::runTSNE(
            sce,
            subset_row = metadata(sce)$HVGs,
            dimred = "PCA",
            assay.type = "logcounts",
            perplexity = perplexity,
            name = tsne_name
        )

        scater::plotReducedDim(sce_tsne, dimred = tsne_name, colour_by = "Louvain") +
            ggplot2::ggtitle(paste0("Perplexity = ", perplexity))
    })

    # ggplot objects ignore par(mfrow): place them side by side on a single
    # page with a grid layout
    pdf(file = file.path(output_dir, "tSNE_perplexity_comparison.pdf"),
        width = 5 * length(perplexities), height = 5)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(perplexities))))
    for (i in seq_along(tsne_plots)) {
        print(tsne_plots[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
    }
    dev.off()

    ## Cell-cell communication analysis

    CCC_results <- CCCanalysis(sce, LRdb)

    write.table(
        CCC_results,
        file = file.path(output_dir, "CCC_results.txt"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat("\nCell types:", levels(sce$cell_annotation), sep = "\n - ")
    cat("\nTop 5 ligand-receptor interactions:\n")
    print(head(CCC_results, 5))
    cat("\nResults written to", output_dir, "\n")

    invisible(sce)
}


# When run as a standalone script (e.g. via Rscript), execute the analysis using
# the command line arguments for input/output directories, if given
if (!interactive() && sys.nframe() == 0) {
    cli_args <- commandArgs(trailingOnly = TRUE)
    input_dir <- if (length(cli_args) >= 1) cli_args[1] else DEFAULT_DATA_DIR
    output_dir <- if (length(cli_args) >= 2) cli_args[2] else DEFAULT_OUTPUT_DIR
    run_scRNA_analysis(input_dir = input_dir, output_dir = output_dir)
}
