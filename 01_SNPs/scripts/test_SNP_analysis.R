# test_SNP_analysis.R
#
# Runs the SNP_association_test() pipeline defined in SNP_analysis.R end to
# end on a genotype data file, writing association_results.txt and
# test_boxplot.pdf to an output directory.
#
# This script lives in 01_SNPs/scripts/. By default it reads its input from
# the sibling 01_SNPs/data/ folder and writes results to 01_SNPs/output/.
#
# Usage:
#   - From the command line:   Rscript test_SNP_analysis.R [input_dir] [output_dir]
#     (both arguments are optional; they default to 01_SNPs/data and 01_SNPs/output)
#   - From an interactive R/RStudio session:
#       source("test_SNP_analysis.R")
#       run_snp_analysis()                                  # use all defaults
#       run_snp_analysis(input_dir = "...", output_dir = "...")

#Locate the directory this script lives in, so input/output paths can default
#to it regardless of the current working directory. Works both when run via
#Rscript (using the --file= command line argument) and when run via source()
#in an interactive session (using the 'ofile' set by source() on the call stack)
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
  
  #Fallback if the script's own path can't be determined (e.g. pasted directly
  #into a console)
  getwd()
}

#Capture this script's own directory now, while the call stack still includes
#the source()/Rscript frame that get_script_dir() inspects. Calling
#get_script_dir() later (e.g. as a lazily-evaluated function default) would no
#longer see that frame once sourcing has finished
THIS_SCRIPT_DIR <- get_script_dir()

#This script lives in 01_SNPs/scripts/; default input/output paths point to
#the sibling data/ and output/ folders (output/ is created if missing)
DEFAULT_DATA_DIR <- normalizePath(file.path(THIS_SCRIPT_DIR, "..", "data"), mustWork = FALSE)
DEFAULT_OUTPUT_DIR <- normalizePath(file.path(THIS_SCRIPT_DIR, "..", "output"), mustWork = FALSE)

#Main entry point: reads a genotype file from input_dir, runs the full SNP
#association analysis, and writes the results/boxplot to output_dir
run_snp_analysis <- function(input_dir = DEFAULT_DATA_DIR,
                              output_dir = DEFAULT_OUTPUT_DIR,
                              data_file = "SNPdata.txt",
                              indCTRL = 1201:2000,
                              MAFth = 0.01,
                              HWEalpha = 0.01) {
  
  source(file.path(THIS_SCRIPT_DIR, "SNP_analysis.R"))
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  #SNP_association_test() hardcodes the HWE boxplot filename to
  #"test_boxplot.pdf" (as required by the assignment). Temporarily switch the
  #working directory to output_dir so that file is written there, then
  #restore the original working directory no matter how this function exits
  old_wd <- setwd(output_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  cat("Running SNP_association_test() on", file.path(input_dir, data_file), "...\n\n")
  
  result <- SNP_association_test(
    filepath = file.path(input_dir, data_file),
    indCTRL = indCTRL,
    MAFth = MAFth,
    HWEalpha = HWEalpha
  )
  
  results_path <- file.path(output_dir, "association_results.txt")
  write.table(result, file = results_path, sep = "\t", quote = FALSE, col.names = NA)
  
  #Quick summary of the results
  cat("\n--- Results summary ---\n")
  cat("SNPs passing MAF/HWE filters:", nrow(result), "\n")
  cat("SNPs significant at qval < 0.05:", sum(result[, "qval"] < 0.05), "\n")
  cat("SNPs significant at qval < 0.01:", sum(result[, "qval"] < 0.01), "\n")
  cat("p-value range: [", signif(min(result[, "pval"]), 3), ",", signif(max(result[, "pval"]), 3), "]\n")
  
  cat("\nTop 5 SNPs by association p-value:\n")
  print(head(result[order(result[, "pval"]), ], 5))
  
  cat("\nWrote:\n -", results_path, "\n -", file.path(output_dir, "test_boxplot.pdf"), "\n")
  
  invisible(result)
}

#When run as a standalone script (e.g. via Rscript), execute the analysis
#automatically using command line arguments for input/output directories, if given
if (!interactive() && sys.nframe() == 0) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  input_dir <- if (length(cli_args) >= 1) cli_args[1] else DEFAULT_DATA_DIR
  output_dir <- if (length(cli_args) >= 2) cli_args[2] else DEFAULT_OUTPUT_DIR
  run_snp_analysis(input_dir = input_dir, output_dir = output_dir)
}
