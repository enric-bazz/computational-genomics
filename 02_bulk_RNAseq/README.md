# Gene Expression Analysis

R implementation for quality assessment, normalization, and differential expression analysis of human bone marrow gene expression count data.

## Contents

```text
02_bulk_RNAseq/
├── README.md
├── scripts/
│   ├── gene_expression_analysis.R
│   └── test_gene_expression_analysis.R
├── data/                     (untracked)
│   ├── raw_count.txt
│   └── gene_annot.txt
└── output/
    ├── MvAplots.pdf
    ├── DE_summary.txt
    └── DE_results.txt
```

`scripts/gene_expression_analysis.R` contains the three required functions:

* `MvAplot()`
* `TMMnorm()`
* `DEbyEdgeR()`


## Functions

### `MvAplot()`

``` r
MvAplot <- function(exprData, pdffilename, pcolor = "black", lcolor = "red")
```

Generates MvA plots comparing each sample with the first sample. Each comparison is saved as a separate page in the specified PDF file.

### `TMMnorm()`

``` r
TMMnorm <- function(exprData, annot, Mtrim = 0.02, Atrim = c(0, 8))
```

Normalizes raw expression counts using a trimmed-mean approach, accounting for sequencing depth and gene length. Returns the normalized expression matrix and scaling factors.

### `DEbyEdgeR()`

``` r
DEbyEdgeR <- function(rawdata, groups, alpha = 0.05)
```

Performs differential expression analysis between patients and controls using `edgeR`. Returns per-gene p-values, Benjamini–Hochberg adjusted q-values, log2 fold changes, and summary estimates of selected genes, false positives, false negatives, and FDR.

## Assumptions

- Gene lengths are taken from the `Length` column of `annot` and matched to the count matrix by gene name.
- In `DEbyEdgeR()`, samples whose label contains `CTRL` are controls and all others are patients, ignoring the Group1/Group2 split as required.
- The instructions mention "a list of three elements" for `DEbyEdgeR()` but describe two; the function returns those two.

## Test and reproducibility

`scripts/test_gene_expression_analysis.R` provides a `run_gene_expression_analysis()` wrapper that runs all three functions and saves the results to the output directory.

From any directory:

``` sh
Rscript 02_bulk_RNAseq/scripts/test_gene_expression_analysis.R [input_dir] [output_dir]
```

Both directories are optional and default to `02_bulk_RNAseq/data` and `02_bulk_RNAseq/output`. The script can also be sourced in R and run with `run_gene_expression_analysis()`.

## Outputs

### MvA plots

`output/MvAplots.pdf` contains one MvA plot per sample, each compared with sample 1.

### Differential expression results

- `output/DE_summary.txt`: number of selected genes, expected false positives, expected false negatives, and FDR at `alpha = 0.05`.
- `output/DE_results.txt`: tab-delimited table with the per-gene p-value, Benjamini–Hochberg q-value, and log2 fold change, with gene names as row names.
