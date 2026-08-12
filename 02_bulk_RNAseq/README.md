# Gene Expression Analysis

R implementation for quality assessment, normalization, and differential expression analysis of human bone marrow gene expression count data.

## Contents

```text
02_bulk_RNAseq/
├── README.md
├── scripts/
│   └── gene_expression_analysis.R
│   └── test_gene_expression_analysis.R
├── data/
│   ├── raw_count.txt
│   └── gene_annot.txt
├── output/
│   ├── MvAplots.pdf
│   └── DE_results.txt
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

## Test and reproducibility

The analysis can be reproduced using the test script from the `02_bulk_RNAseq` directory:

``` bash
Rscript test_gene_expression_analysis.R
```

Input and output directories can also be specified:

``` bash
Rscript test_gene_expression_analysis.R data output
```

The script runs all three functions and stores the generated results in `output/`.

## Outputs

### MvA plots

`output/MvAplots.pdf` contains one MvA plot per sample comparison.

### Differential expression results

`output/DE_results.txt` contains the per-gene p-values, Benjamini–Hochberg q-values, and log2 fold changes, together with the summary statistics returned by `DEbyEdgeR()`.
