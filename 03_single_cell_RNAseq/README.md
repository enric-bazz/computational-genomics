# Single-Cell RNA-seq Analysis

This folder contains an R-based single-cell RNA-seq analysis pipeline for mouse artery data. The workflow starts from raw-count `SingleCellExperiment` data, performs preprocessing and dimensionality reduction, identifies cell populations and assigns cell-type annotations, and evaluates potential cell–cell communication using a ligand–receptor database.

## Context and analysis goal

The input consists of raw single-cell RNA-seq counts from mouse arteries. A separate reference dataset provides cell-type annotations, while a ligand–receptor database provides known ligand–receptor gene pairs.

The analysis comprises three main stages:

1.  **Preprocessing and dimensionality reduction** — filtering, normalization, highly variable gene selection, PCA, and t-SNE.
2.  **Cell-type annotation** — Louvain clustering followed by cluster-level annotation with `SingleR`.
3.  **Cell–cell communication analysis** — expression-based scoring of ligand–receptor pairs between annotated cell populations.

The runner script also performs a separate t-SNE perplexity comparison using values of 50, 100, and 400.

## Contents

``` text
03_single_cell_RNAseq/
├── README.md
├── scripts/
│   ├── scRNA_analysis.R
│   └── test_scRNA_analysis.R
├── data/
│   ├── sce_count.RData
│   ├── ref.RData
│   └── LRdb.RData
└── outputs/
    ├── tSNE_annotation.pdf
    ├── tSNE_perplexity_comparison.pdf
    └── CCC_results.txt
```

## Functions

### `process_sce()`

``` r
process_sce <- function(sce, n_HVG = 1000, percVar = 70)
```

The function:

- removes genes with zero expression across all cells;
- normalizes the raw counts and stores the result in `logcounts`;
- identifies the top highly variable genes using variance modelling;
- initially computes up to 100 PCs using the selected HVGs;
- determines the minimum number of PCs required to explain at least `percVar` percent of the variance;
- recomputes PCA using only the required number of PCs;
- computes a t-SNE representation using the selected HVGs.

The selected HVGs are stored in the `HVGs` metadata field.

### `annotation()`

``` r
annotation <- function(sce, ref, pdffilename)
```

Cells are clustered using a 10-nearest-neighbor SNN graph and the Louvain algorithm. The resulting clusters are annotated with `SingleR` using the supplied reference dataset.

The function adds the following fields to the `SingleCellExperiment`:

- `Louvain` — cluster assignments;
- `cell_annotation` — assigned cell types.

It also generates a PDF containing t-SNE plots colored by Louvain cluster and cell-type annotation.

### `CCCanalysis()`

``` r
CCCanalysis <- function(sce, LRdb)
```

For every sender–receiver cell-type combination and ligand–receptor pair, the function calculates:

- mean log-normalized ligand expression in the sender population;
- mean log-normalized receptor expression in the receiver population;
- `LRscore`, defined as their product.

The output contains:

``` text
ligand | receptor | sender | receiver | LRscore
```

and is sorted in descending order of `LRscore`.

## Testing and Reproducing

The analysis is implemented in `scripts/scRNA_analysis.R` and executed through `scripts/test_scRNA_analysis.R`. The runner accepts the input and output directories as command-line arguments.

The required input files are:

- `data/sce_count.RData`
- `data/ref.RData`
- `data/LRdb.RData`

Run the complete analysis from the `03_single_cell_RNAseq/` directory with:

``` text
Rscript scripts/test_scRNA_analysis.R data outputs
```

## t-SNE perplexity analysis

The runner additionally evaluates the effect of t-SNE perplexity by computing independent embeddings with:

``` text
50, 100, 400
```

The same processed data, PCA representation, and Louvain assignments are used for all three embeddings. Only the t-SNE perplexity is changed, allowing its effect on the visualization of the inferred cell populations to be compared directly.

## Interpretation

The cell–cell communication analysis provides an **expression-based prioritization of potential interactions**, not direct evidence of cell–cell communication.

A high `LRscore` means that the corresponding ligand has relatively high mean expression in the proposed sender population and the receptor has relatively high mean expression in the proposed receiver population. It does **not** demonstrate that the ligand–receptor interaction occurs biologically, that the proteins are active, or that the two populations actually communicate.

Similarly, the t-SNE perplexity comparison is intended to assess the stability and visual organization of the low-dimensional representation. Differences between embeddings should not be interpreted as evidence for or against biological populations by themselves.
