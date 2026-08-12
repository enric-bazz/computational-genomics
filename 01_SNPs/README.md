# SNP Analysis

This folder contains an R implementation for basic SNP-level genetic analysis, together with a test run on the provided genotype dataset.

## Contents

```text
01_SNPs/
├── README.md
├── Homework1_Instructions.pdf
├── scripts/
│   ├── SNP_analysis.R       — the three analysis functions (the graded submission)
│   └── test_SNP_analysis.R  — runner script that reproduces the test run below
├── data/
│   └── SNPdata.txt          — genotype data used for the test run
└── output/
    ├── association_results.txt  — output of the SNP association analysis
    └── test_boxplot.pdf         — boxplot of the HWE-test p-values from the test run
```

`scripts/SNP_analysis.R` contains the three required functions:

* `compute_MAF()`
* `HWE_test()`
* `SNP_association_test()`

`data/` and `output/` are separate from `scripts/` so that `SNP_analysis.R` remains the single, self-contained file expected for submission — it takes a `filepath` argument rather than assuming any particular folder layout.

## Data

The input dataset is a tab-delimited genotype matrix in which rows represent SNPs and columns represent subjects. Genotypes are encoded as:

* `0` — homozygous reference genotype (`AA`)
* `1` — heterozygous genotype (`Aa`)
* `2` — homozygous alternative genotype (`aa`)

The implementation does not depend on a fixed number of SNPs or subjects, nor on the names or ordering of the samples. The control samples are identified through the `indCTRL` argument of `SNP_association_test()`.

## Functions

### `compute_MAF()`

```r
compute_MAF <- function(SNPdata)
```

Computes the minor allele frequency for every SNP in the supplied genotype matrix.

The function:

* determines the number of SNPs and subjects from the input matrix;
* reports these values to the terminal;
* calculates the MAF for each SNP;
* preserves the SNP identifiers and their original order in the returned vector.

### `HWE_test()`

```r
HWE_test <- function(SNPdata, pdf_file_name)
```

Performs a Hardy–Weinberg equilibrium test independently for each SNP.

For each SNP, the function calculates the observed chi-square statistic and derives its p-value using `pchisq()`. It reports the SNP(s) with the lowest p-value and returns a named vector containing the HWE p-values.

The distribution of the resulting p-values is also saved as a boxplot to the path specified by `pdf_file_name`.

### `SNP_association_test()`

```r
SNP_association_test <- function(filepath, indCTRL, MAFth=0.01, HWEalpha=0.01)
```

Performs SNP-level case-control association analysis starting from a genotype data file.

The function:

1. reads the genotype matrix from `filepath`;
2. calculates MAF values and HWE p-values using the corresponding functions;
3. filters SNPs according to the specified MAF and HWE thresholds;
4. performs a chi-square association test between genotype and case/control status;
5. calculates multiple-testing-adjusted q-values using the Benjamini–Hochberg procedure;
6. returns the genotype counts and statistical results for the SNPs passing the filters.

The returned matrix contains the following columns:

```text
AA_ctrl
Aa_ctrl
aa_ctrl
AA_case
Aa_case
aa_case
pval
qval
```

SNP identifiers are retained as row names.

## Test run

The functions were tested using the provided genotype dataset, with the control samples specified through their column indices. The default filtering thresholds were used:

```r
SNP_association_test(
    filepath = ".../SNPdata.txt",
    indCTRL = 1201:2000,
    MAFth = 0.01,
    HWEalpha = 0.01
)
```

The test run also generates `test_boxplot.pdf`, containing the distribution of HWE-test p-values.

The resulting association matrix is saved in `association_results.txt` for inspection and reproducibility.

## How to reproduce

`scripts/test_SNP_analysis.R` wraps the test run above in a `run_snp_analysis()` entry point that sources `SNP_analysis.R`, runs `SNP_association_test()` on a genotype file, and writes `association_results.txt` and `test_boxplot.pdf` to an output directory. Both the input and output directories default to the sibling `data/` and `output/` folders, so it works out of the box from this repository layout.

From the command line (run from `01_SNPs/scripts/`):

```sh
Rscript test_SNP_analysis.R [input_dir] [output_dir]
```

`input_dir` and `output_dir` are optional and default to `01_SNPs/data` and `01_SNPs/output`; `input_dir` must contain `SNPdata.txt`.

From an interactive R/RStudio session:

```r
source("01_SNPs/scripts/test_SNP_analysis.R")
run_snp_analysis()                                          # use all defaults (data/ -> output/)
run_snp_analysis(input_dir = "path/to/data", output_dir = "path/to/results")
```

`run_snp_analysis()` also accepts `data_file`, `indCTRL`, `MAFth`, and `HWEalpha` to rerun the analysis with different inputs or thresholds, and prints a short summary (number of SNPs passing filters, number significant at common q-value cutoffs, and the top hits by p-value) to the console.

## Implementation notes

The code operates on the structure of the genotype matrix rather than assuming particular SNP or sample names. Consequently, the functions can be applied to datasets with different numbers of SNPs and subjects and with controls and cases occurring in arbitrary column positions, provided that the corresponding control indices are supplied.
