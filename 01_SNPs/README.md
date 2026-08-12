# SNP Analysis

This folder contains an R implementation for basic SNP-level genetic analysis, together with a test run on the provided genotype dataset.

## Contents

``` text
01_SNPs/
├── README.md
├── scripts/
│   ├── SNP_analysis.R
│   └── test_SNP_analysis.R
├── data/
│   └── SNPdata.txt (untracked)
└── output/
    ├── association_results.txt
    └── test_boxplot.pdf
```

`scripts/SNP_analysis.R` contains the three required functions:

- `compute_MAF()`
- `HWE_test()`
- `SNP_association_test()`

## Data

The input dataset is a tab-delimited genotype matrix in which rows represent SNPs and columns represent subjects. Genotypes are encoded as:

- `0` — homozygous reference genotype (`AA`)
- `1` — heterozygous genotype (`Aa`)
- `2` — homozygous alternative genotype (`aa`)

The implementation does not depend on a fixed number of SNPs or subjects, nor on the names or ordering of the samples. The control samples are identified through the `indCTRL` argument of `SNP_association_test()`.

## Functions

### `compute_MAF()`

``` r
compute_MAF <- function(SNPdata)
```

Computes the minor allele frequency for every SNP in the supplied genotype matrix.

The function:

- determines the number of SNPs and subjects from the input matrix;
- reports these values to the terminal;
- calculates the MAF for each SNP;
- preserves the SNP identifiers and their original order in the returned vector.

### `HWE_test()`

``` r
HWE_test <- function(SNPdata, pdf_file_name)
```

Performs a Hardy–Weinberg equilibrium test independently for each SNP.

For each SNP, the function calculates the observed chi-square statistic and derives its p-value using `pchisq()`. It reports the SNP(s) with the lowest p-value and returns a named vector containing the HWE p-values.

The distribution of the resulting p-values is also saved as a boxplot to the path specified by `pdf_file_name`.

### `SNP_association_test()`

``` r
SNP_association_test <- function(filepath, indCTRL, MAFth=0.01, HWEalpha=0.01)
```

Performs SNP-level case-control association analysis starting from a genotype data file.

The function:

1.  reads the genotype matrix from `filepath`;
2.  calculates MAF values and HWE p-values using the corresponding functions;
3.  filters SNPs according to the specified MAF and HWE thresholds;
4.  performs a chi-square association test between genotype and case/control status;
5.  calculates multiple-testing-adjusted q-values using the Benjamini–Hochberg procedure;
6.  returns the genotype counts and statistical results for the SNPs passing the filters.

The returned matrix contains the following columns:

``` text
AA_ctrl | Aa_ctrl | aa_ctrl | AA_case | Aa_case | aa_case | pval | qval
```

SNP identifiers are retained as row names.

## Test and reproducibility

The analysis was tested on the provided genotype dataset using control samples `1201:2000` and the default thresholds (`MAFth = 0.01`, `HWEalpha = 0.01`).

``` r
SNP_association_test(
    filepath = ".../SNPdata.txt",
    indCTRL = 1201:2000,
    MAFth = 0.01,
    HWEalpha = 0.01
)
```

The test run also generates the resulting association matrix, saved in `association_results.txt`, and `test_boxplot.pdf`, containing the distribution of HWE-test p-values.

`scripts/test_SNP_analysis.R` provides a `run_snp_analysis()` wrapper that runs the analysis and saves `association_results.txt` and `test_boxplot.pdf` to the output directory.

From `01_SNPs/scripts/`:

``` sh
Rscript test_SNP_analysis.R [input_dir] [output_dir]
```

or from R/RStudio:

``` r
source("01_SNPs/scripts/test_SNP_analysis.R")
run_snp_analysis()
run_snp_analysis(input_dir = "path/to/data", output_dir = "path/to/results")
```

Both directories are optional and default to `01_SNPs/data` and `01_SNPs/output`; `input_dir` must contain `SNPdata.txt`. The wrapper also accepts `data_file`, `indCTRL`, `MAFth`, and `HWEalpha` for custom runs.

## Implementation notes

The code operates on the structure of the genotype matrix rather than assuming particular SNP or sample names. Consequently, the functions can be applied to datasets with different numbers of SNPs and subjects and with controls and cases occurring in arbitrary column positions, provided that the corresponding control indices are supplied.
