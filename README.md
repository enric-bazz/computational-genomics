# Computational Genomics Analysis — R Coursework

This repository contains R-based analyses developed as part of the **Computational Genomics** course at **Department of Information Engineering (DEI), University of Padua**.

The repository collects a series of practical computational analyses covering different stages and types of biological data analysis, with a particular focus on statistical methods and reproducible workflows in R. The exercises progress from classical genotype-based analysis to bulk gene-expression analysis and, finally, single-cell RNA-seq analysis.

The code is accompanied by test runs on the provided datasets and by the main outputs generated during the analyses.

## Repository structure

The repository is organized into separate folders, with each folder containing the code and results associated with one analysis:

```text
.
├── Homework1/
│   ├── R scripts
│   ├── input data
│   └── analysis outputs
│
├── Homework2/
│   ├── R scripts
│   ├── input data
│   └── analysis outputs
│
└── Homework3/
    ├── R scripts
    ├── input data
    └── analysis outputs
```

Each folder is self-contained and includes a dedicated README describing its analysis, implementation, inputs, and outputs.

## Analysis overview

### 1. SNP analysis

The first analysis focuses on genotype data from case and control subjects. The implemented functions calculate minor allele frequencies, assess Hardy–Weinberg equilibrium, and perform SNP-level case-control association analysis.

The analysis includes multiple-testing correction using the Benjamini–Hochberg procedure and produces summary statistics for the SNPs passing the specified quality filters.

The corresponding input consists of a genotype matrix in which SNPs are represented by rows and subjects by columns.

### 2. Bulk gene-expression analysis

The second analysis focuses on raw gene-expression count data from human bone marrow samples, including control samples and leukemia patients belonging to two groups.

The workflow covers:

* MvA-based inspection of sample expression profiles;
* normalization of raw count data using a trimmed-mean approach;
* differential expression analysis using `edgeR`;
* statistical significance and multiple-testing assessment;
* calculation of fold changes and false-positive/false-negative estimates.

Gene annotation information, including gene lengths, is incorporated into the normalization procedure.

### 3. Single-cell RNA-seq analysis

The third analysis moves to single-cell data, using raw counts from a mouse artery single-cell RNA-seq experiment stored as a `SingleCellExperiment` object.

The workflow consists of:

* filtering and log-normalization;
* highly variable gene selection;
* PCA and t-SNE dimensionality reduction;
* cell clustering;
* reference-based cell-type annotation using `SingleR`;
* ligand–receptor-based cell–cell communication analysis.

A reference `SingleCellExperiment` is used for cell-type annotation, while a ligand–receptor database provides the interaction pairs evaluated in the final analysis.

## General workflow

Across the repository, the analyses follow a common structure:

```text
Input biological data
        │
        ▼
Data processing and quality assessment
        │
        ▼
Statistical / computational analysis
        │
        ▼
Biological interpretation
        │
        ▼
Saved numerical and graphical outputs
```

The specific methods differ between datasets, reflecting the different biological questions and data types:

```text
Genotype data
    │
    └── SNP frequencies
         └── HWE testing
              └── SNP association analysis


Bulk RNA-seq counts
    │
    └── MvA assessment
         └── normalization
              └── differential expression


Single-cell RNA-seq counts
    │
    └── preprocessing
         └── dimensionality reduction
              └── clustering
                   └── cell-type annotation
                        └── cell–cell communication
```

## Code organization

The main R scripts contain the analysis functions for their respective datasets. The functions are designed to operate on the supplied data structures rather than relying on fixed dataset dimensions or hard-coded sample identities where the corresponding information is provided as an input parameter.

Separate test scripts are used to load the input data, invoke the analysis functions, and save the resulting outputs. This keeps the reusable analysis code separate from the code used to execute a particular test run.

## Outputs

Depending on the analysis, the repository contains both numerical and graphical outputs, including:

* SNP-level statistical results;
* HWE and MAF summaries;
* MvA plots;
* normalized expression data;
* differential expression statistics;
* PCA and t-SNE representations;
* cell clustering and annotation plots;
* ranked ligand–receptor interaction results.

The individual README files inside each analysis folder provide the details of the corresponding functions, input data, output files, and analysis workflow.

## Software and methods

The repository is primarily implemented in **R** and uses established bioinformatics packages and statistical procedures appropriate to the different data types.

The analyses include methods for:

* genotype and SNP analysis;
* gene-expression normalization and differential expression;
* dimensionality reduction and clustering;
* reference-based single-cell annotation;
* ligand–receptor-based cell–cell communication analysis.

The purpose of the repository is to collect these analyses in a coherent, reproducible structure while keeping each workflow independently executable and documented.
