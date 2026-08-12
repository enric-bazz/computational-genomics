# Computational Genomics Analysis

This repository contains R-based analyses developed as part of the **Computational Genomics** course at **Department of Information Engineering (DEI), University of Padua**.

The repository contains a series of self-contained R workflows covering genotype, bulk gene-expression, and single-cell RNA-seq analyses. Each analysis folder includes the main analysis functions, test scripts for running them on the provided datasets, and a README documenting the inputs, outputs, and workflow. The analyses are designed to be reusable and reproducible, avoiding hard-coded dataset-specific assumptions where possible.

## Analysis overview

### 1. SNP analysis

The first analysis focuses on genotype data from case and control subjects. The implemented functions calculate minor allele frequencies, assess Hardy–Weinberg equilibrium, and perform SNP-level case-control association analysis.

The analysis includes multiple-testing correction using the Benjamini–Hochberg procedure and produces summary statistics for the SNPs passing the specified quality filters.

The corresponding input consists of a genotype matrix in which SNPs are represented by rows and subjects by columns.

### 2. Bulk gene-expression analysis

The second analysis focuses on raw gene-expression count data from human bone marrow samples, including control samples and leukemia patients belonging to two groups.

The workflow covers:

- MvA-based inspection of sample expression profiles;
- normalization of raw count data using a trimmed-mean approach;
- differential expression analysis using `edgeR`;
- statistical significance and multiple-testing assessment;
- calculation of fold changes and false-positive/false-negative estimates.

Gene annotation information, including gene lengths, is incorporated into the normalization procedure.

### 3. Single-cell RNA-seq analysis

The third analysis moves to single-cell data, using raw counts from a mouse artery single-cell RNA-seq experiment stored as a `SingleCellExperiment` object.

The workflow consists of:

- filtering and log-normalization;
- highly variable gene selection;
- PCA and t-SNE dimensionality reduction;
- cell clustering;
- reference-based cell-type annotation using `SingleR`;
- ligand–receptor-based cell–cell communication analysis.

A reference `SingleCellExperiment` is used for cell-type annotation, while a ligand–receptor database provides the interaction pairs evaluated in the final analysis.
