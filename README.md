---
editor_options: 
  markdown: 
    wrap: 72
---

# bulk-rna-seq-limma

Repository for bulk RNA-seq differential expression analysis using
**edgeR + limma-voom**.

This repository is designed for **reproducible analysis and teaching**.
No expression data are included.

------------------------------------------------------------------------

## Repository structure

├── R/ \# analysis functions\
├── scripts/ \# executable pipeline steps\
├── data/\
├── results/ \# output (git-ignored)\
├── bulk-rna-seq-limma.Rproj

------------------------------------------------------------------------

## Software environment

This pipeline was tested with the following environment. Exact session
information is written to `results/*_sessionInfo.txt` when the pipeline
is run.

### R

-   **R version:** 4.4.3
-   **Bioconductor:** 3.20
-   **Platform tested:** Windows 11 (x86_64)

------------------------------------------------------------------------

## Required packages

### Core Bioconductor packages

-   **edgeR** (≥ 4.4.2)
-   **limma** (≥ 3.62.2)

### Tidyverse / plotting

-   dplyr (1.1.4)
-   readr (2.1.6)
-   stringr (1.6.0)
-   tidyr (1.3.2)
-   tibble (3.3.1)
-   ggplot2 (4.0.1)

------------------------------------------------------------------------

## Installation

``` r
install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install(c("edgeR", "limma"))
```

------------------------------------------------------------------------

## Usage

### 1. Configure inputs

Create or edit:

``` r
scripts/99_user_config.R
```

Example:

``` r
rsem_dir <- "/path/to/rsem_gene_results"
meta_csv <- "data/sample_meta.csv"
```

This file is intentionally git-ignored to avoid committing
machine-specific paths.

### 2. Run the full pipeline

From the project root:

``` r
source("scripts/00_run_all.R")
```

#### Pipeline steps:

1.  Build `DGEList` from RSEM gene counts
2.  Perform limma-voom differential expression (DMG vs HGG)
3.  Write results, plots, and session information to `results/`

------------------------------------------------------------------------

## Input requirements

### RSEM gene result files

-   One file per sample
-   Naming convention **required**: `BS_XXXXXXXX.rsem.genes.results.gz`

### Metadata

Required columns:

-   `patient_id`
-   `biospecimen_id`
-   `tumor_status`
-   `subtype` (e.g., DMG, HGG)

------------------------------------------------------------------------

## Outputs

Written to `results/` (ignored by Git):

-   Filtered and normalized DGE objects (`.rds`)
-   Differential expression table (`.csv`)
-   MDS plot (`.png`)
-   Session information (`*_sessionInfo.txt`)

------------------------------------------------------------------------

## Notes

-   Results are excluded from version control
-   User-specific configuration is excluded
-   Designed as ***template*** for reuse, extension, and teaching
