# scripts/01_make_dge.R

# user inputs -------------------------------------------------------------

# path to directory containing:
# BS_XXXXXXXX.rsem.genes.results.gz files

if (!exists("rsem_dir", inherits = FALSE)) {
  stop(
    "You must define `rsem_dir` before running the pipeline.\n",
    "Example:\n",
    "  rsem_dir <- '/path/to/rsem_gene_results'\n"
  )
}

if (!exists("meta_csv", inherits = FALSE)) {
  stop(
    "You must define `meta_csv` before running the pipeline.\n",
    "Example:\n",
    "  meta_csv <- 'data/sample_meta.csv'\n"
  )
}

# run ---------------------------------------------------------------------

# source the function(s) for this step
source("R/01_make_dge_from_rsem.R")

# run analysis
make_dge_from_rsem(
  rsem_dir = rsem_dir,
  meta_csv = meta_csv,
  out_rds  = "results/01_dge.rds"
)
