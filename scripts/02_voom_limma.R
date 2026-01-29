# scripts/02_voom_limma.R


# user inputs -------------------------------------------------------------

in_rds <- "results/01_dge.rds"

# where to write outputs (prefix used by the function)
out_prefix <- "results/02"

# optional: to change filtering thresholds for genes
min_count <- 10
min_prop  <- 0.7


# run ---------------------------------------------------------------------

if (!file.exists(in_rds)) {
  stop(
    "Input not found: ",
    in_rds,
    "\n",
    "Should run scripts/01_make_dge.R first"
  )
}

# source the function(s) for this step
source("R/02_voom_limma_dmg_vs_hgg.R")

# run analysis
run_voom_limma_dmg_vs_hgg(
  in_rds    = in_rds,
  out_prefix = out_prefix,
  min_count  = min_count,
  min_prop   = min_prop
)