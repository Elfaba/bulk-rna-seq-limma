# scripts/00_run_all.R

if (!dir.exists("R") || !dir.exists("scripts")) {
  stop("Run this from the project root (MUST contain R/ and scripts/).")
}

message("Step 0: Load user config")
source("scripts/99_user_config.R")

message("Starting bulk RNA-seq analysis pipeline")

message("Step 1: Creating DGE object")
source("scripts/01_make_dge.R")

message("Step 2: Running limma-voom (DMG vs HGG)")
source("scripts/02_voom_limma.R")

message("Pipeline completed successfully")
