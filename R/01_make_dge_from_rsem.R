# R/01_make_dge_from_rsem.R

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
  library(stringr)
})


# helpers -----------------------------------------------------------------

# is `meta.csv` missing any required columns?
stop_if_missing_cols <- function(df, cols, df_name = "data.frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      df_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

# does the IDs in `meta.csv` follow the naming rules?
# patient_id: PT_ followed by exactly 8 alphanumeric characters
# biospecimen_id: BS_ followed by exactly 8 alphanumeric characters
validate_ids <- function(meta) {
  bad_pt <- meta$patient_id[
    !str_detect(meta$patient_id, "^PT_[A-Za-z0-9]{8}$")
  ]
  bad_bs <- meta$biospecimen_id[
    !str_detect(meta$biospecimen_id, "^BS_[A-Za-z0-9]{8}$")
  ]
  
  if (length(bad_pt) > 0) {
    stop(
      "Invalid patient ID format. Example bad value: ",
      bad_pt[1],
      "\nExpected: PT_ followed by exactly 8 alphanumeric characters."
    )
  }
  if (length(bad_bs) > 0) {
    stop(
      "Invalid biospecimen ID format. Example bad value: ",
      bad_bs[1],
      "\nExpected: BS_ followed by exactly 8 alphanumeric characters."
    )
  }
}

# list all files in rsem_dir that matches the naming convention suffix
# note: do not print the file list to avoid bloating output
list_rsem_files <- function(rsem_dir) {
  files <- list.files(
    rsem_dir,
    pattern = "\\.rsem\\.genes\\.results\\.gz$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop(
      "No files ending with '.rsem.genes.results.gz' found in: ",
      rsem_dir
    )
  }
  files
}

# extract the BS_XXXXXXXX prefix from each filename
extract_bsid_from_filename <- function(path) {
  fname <- basename(path)
  bsid <- str_match(
    fname,
    "^(BS_[A-Za-z0-9]{8})\\.rsem\\.genes\\.results\\.gz$"
  )[, 2]
  bsid
}

# format a short ID preview for messages (avoid huge console output)
preview_ids <- function(x, n = 20) {
  x <- unique(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return("(none)")
  out <- head(x, n)
  if (length(x) > n) {
    paste0(
      paste(out, collapse = ", "),
      ", ... (",
      length(x) - n,
      " more)"
    )
  } else {
    paste(
      out,
      collapse = ", "
    )
  }
}


# main function to make dge -----------------------------------------------

make_dge_from_rsem <- function(
    rsem_dir,
    meta_csv,
    out_rds = "results/01_dge.rds",
    rsem_columns = c(1, 5)  # 1=gene_id, 5=expected_count
) {
  if (!dir.exists(rsem_dir)) stop("rsem_dir not found: ", rsem_dir)
  if (!file.exists(meta_csv)) stop("meta.csv not found: ", meta_csv)
  
  meta <- readr::read_csv(meta_csv, show_col_types = FALSE)
  
  required <- c("patient_id", "biospecimen_id", "tumor_status", "subtype")
  stop_if_missing_cols(meta, required, "meta.csv")
  
  meta <- meta %>%
    mutate(
      patient_id = as.character(patient_id),
      biospecimen_id = as.character(biospecimen_id),
      tumor_status = as.character(tumor_status),
      subtype = as.character(subtype)
    )
  
  validate_ids(meta)
  
  dup_bs <- unique(meta$biospecimen_id[duplicated(meta$biospecimen_id)])
  
  if (length(dup_bs) > 0) {
    stop(
      "Duplicated biospecimen_id found in meta.csv: ",
      length(dup_bs),
      " (showing up to 20): ",
      preview_ids(dup_bs, n = 20)
    )
  }
  
  # list rsem files and map them to bsids
  paths <- list_rsem_files(rsem_dir)
  bsids <- vapply(paths, extract_bsid_from_filename, character(1))
  
  if (any(is.na(bsids))) {
    bad <- basename(paths[is.na(bsids)])[1]
    stop(
      "Some files do not match expected naming:\n",
      "BS_XXXXXXXX.rsem.genes.results.gz\n",
      "Example: ",
      bad
    )
  }
  
  # report mismatches
  bs_in_files <- unique(bsids)
  bs_in_meta <- unique(meta$biospecimen_id)
  
  bs_missing_in_meta <- setdiff(bs_in_files, bs_in_meta)
  bs_missing_in_files <- setdiff(bs_in_meta, bs_in_files)
  
  if (length(bs_missing_in_meta) > 0) {
    message(
      "BSID(s) present in RSEM folder but missing from meta.csv: ",
      length(bs_missing_in_meta),
      " (showing up to 20): ",
      preview_ids(bs_missing_in_meta, n = 20)
    )
  }
  
  if (length(bs_missing_in_files) > 0) {
    message(
      "BSID(s) present in meta.csv but missing from RSEM folder: ",
      length(bs_missing_in_files),
      " (showing up to 20): ",
      preview_ids(bs_missing_in_files, n = 20)
    )
  }
  
  # keep only files that have metadata
  keep <- bsids %in% meta$biospecimen_id
  if (!any(keep)) {
    stop(
      "None of the BSIDs from filenames match meta$biospecimen_id.\n",
      "Example BSID from files: ", bsids[1], "\n",
      "Example BSID from meta: ", meta$biospecimen_id[1]
    )
  }
  
  paths <- paths[keep]
  bsids <- bsids[keep]
  
  # read DGE
  y <- edgeR::readDGE(paths, columns = rsem_columns)
  
  # ensure sample names are BSIDs
  colnames(y$counts) <- bsids
  
  stopifnot(ncol(y$counts) == length(bsids))
  
  # align meta to counts column order
  meta2 <- meta %>%
    filter(biospecimen_id %in% bsids) %>%
    distinct(biospecimen_id, .keep_all = TRUE)
  
  meta2 <- meta2[match(bsids, meta2$biospecimen_id), ]
  stopifnot(all(meta2$biospecimen_id == bsids))
  
  # attach metadata to y$samples (edgeR convention)
  y$samples$patient_id     <- meta2$patient_id
  y$samples$biospecimen_id <- meta2$biospecimen_id
  y$samples$tumor_status   <- meta2$tumor_status
  y$samples$subtype        <- meta2$subtype
  
  # save for next step
  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(dge = y, meta = meta2), out_rds)
  
  message("Saved DGE and aligned metadata to: ", out_rds)
  
  session_file <- file.path(dirname(out_rds), "01_sessionInfo.txt")
  
  capture.output(
    if (requireNamespace("BiocManager", quietly = TRUE)) BiocManager::version(),
    sessionInfo(),
    file = session_file
  )
  
  invisible(list(dge = y, meta = meta2))
}