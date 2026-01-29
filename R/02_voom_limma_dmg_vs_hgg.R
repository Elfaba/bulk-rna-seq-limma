# R/02_voom_limma_dmg_vs_hgg.R

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})


# helpers -----------------------------------------------------------------

# format a short ID preview for messages
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

# build gene_id to gene_name reference
# from the RSEM gene_id strings
make_gene_ref <- function(rn) {
  tibble(raw = rn) %>%
    separate(
      raw,
      into = c("gene_id", "gene_name"),
      sep = "_",
      extra = "merge",
      remove = FALSE
    ) %>%
    mutate(gene_id = str_replace(gene_id, "\\..*$", "")) %>%
    distinct(raw, .keep_all = TRUE)
}

# choose one biospecimen per patient based on mean expression
pick_one_sample_per_patient <- function(y, patient_col = "patient_id") {
  stopifnot(inherits(y, "DGEList"))
  if (!patient_col %in% names(y$samples)) {
    stop(
      "patient_col not found in y$samples: ",
      patient_col
    )
  }
  
  # mean expression per sample on raw counts
  mean_expr <- colMeans(y$counts)
  
  df <- tibble(
    bsid = colnames(y$counts),
    patient_id = y$samples[[patient_col]],
    mean_expr = mean_expr
  )
  
  dup_patients <- df %>%
    count(patient_id) %>%
    filter(n > 1) %>%
    pull(patient_id)
  
  if (length(dup_patients) > 0) {
    message(
      "Duplicate patient_id detected: ",
      length(dup_patients),
      " (showing up to 20): ",
      preview_ids(dup_patients, n = 20),
      "\nResolve to one sample per patient using mean expression"
    )
  } else {
    message("No duplicate patient_id detected.")
  }
  
  keep_bsid <- df %>%
    group_by(patient_id) %>%
    arrange(desc(mean_expr)) %>%
    slice(1) %>%
    ungroup() %>%
    pull(bsid)
  
  y[, keep_bsid, keep.lib.sizes = FALSE]
}


# main function to run limma-voom -----------------------------------------

run_voom_limma_dmg_vs_hgg <- function(
    in_rds = "results/01_dge.rds",
    out_prefix = "results/02",
    min_count = 10,
    min_prop  = 0.7
) {
  if (!file.exists(in_rds)) stop("Input RDS not found: ", in_rds)
  
  obj <- readRDS(in_rds)
  y <- obj$dge

  # expectation: subtype to exist and be only DMG/HGG
  if (!"subtype" %in% names(y$samples)) stop("subtype not found in y$samples.")
  
  y$samples$subtype <- as.character(y$samples$subtype)
  
  bad <- setdiff(unique(y$samples$subtype), c("DMG", "HGG"))
  
  if (length(bad) > 0) {
    stop(
      "Unexpected subtype value(s): ",
      paste(bad, collapse = ", "),
      "\nExpecting only: DMG, HGG"
    )
  }
  
  # gene_id -> gene_name reference (from rownames)
  gene_ref <- make_gene_ref(rownames(y))
  
  dir.create(
    dirname(out_prefix),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    gene_ref,
    paste0(out_prefix, "_gene_ref.csv"),
    row.names = FALSE
  )
  
  # resolve duplicate patients by choosing highest-mean-expression sample
  y0 <- y
  y  <- pick_one_sample_per_patient(y, patient_col = "patient_id")
  
  message("Samples before duplicate resolution: ", ncol(y0$counts))
  message("Samples after  duplicate resolution: ", ncol(y$counts))
  
  # begin limma workflow
  group <- factor(
    y$samples$subtype,
    levels = c("DMG", "HGG")
  )
  
  stopifnot(all(!is.na(group)))
  
  keep <- filterByExpr(
    y,
    group = group,
    min.count = min_count,
    min.prop  = min_prop
  )
  
  y2 <- y[keep, , keep.lib.sizes = FALSE]
  
  message("Genes before filtering: ", nrow(y$counts))
  message("Genes after  filtering: ", nrow(y2$counts))
  
  y3 <- calcNormFactors(y2, method = "TMM")
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # voom:
  # transforms counts to log2-CPM
  # estimates a mean–variance relationship
  #   from model residuals if design is supplied
  #   assigns precision weights
  # voom DOES NOT adjust expression values by group
  # v$E stores log2CPM
  # v$E is identical regardless of whether a design matrix is supplied
  # v$weights stores precision weights
  # v$weights affect statistical testing:
  #   - standard errors
  #   - moderated t-statistics (via empirical Bayes shrinkage)
  #   - p-values
  #   - FDR
  #   - power and stability
  
  v <- voom(y3, design = design, plot = TRUE)
  
  fit <- lmFit(v, design)
  
  contr <- makeContrasts(
    DMG_vs_HGG = DMG - HGG,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contr)
  
  fit2 <- eBayes(fit2)
  
  tt <- topTable(
    fit2,
    coef = "DMG_vs_HGG",
    number = Inf,
    sort.by = "P"
  )
  
  tt$gene_id_raw <- rownames(tt)
  
  # add gene names and gene ids to result
  ann <- make_gene_ref(tt$gene_id_raw)
  
  tt <- tt %>%
    tibble::rownames_to_column("rowname") %>%
    left_join(ann, by = c("rowname" = "raw")) %>%
    select(-rowname)
  
  write.csv(
    tt,
    paste0(out_prefix, "_limma_DMG_vs_HGG.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    list(dge = y3, voom = v, fit = fit2, topTable = tt),
    paste0(out_prefix, "_voom_limma_DMG_vs_HGG.rds")
  )
  
  # mds plot
  mds <- plotMDS(v$E, plot = FALSE, dim.plot = c(1, 2))
  
  df_mds <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    subtype = group
  )
  
  p <- ggplot(
    df_mds,
    aes(
      Dim1,
      Dim2,
      fill = subtype
    )
  ) +
    geom_point(
      shape = 21,
      size = 2,
      alpha = 0.8,
      color = "grey15",
      stroke = 0.25
    ) +
    coord_fixed() +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    labs(fill = "Subtype")
  
  ggsave(
    filename = paste0(out_prefix, "_plot_MDS_color_subtype.png"),
    plot = p,
    width = 4,
    height = 4,
    dpi = 300
  )
  
  message("Wrote outputs with prefix: ", out_prefix)
  
  capture.output(
    if (requireNamespace("BiocManager", quietly = TRUE)) BiocManager::version(),
    sessionInfo(),
    file = paste0(out_prefix, "_sessionInfo.txt")
  )
  
  invisible(
    list(
      dge = y3,
      voom = v,
      fit = fit2,
      topTable = tt,
      gene_ref = gene_ref
    )
  )
}