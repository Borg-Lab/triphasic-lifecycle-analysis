run_module_inference <- function(counts_df, analysis_cfg) {
  schema <- data.frame(
    gene_id = character(), module_id = character(), kME = numeric(), hub_flag = logical(),
    stringsAsFactors = FALSE
  )

  if (!isTRUE(analysis_cfg$wgcna_enable)) return(schema)
  if (nrow(counts_df) < 2) return(schema)
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required when wgcna_enable is TRUE.", call. = FALSE)
  }

  mat <- as.matrix(counts_df[, setdiff(colnames(counts_df), "gene_id"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- counts_df$gene_id
  mat <- round(mat)
  keep <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 50 || ncol(mat) < 4) return(schema)

  # genes in columns, samples in rows for WGCNA
  datExpr <- t(log2(mat + 1))

  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  if (!isTRUE(gsg$allOK)) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  if (ncol(datExpr) < 50 || nrow(datExpr) < 4) return(schema)

  if (is.function(WGCNA::allowWGCNAThreads)) {
    suppressWarnings(WGCNA::allowWGCNAThreads(nThreads = 8))
  }

  power <- if (!is.null(analysis_cfg$wgcna_power)) as.integer(analysis_cfg$wgcna_power) else 12L
  min_module_size <- if (!is.null(analysis_cfg$wgcna_min_module_size)) as.integer(analysis_cfg$wgcna_min_module_size) else 30L
  merge_cut_height <- if (!is.null(analysis_cfg$wgcna_merge_cut_height)) as.numeric(analysis_cfg$wgcna_merge_cut_height) else 0.25
  deep_split <- if (!is.null(analysis_cfg$wgcna_deep_split)) as.integer(analysis_cfg$wgcna_deep_split) else 2L
  max_block_size <- if (!is.null(analysis_cfg$wgcna_max_block_size)) as.integer(analysis_cfg$wgcna_max_block_size) else 5000L

  # Compatibility fix for some WGCNA builds expecting weighted cor in scope.
  temp_cor <- cor
  cor <- WGCNA::cor
  on.exit({
    cor <- temp_cor
  }, add = TRUE)

  net <- WGCNA::blockwiseModules(
    datExpr = datExpr,
    power = power,
    networkType = "signed",
    TOMType = "signed",
    deepSplit = deep_split,
    pamRespectsDendro = FALSE,
    minModuleSize = min_module_size,
    maxBlockSize = max_block_size,
    reassignThreshold = 0,
    mergeCutHeight = merge_cut_height,
    saveTOMs = FALSE,
    numericLabels = TRUE,
    verbose = 0
  )

  module_colors <- WGCNA::labels2colors(net$colors)
  names(module_colors) <- colnames(datExpr)
  MEs <- WGCNA::moduleEigengenes(datExpr, colors = module_colors)$eigengenes
  MEs <- WGCNA::orderMEs(MEs)
  kme_mat <- WGCNA::signedKME(datExpr, MEs, outputColumnName = "kME")

  genes <- colnames(datExpr)
  out <- data.frame(
    gene_id = genes,
    module_id = paste0("module_", module_colors),
    kME = NA_real_,
    hub_flag = FALSE,
    stringsAsFactors = FALSE
  )

  for (m in unique(module_colors)) {
    idx <- which(module_colors == m)
    kcol <- paste0("kME", m)
    if (!kcol %in% colnames(kme_mat)) next
    out$kME[idx] <- as.numeric(kme_mat[idx, kcol])
    cutoff <- stats::quantile(abs(out$kME[idx]), probs = 0.9, na.rm = TRUE)
    out$hub_flag[idx] <- abs(out$kME[idx]) >= cutoff
  }

  out
}
