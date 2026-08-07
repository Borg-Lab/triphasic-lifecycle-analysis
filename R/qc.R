build_qc_metrics <- function(counts_df, tpm_df, design_df, analysis_cfg) {
  out <- data.frame(
    sample_id = design_df$sample_id,
    generation = design_df$generation,
    stage = design_df$stage,
    tissue = design_df$tissue,
    sex = design_df$sex,
    replicate = design_df$replicate,
    batch = design_df$batch,
    library_size = NA_real_,
    detected_genes = NA_integer_,
    stringsAsFactors = FALSE
  )

  if (nrow(counts_df) > 0) {
    for (i in seq_len(nrow(out))) {
      s <- out$sample_id[i]
      if (s %in% colnames(counts_df)) {
        out$library_size[i] <- sum(as.numeric(counts_df[[s]]), na.rm = TRUE)
      }
    }
  }

  if (nrow(tpm_df) > 0) {
    cutoff <- as.numeric(analysis_cfg$tpm_detect_threshold)
    for (i in seq_len(nrow(out))) {
      s <- out$sample_id[i]
      if (s %in% colnames(tpm_df)) {
        out$detected_genes[i] <- sum(as.numeric(tpm_df[[s]]) >= cutoff, na.rm = TRUE)
      }
    }
  }

  out
}

compute_group_tpm_means <- function(tpm_df, design_df, group_col = "stage") {
  if (nrow(tpm_df) == 0) return(data.frame(gene_id = character()))
  if (!group_col %in% colnames(design_df)) stop(sprintf("group_col not in design: %s", group_col), call. = FALSE)

  groups <- unique(design_df[[group_col]])
  out <- data.frame(gene_id = tpm_df$gene_id, stringsAsFactors = FALSE)
  for (g in groups) {
    s <- design_df$sample_id[design_df[[group_col]] == g & as.logical(design_df$include)]
    s <- intersect(s, colnames(tpm_df))
    cname <- paste0("mean_tpm_", make.names(g))
    if (length(s) == 0) {
      out[[cname]] <- NA_real_
    } else if (length(s) == 1) {
      out[[cname]] <- as.numeric(tpm_df[[s]])
    } else {
      out[[cname]] <- rowMeans(as.matrix(tpm_df[, s, drop = FALSE]), na.rm = TRUE)
    }
  }
  out
}
