run_differential_expression <- function(counts_df, design_df, contrast_df) {
  schema <- data.frame(
    gene_id = character(), contrast_id = character(),
    log2fc = numeric(), p_value = numeric(), fdr = numeric(), base_mean = numeric(),
    stringsAsFactors = FALSE
  )

  if (nrow(counts_df) == 0 || nrow(contrast_df) == 0) return(schema)
  if (!"gene_id" %in% colnames(counts_df)) stop("counts_df must include gene_id", call. = FALSE)
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for differential expression.", call. = FALSE)
  }

  sample_cols <- setdiff(colnames(counts_df), "gene_id")
  design_use <- design_df[as.logical(design_df$include) & design_df$sample_id %in% sample_cols, , drop = FALSE]
  if (nrow(design_use) == 0) return(schema)
  design_use <- design_use[match(design_use$sample_id, sample_cols, nomatch = 0) > 0, , drop = FALSE]
  design_use <- design_use[!duplicated(design_use$sample_id), , drop = FALSE]
  if (nrow(design_use) == 0) return(schema)

  mat <- as.matrix(counts_df[, design_use$sample_id, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat <- round(mat)
  rownames(mat) <- counts_df$gene_id
  keep_rows <- rowSums(mat, na.rm = TRUE) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  if (nrow(mat) == 0) return(schema)

  out <- list()

  # Fit one DESeq2 model per factor referenced in contrasts.
  for (fac in unique(contrast_df$factor)) {
    if (!fac %in% colnames(design_use)) next
    idx <- !is.na(design_use[[fac]]) & nzchar(as.character(design_use[[fac]]))
    coldata <- design_use[idx, c("sample_id", fac), drop = FALSE]
    if (nrow(coldata) < 2) next
    coldata <- coldata[!duplicated(coldata$sample_id), , drop = FALSE]
    cond <- factor(as.character(coldata[[fac]]))
    if (length(levels(cond)) < 2) next

    mat_fac <- mat[, coldata$sample_id, drop = FALSE]
    rownames(coldata) <- coldata$sample_id
    coldata$condition <- cond

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = mat_fac,
      colData = coldata["condition"],
      design = ~ condition
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)

    ctr_rows <- contrast_df[contrast_df$factor == fac, , drop = FALSE]
    for (i in seq_len(nrow(ctr_rows))) {
      ctr <- ctr_rows[i, , drop = FALSE]
      num <- as.character(ctr$numerator)
      den <- as.character(ctr$denominator)
      if (!num %in% levels(coldata$condition) || !den %in% levels(coldata$condition)) next

      res <- DESeq2::results(dds, contrast = c("condition", num, den))
      out[[length(out) + 1]] <- data.frame(
        gene_id = rownames(res),
        contrast_id = ctr$contrast_id,
        log2fc = as.numeric(res$log2FoldChange),
        p_value = as.numeric(res$pvalue),
        fdr = as.numeric(res$padj),
        base_mean = as.numeric(res$baseMean),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) return(schema)
  do.call(rbind, out)
}
