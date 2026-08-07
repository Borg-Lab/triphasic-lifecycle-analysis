#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
contrasts <- read_contrast_table(Sys.getenv("CONTRASTS_CONFIG", unset = "config/contrasts.tsv"))
expr <- prepare_expression_inputs(paths, design, analysis)

out_dir <- file.path(paths$output_root, "review_robustness")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

stage_order <- c(
  "Bm3235f_vegtissue", "Bm3235f_procarps",
  "Bm3235m_vegtissue", "Bm3235m_spermatangia", "Bm3235m_spermatia",
  "Bm3235f_pericarp", "Bm3235f_cystocarps_imm", "Bm3235f_cystocarps_mat",
  "Bm3235c_carposporophyte", "Bm3235t_carpospores",
  "Bm3235t_vegtissue", "Bm3235t_stichidia", "Bm3235t_tetraspores"
)

counts <- expr$counts
sample_ids <- intersect(design$sample_id[as.logical(design$include)], colnames(counts))
counts <- counts[, c("gene_id", sample_ids), drop = FALSE]
mat <- as.matrix(counts[, sample_ids, drop = FALSE])
storage.mode(mat) <- "numeric"
mat <- round(mat)
rownames(mat) <- counts$gene_id
mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]

design_use <- design[match(colnames(mat), design$sample_id), , drop = FALSE]
rownames(design_use) <- design_use$sample_id

# Stable/low-variance control genes: reproducible from current group TPM means.
tpm <- read.delim(file.path(paths$output_root, "tables", "group_tpm_means.tsv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tpm_cols <- grep("^mean_tpm_", colnames(tpm), value = TRUE)
tpm_mat <- as.matrix(tpm[, tpm_cols, drop = FALSE])
storage.mode(tpm_mat) <- "numeric"
rownames(tpm_mat) <- tpm$gene_id
stage_ids <- sub("^mean_tpm_", "", tpm_cols)

common_tpm <- intersect(rownames(tpm_mat), rownames(mat))
tpm_mat <- tpm_mat[common_tpm, , drop = FALSE]
log_tpm <- log2(tpm_mat + 1)
var_log_tpm <- apply(log_tpm, 1, stats::var, na.rm = TRUE)
iqr_log_tpm <- apply(log_tpm, 1, stats::IQR, na.rm = TRUE)
max_tpm <- apply(tpm_mat, 1, max, na.rm = TRUE)
min_tpm <- apply(tpm_mat, 1, min, na.rm = TRUE)
detected_all <- rowSums(tpm_mat > analysis$tpm_detect_threshold, na.rm = TRUE) == ncol(tpm_mat)

stable_candidates <- data.frame(
  gene_id = rownames(tpm_mat),
  max_tpm = max_tpm,
  min_tpm = min_tpm,
  variance_log2_tpm = var_log_tpm,
  iqr_log2_tpm = iqr_log_tpm,
  detected_all_stages = detected_all,
  stringsAsFactors = FALSE
)
stable_candidates <- stable_candidates[stable_candidates$detected_all_stages & is.finite(stable_candidates$variance_log2_tpm), , drop = FALSE]
stable_candidates <- stable_candidates[order(stable_candidates$variance_log2_tpm, stable_candidates$iqr_log2_tpm, -stable_candidates$min_tpm), , drop = FALSE]

# Match the current run by default: use the current stable non-DE count as the
# target size, then select that many lowest-variance detected genes.
profile_path <- file.path(paths$output_root, "tables", "gene_de_dynamic_profile.tsv")
if (file.exists(profile_path)) {
  prof <- read.delim(profile_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  nonde_stable <- prof$gene_id[isTRUE("tested_any" %in% names(prof)) & prof$tested_any & !prof$de_any]
} else {
  nonde_stable <- character()
}
default_target_n <- if (length(nonde_stable) > 0) length(nonde_stable) else 1157L
stable_target_n <- as.integer(Sys.getenv("DE_ROBUST_STABLE_N", as.character(default_target_n)))
lowvar_target_n <- as.integer(Sys.getenv("DE_ROBUST_LOWVAR_N", as.character(default_target_n)))
stable_target_n <- max(1L, stable_target_n)
lowvar_target_n <- max(1L, lowvar_target_n)

stable_genes <- head(stable_candidates$gene_id, min(stable_target_n, nrow(stable_candidates)))
lowvar_genes <- stable_genes[seq_len(min(lowvar_target_n, length(stable_genes)))]
control_label <- sprintf("lowest_variance_%d_control_size_factors", length(stable_genes))
lowvar_label <- sprintf("lowest_variance_%d_gene_subset", length(lowvar_genes))

stable_table <- stable_candidates
stable_table$selected_as_lowest_variance_control <- stable_table$gene_id %in% stable_genes
stable_table$current_stable_nonDE <- stable_table$gene_id %in% nonde_stable
write_tsv(stable_table, file.path(out_dir, "de_robustness_lowest_variance_gene_candidates.tsv"))

summarise_de <- function(res_df, label) {
  if (nrow(res_df) == 0) return(data.frame())
  ctab <- merge(res_df, contrasts[, c("contrast_id", "lfc_threshold", "fdr_threshold")], by = "contrast_id", all.x = TRUE, sort = FALSE)
  ctab$de_hit <- with(ctab, is.finite(fdr) & is.finite(log2fc) & fdr <= fdr_threshold & abs(log2fc) >= lfc_threshold)
  gene_any <- aggregate(ctab$de_hit, list(gene_id = ctab$gene_id), function(x) any(x, na.rm = TRUE))
  colnames(gene_any)[2] <- "de_any"
  gene_test <- aggregate(is.finite(ctab$fdr), list(gene_id = ctab$gene_id), function(x) any(x, na.rm = TRUE))
  colnames(gene_test)[2] <- "tested_any"
  g <- merge(gene_any, gene_test, by = "gene_id", all = TRUE)
  n_tested <- sum(g$tested_any, na.rm = TRUE)
  n_de <- sum(g$tested_any & g$de_any, na.rm = TRUE)
  by_contrast <- do.call(rbind, lapply(split(ctab, ctab$contrast_id), function(d) {
    data.frame(
      analysis = label,
      contrast_id = unique(d$contrast_id),
      tested_genes = sum(is.finite(d$fdr)),
      up_sig = sum(d$de_hit & d$log2fc > 0, na.rm = TRUE),
      down_sig = sum(d$de_hit & d$log2fc < 0, na.rm = TRUE),
      total_sig = sum(d$de_hit, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  list(
    summary = data.frame(
      analysis = label,
      n_genes_tested = n_tested,
      n_genes_stage_dependent = n_de,
      pct_stage_dependent = 100 * n_de / max(1, n_tested),
      stringsAsFactors = FALSE
    ),
    by_contrast = by_contrast,
    gene_any = transform(g, analysis = label)
  )
}

run_deseq_analysis <- function(label, control_genes = NULL, restrict_genes = NULL) {
  mat_use <- mat
  if (!is.null(restrict_genes)) {
    keep <- intersect(rownames(mat_use), restrict_genes)
    mat_use <- mat_use[keep, , drop = FALSE]
  }
  out <- list()
  norm_counts_list <- list()
  sf_rows <- list()
  for (fac in unique(contrasts$factor)) {
    if (!fac %in% colnames(design_use)) next
    idx <- !is.na(design_use[[fac]]) & nzchar(as.character(design_use[[fac]]))
    coldata <- design_use[idx, c("sample_id", fac), drop = FALSE]
    cond <- factor(as.character(coldata[[fac]]))
    if (length(levels(cond)) < 2) next
    rownames(coldata) <- coldata$sample_id
    coldata$condition <- cond
    mat_fac <- mat_use[, coldata$sample_id, drop = FALSE]
    dds <- DESeqDataSetFromMatrix(countData = mat_fac, colData = coldata["condition"], design = ~ condition)
    if (!is.null(control_genes)) {
      ctrl <- rownames(dds) %in% control_genes
      if (sum(ctrl) < 10) stop("Too few control genes available for size factor estimation.", call. = FALSE)
      dds <- estimateSizeFactors(dds, controlGenes = ctrl)
    }
    dds <- DESeq(dds, quiet = TRUE)
    sf_rows[[fac]] <- data.frame(analysis = label, factor = fac, sample_id = names(sizeFactors(dds)), size_factor = as.numeric(sizeFactors(dds)), stringsAsFactors = FALSE)
    norm_counts_list[[fac]] <- as.data.frame(counts(dds, normalized = TRUE))
    norm_counts_list[[fac]]$gene_id <- rownames(dds)
    ctr_rows <- contrasts[contrasts$factor == fac, , drop = FALSE]
    for (i in seq_len(nrow(ctr_rows))) {
      ctr <- ctr_rows[i, , drop = FALSE]
      num <- as.character(ctr$numerator); den <- as.character(ctr$denominator)
      if (!num %in% levels(coldata$condition) || !den %in% levels(coldata$condition)) next
      res <- results(dds, contrast = c("condition", num, den))
      out[[length(out) + 1L]] <- data.frame(
        analysis = label,
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
  res_df <- if (length(out)) do.call(rbind, out) else data.frame()
  sf_df <- if (length(sf_rows)) do.call(rbind, sf_rows) else data.frame()
  # Keep one normalized matrix for the stage factor, falling back to first available.
  nc <- if ("stage" %in% names(norm_counts_list)) norm_counts_list[["stage"]] else norm_counts_list[[1]]
  list(results = res_df, size_factors = sf_df, norm_counts = nc)
}

calc_stage_metrics <- function(norm_counts, label, gene_subset_label) {
  if (is.null(norm_counts) || nrow(norm_counts) == 0) return(data.frame())
  sample_cols <- intersect(design_use$sample_id, colnames(norm_counts))
  nc <- as.matrix(norm_counts[, sample_cols, drop = FALSE])
  storage.mode(nc) <- "numeric"
  rownames(nc) <- norm_counts$gene_id
  stage <- design_use$stage[match(sample_cols, design_use$sample_id)]
  stage_levels <- c(intersect(stage_order, unique(stage)), setdiff(unique(stage), stage_order))
  stage_mean <- sapply(stage_levels, function(st) rowMeans(nc[, stage == st, drop = FALSE], na.rm = TRUE))
  if (is.null(dim(stage_mean))) stage_mean <- matrix(stage_mean, ncol = 1, dimnames = list(rownames(nc), stage_levels))
  log_stage <- log2(stage_mean + 1)
  dynamic <- apply(log_stage, 1, max, na.rm = TRUE) - apply(log_stage, 1, min, na.rm = TRUE)
  stage_median <- apply(log_stage, 2, median, na.rm = TRUE)
  all_except <- sapply(stage_levels, function(st) median(log_stage[, st] - rowMeans(log_stage[, setdiff(stage_levels, st), drop = FALSE], na.rm = TRUE), na.rm = TRUE))
  data.frame(
    analysis = label,
    gene_subset = gene_subset_label,
    n_genes = nrow(log_stage),
    median_dynamic_range_log2_norm_counts = median(dynamic, na.rm = TRUE),
    q25_dynamic_range_log2_norm_counts = as.numeric(quantile(dynamic, 0.25, na.rm = TRUE)),
    q75_dynamic_range_log2_norm_counts = as.numeric(quantile(dynamic, 0.75, na.rm = TRUE)),
    spermatia_median_log2_norm = unname(stage_median["Bm3235m_spermatia"]),
    spermatia_shift_vs_other_median = unname(all_except["Bm3235m_spermatia"]),
    tetraspores_shift_vs_other_median = unname(all_except["Bm3235t_tetraspores"]),
    carpospores_shift_vs_other_median = unname(all_except["Bm3235t_carpospores"]),
    stringsAsFactors = FALSE
  )
}

analyses <- list(
  standard = run_deseq_analysis("standard_DESeq2_median_ratio"),
  stable_sf = run_deseq_analysis(control_label, control_genes = stable_genes),
  lowvar_only = run_deseq_analysis(lowvar_label, restrict_genes = lowvar_genes)
)

all_de <- do.call(rbind, lapply(analyses, `[[`, "results"))
all_sf <- do.call(rbind, lapply(analyses, `[[`, "size_factors"))
write_tsv(all_de, file.path(out_dir, "de_robustness_de_results.tsv"))
write_tsv(all_sf, file.path(out_dir, "de_robustness_size_factors.tsv"))
write_tsv(data.frame(gene_id = stable_genes, lowest_variance_rank = seq_along(stable_genes), stringsAsFactors = FALSE), file.path(out_dir, "de_robustness_lowest_variance_control_genes.tsv"))

sum_list <- lapply(names(analyses), function(nm) summarise_de(analyses[[nm]]$results, analyses[[nm]]$results$analysis[1]))
summary_df <- do.call(rbind, lapply(sum_list, `[[`, "summary"))
contrast_df <- do.call(rbind, lapply(sum_list, `[[`, "by_contrast"))
gene_any_df <- do.call(rbind, lapply(sum_list, `[[`, "gene_any"))
write_tsv(summary_df, file.path(out_dir, "de_robustness_stage_dependent_summary.tsv"))
write_tsv(contrast_df, file.path(out_dir, "de_robustness_contrast_summary.tsv"))
write_tsv(gene_any_df, file.path(out_dir, "de_robustness_gene_any_de.tsv"))

metrics <- do.call(rbind, list(
  calc_stage_metrics(analyses$standard$norm_counts, "standard_DESeq2_median_ratio", "all_tested_genes"),
  calc_stage_metrics(analyses$stable_sf$norm_counts, control_label, "all_tested_genes"),
  calc_stage_metrics(analyses$lowvar_only$norm_counts, lowvar_label, sprintf("lowest_variance_%d_genes", length(lowvar_genes)))
))
write_tsv(metrics, file.path(out_dir, "de_robustness_dynamic_global_metrics.tsv"))

audit <- data.frame(
  metric = c("default_target_from_current_stable_nonDE", "lowest_variance_target_n", "stable_candidates_detected_all_tpm_gt_threshold", "lowest_variance_control_genes_selected", "current_stable_nonDE_available", "overlap_control_with_current_stable_nonDE", "lowest_variance_subset_n"),
  value = c(default_target_n, stable_target_n, nrow(stable_candidates), length(stable_genes), length(nonde_stable), length(intersect(stable_genes, nonde_stable)), length(lowvar_genes)),
  stringsAsFactors = FALSE
)
write_tsv(audit, file.path(out_dir, "de_robustness_audit.tsv"))

analysis_label_map <- setNames(
  c(
    "Standard DESeq2\nall tested",
    sprintf("Lowest-variance %d\ncontrol size factors\nall tested", length(stable_genes)),
    sprintf("Lowest-variance %d\ngene subset\nn=%d", length(lowvar_genes), length(lowvar_genes))
  ),
  c("standard_DESeq2_median_ratio", control_label, lowvar_label)
)
analysis_colors <- setNames(c("#8E2F3F", "#B8860B", "#4C2A85"), names(analysis_label_map))
summary_df$analysis_label <- factor(analysis_label_map[summary_df$analysis], levels = unname(analysis_label_map))
metrics$analysis_label <- factor(analysis_label_map[metrics$analysis], levels = unname(analysis_label_map))

p <- ggplot(summary_df, aes(x = analysis_label, y = pct_stage_dependent, fill = analysis)) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = sprintf("%.1f%%", pct_stage_dependent)), vjust = -0.25, size = 3) +
  scale_fill_manual(values = analysis_colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, max(summary_df$pct_stage_dependent, na.rm = TRUE) * 1.15)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank()) +
  labs(title = "DE robustness to normalization / low-variance restriction", x = NULL, y = "Stage-dependent genes (% tested)")
ggsave(file.path(fig_dir, "de_robustness_stage_dependent_fraction.pdf"), p, width = 7.4, height = 4.8)

p2 <- ggplot(metrics, aes(x = analysis_label, y = median_dynamic_range_log2_norm_counts, fill = analysis)) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = sprintf("%.2f", median_dynamic_range_log2_norm_counts)), vjust = -0.25, size = 3) +
  scale_fill_manual(values = analysis_colors, drop = FALSE) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank()) +
  labs(title = "Dynamic range robustness", x = NULL, y = "Median log2 normalized-count range")
ggsave(file.path(fig_dir, "de_robustness_dynamic_range_summary.pdf"), p2, width = 7.4, height = 4.8)

cat("Wrote DE robustness outputs to ", out_dir, "\n", sep = "")
