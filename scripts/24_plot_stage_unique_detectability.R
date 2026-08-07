#!/usr/bin/env Rscript

source("R/io.R")

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
})

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
in_file <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")
out_fig <- file.path(paths$output_root, "figures")
out_tab <- file.path(paths$output_root, "tables")
ensure_dir(out_fig)
ensure_dir(out_tab)

if (!file.exists(in_file)) stop(sprintf("Missing input: %s", in_file), call. = FALSE)
dat <- read.delim(in_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(dat) == 0) stop("group_tpm_means.tsv is empty.", call. = FALSE)

stage_cols <- grep("^mean_tpm_", colnames(dat), value = TRUE)
if (length(stage_cols) == 0) stop("No mean_tpm_ columns found.", call. = FALSE)

keep <- vapply(stage_cols, function(cc) any(is.finite(as.numeric(dat[[cc]]))), logical(1))
stage_cols <- stage_cols[keep]

pref <- c(
  "mean_tpm_Bm3235f_vegtissue",
  "mean_tpm_Bm3235f_procarps",
  "mean_tpm_Bm3235m_vegtissue",
  "mean_tpm_Bm3235m_spermatangia",
  "mean_tpm_Bm3235m_spermatia",
  "mean_tpm_Bm3235f_cystocarps_imm",
  "mean_tpm_Bm3235f_cystocarps_mat",
  "mean_tpm_Bm3235t_carpospores",
  "mean_tpm_Bm3235t_vegtissue",
  "mean_tpm_Bm3235t_stichidia",
  "mean_tpm_Bm3235t_tetraspores"
)
stage_cols <- c(intersect(pref, stage_cols), setdiff(stage_cols, pref))

mat <- as.matrix(dat[, stage_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
n_genes <- nrow(mat)

thresh <- c(1, 3, 5)
det_list <- lapply(thresh, function(t) colSums(mat > t, na.rm = TRUE))
names(det_list) <- paste0("tpm_gt_", thresh)

# 1) Percentage detected per stage
pct_df <- data.frame(
  stage = gsub("^mean_tpm_", "", stage_cols),
  pct_tpm_gt_1 = 100 * det_list$tpm_gt_1[stage_cols] / n_genes,
  pct_tpm_gt_3 = 100 * det_list$tpm_gt_3[stage_cols] / n_genes,
  pct_tpm_gt_5 = 100 * det_list$tpm_gt_5[stage_cols] / n_genes,
  stringsAsFactors = FALSE
)
write_tsv(pct_df, file.path(out_tab, "expr_detectability_percent_summary.tsv"))

pct_long <- rbind(
  data.frame(stage = pct_df$stage, threshold = "TPM > 1", percent = pct_df$pct_tpm_gt_1),
  data.frame(stage = pct_df$stage, threshold = "TPM > 3", percent = pct_df$pct_tpm_gt_3),
  data.frame(stage = pct_df$stage, threshold = "TPM > 5", percent = pct_df$pct_tpm_gt_5)
)
pct_long$stage <- factor(pct_long$stage, levels = pct_df$stage)
pct_long$threshold <- factor(pct_long$threshold, levels = c("TPM > 1", "TPM > 3", "TPM > 5"))

p_pct <- ggplot2::ggplot(pct_long, ggplot2::aes(x = stage, y = percent, color = threshold, group = threshold)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_color_manual(values = c("TPM > 1" = "gray40", "TPM > 3" = "black", "TPM > 5" = "#a50f15")) +
  ggplot2::labs(
    title = "Detectability Percent by Stage",
    x = "Stage",
    y = "Detected genes (%)",
    color = "Threshold"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
# Omitted non-manuscript PDF export: ggplot2::ggsave(file.path(out_fig, "expr_detectability_percent_lineplot.pdf"), p_pct, width = 10, height = 5.5)

# 2) Drop-off curve per stage (within-stage threshold progression)
drop_df <- data.frame(
  stage = rep(gsub("^mean_tpm_", "", stage_cols), each = length(thresh)),
  threshold = rep(thresh, times = length(stage_cols)),
  n_detected = as.numeric(c(rbind(
    det_list$tpm_gt_1[stage_cols],
    det_list$tpm_gt_3[stage_cols],
    det_list$tpm_gt_5[stage_cols]
  ))),
  stringsAsFactors = FALSE
)
write_tsv(drop_df, file.path(out_tab, "expr_detectability_dropoff_summary.tsv"))

p_drop <- ggplot2::ggplot(drop_df, ggplot2::aes(x = threshold, y = n_detected / 1000, group = stage, color = stage)) +
  ggplot2::geom_line(alpha = 0.8) +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::labs(
    title = "Detectability Drop-off by Stage",
    subtitle = "How many genes remain detected as TPM threshold increases",
    x = "TPM threshold",
    y = "Detected genes (x1000)"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "right")
# Omitted non-manuscript PDF export: ggplot2::ggsave(file.path(out_fig, "expr_detectability_dropoff_curves.pdf"), p_drop, width = 10, height = 6.5)

# 3) Shared vs unique prevalence (at TPM > 1)
detect_matrix <- mat > 1
detect_n_stages <- rowSums(detect_matrix, na.rm = TRUE)
prev_tbl <- as.data.frame(table(detect_n_stages), stringsAsFactors = FALSE)
colnames(prev_tbl) <- c("n_stages_detected", "n_genes")
prev_tbl$n_stages_detected <- as.integer(as.character(prev_tbl$n_stages_detected))
prev_tbl <- prev_tbl[order(prev_tbl$n_stages_detected), , drop = FALSE]
write_tsv(prev_tbl, file.path(out_tab, "expr_detectability_prevalence_distribution.tsv"))

unique_stage <- rep(NA_character_, nrow(mat))
for (i in seq_len(nrow(mat))) {
  idx <- which(detect_matrix[i, ])
  if (length(idx) == 1) unique_stage[i] <- gsub("^mean_tpm_", "", stage_cols[idx])
}
unique_tbl <- sort(table(unique_stage), decreasing = TRUE)
unique_df <- data.frame(
  stage = names(unique_tbl),
  n_unique_genes_tpm_gt_1 = as.integer(unique_tbl),
  stringsAsFactors = FALSE
)
unique_df <- unique_df[!is.na(unique_df$stage) & unique_df$stage != "NA", , drop = FALSE]
write_tsv(unique_df, file.path(out_tab, "expr_detectability_unique_genes_by_stage.tsv"))

p_prev <- ggplot2::ggplot(prev_tbl, ggplot2::aes(x = n_stages_detected, y = n_genes)) +
  ggplot2::geom_col(fill = "#3182bd") +
  ggplot2::labs(
    title = "Gene Prevalence Across Stages (TPM > 1)",
    x = "Number of stages where gene is detected",
    y = "Number of genes"
  ) +
  ggplot2::theme_bw(base_size = 10)
# Omitted non-manuscript PDF export: ggplot2::ggsave(file.path(out_fig, "expr_detectability_prevalence_histogram.pdf"), p_prev, width = 8.5, height = 5.5)

if (nrow(unique_df) > 0) {
  unique_df$stage <- factor(unique_df$stage, levels = unique_df$stage)
  p_unique <- ggplot2::ggplot(unique_df, ggplot2::aes(x = stage, y = n_unique_genes_tpm_gt_1)) +
    ggplot2::geom_col(fill = "#6a51a3") +
    ggplot2::labs(
      title = "Stage-Unique Genes (TPM > 1)",
      x = "Stage",
      y = "Number of genes detected in only this stage"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(out_fig, "expr_detectability_stage_unique_barplot.pdf"), p_unique, width = 10, height = 5.5)
}

cat("Wrote detectability extras to ", out_fig, " and ", out_tab, "\n", sep = "")
