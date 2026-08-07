#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
design <- read_design_table("config/design.tsv", include_all = TRUE)

out_fig <- file.path(paths$output_root, "figures")
out_tab <- file.path(paths$output_root, "tables")
ensure_dir(out_fig)
ensure_dir(out_tab)

tpm <- safe_read_table(paths$tpm_tsv)
counts <- safe_read_table(paths$counts_tsv)

if (nrow(tpm) == 0 || nrow(counts) == 0) stop("Missing tpm/count matrices in paths.yml", call. = FALSE)

colnames(tpm)[1] <- "gene_id"
colnames(counts)[1] <- "gene_id"
tpm$gene_id <- normalize_gene_id(tpm$gene_id)
counts$gene_id <- normalize_gene_id(counts$gene_id)

# Keep all design samples here (including include==FALSE thallus), because this is a prep-method comparison.
samples <- intersect(design$sample_id, intersect(colnames(tpm), colnames(counts)))
if (length(samples) < 2) stop("Not enough overlapping samples in design and matrices.", call. = FALSE)

meta <- design[match(samples, design$sample_id), , drop = FALSE]
meta$prep_type <- ifelse(meta$tissue == "thallus", "NEBNext_thallus", "SmartSeq3_nonthallus")

if (length(unique(meta$prep_type)) < 2) {
  stop("Need both thallus and non-thallus samples for SmartSeq3 vs NEBNext comparison.", call. = FALSE)
}

tpm_m <- as.matrix(tpm[, samples, drop = FALSE]); storage.mode(tpm_m) <- "numeric"
count_m <- as.matrix(counts[, samples, drop = FALSE]); storage.mode(count_m) <- "numeric"
rownames(tpm_m) <- tpm$gene_id
rownames(count_m) <- counts$gene_id

detect_thresholds <- c(0.1, 1, 3)
det_list <- lapply(seq_along(samples), function(i) {
  s <- samples[i]
  data.frame(
    sample_id = s,
    prep_type = meta$prep_type[i],
    stage = meta$stage[i],
    threshold = detect_thresholds,
    n_detected = vapply(detect_thresholds, function(th) sum(tpm_m[, s] > th, na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  )
})
det_df <- do.call(rbind, det_list)

write_tsv(det_df, file.path(out_tab, "smartseq3_vs_nebnext_detectability_by_sample.tsv"))

p_det <- ggplot(det_df, aes(x = prep_type, y = n_detected, fill = prep_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.9, alpha = 0.8) +
  facet_wrap(~threshold, scales = "free_y", nrow = 1, labeller = label_bquote(TPM > .(threshold))) +
  scale_fill_manual(values = c("SmartSeq3_nonthallus" = "#7b3294", "NEBNext_thallus" = "#af8e67")) +
  labs(
    title = "Detected genes by library prep",
    subtitle = "Thallus (NEBNext) vs non-thallus (Smart-seq3)",
    x = NULL,
    y = "Number of detected genes",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.spacing = grid::unit(1.1, "lines"),
    aspect.ratio = 1
  )

# Explicit panel sizing for faceted plot (box region, not only canvas).
det_panel_w <- 3.2
det_panel_h <- 3.2
det_npanels <- length(unique(det_df$threshold))
det_pdf_w <- det_panel_w * det_npanels + 2.2
det_pdf_h <- det_panel_h + 2.2
# Omitted non-manuscript PDF export: ggsave(file.path(out_fig, "smartseq3_vs_nebnext_detectability_boxplot.pdf"), p_det, width = det_pdf_w, height = det_pdf_h)

# Detectability drop-off curve over TPM thresholds
thr_grid <- c(0, 0.1, seq(0.25, 10, by = 0.25))
curve_list <- lapply(seq_along(samples), function(i) {
  s <- samples[i]
  nd <- vapply(thr_grid, function(th) sum(tpm_m[, s] > th, na.rm = TRUE), numeric(1))
  data.frame(
    sample_id = s,
    prep_type = meta$prep_type[i],
    threshold = thr_grid,
    n_detected = nd,
    stringsAsFactors = FALSE
  )
})
curve_df <- do.call(rbind, curve_list)
curve_sum <- do.call(rbind, lapply(split(curve_df, interaction(curve_df$prep_type, curve_df$threshold, drop = TRUE)), function(d) {
  data.frame(
    prep_type = d$prep_type[1],
    threshold = d$threshold[1],
    mean_detected = mean(d$n_detected),
    sd_detected = sd(d$n_detected),
    stringsAsFactors = FALSE
  )
}))

write_tsv(curve_sum, file.path(out_tab, "smartseq3_vs_nebnext_detectability_curve_summary.tsv"))

p_curve <- ggplot(curve_sum, aes(x = threshold, y = mean_detected, color = prep_type, fill = prep_type)) +
  geom_ribbon(aes(ymin = pmax(mean_detected - sd_detected, 0), ymax = mean_detected + sd_detected), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("SmartSeq3_nonthallus" = "#7b3294", "NEBNext_thallus" = "#af8e67")) +
  scale_fill_manual(values = c("SmartSeq3_nonthallus" = "#7b3294", "NEBNext_thallus" = "#af8e67")) +
  labs(
    title = "Detectability drop-off by library prep",
    x = "TPM threshold",
    y = "Mean detected genes",
    color = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(aspect.ratio = 0.75)

curve_panel_w <- 6.6
curve_panel_h <- 5.0
ggsave(file.path(out_fig, "smartseq3_vs_nebnext_detectability_curves.pdf"), p_curve, width = curve_panel_w + 2.0, height = curve_panel_h + 2.0)

# Expression floor comparison using per-sample median TPM among detected genes (>0)
floor_df <- do.call(rbind, lapply(seq_along(samples), function(i) {
  s <- samples[i]
  v <- tpm_m[, s]
  det <- v[v > 0 & is.finite(v)]
  data.frame(
    sample_id = s,
    prep_type = meta$prep_type[i],
    stage = meta$stage[i],
    median_tpm_detected = if (length(det) > 0) median(det) else NA_real_,
    q10_tpm_detected = if (length(det) > 0) as.numeric(quantile(det, 0.1)) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write_tsv(floor_df, file.path(out_tab, "smartseq3_vs_nebnext_expression_floor.tsv"))

p_floor <- ggplot(floor_df, aes(x = prep_type, y = median_tpm_detected, fill = prep_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 2.0, alpha = 0.8) +
  scale_fill_manual(values = c("SmartSeq3_nonthallus" = "#7b3294", "NEBNext_thallus" = "#af8e67")) +
  labs(
    title = "Expression floor comparison",
    subtitle = "Per-sample median TPM among detected genes (TPM > 0)",
    x = NULL,
    y = "Median TPM (detected genes)",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    aspect.ratio = 0.95
  )

floor_panel_w <- 5.8
floor_panel_h <- 5.5
ggsave(file.path(out_fig, "smartseq3_vs_nebnext_expression_floor_boxplot.pdf"), p_floor, width = floor_panel_w + 2.0, height = floor_panel_h + 2.0)

summary_df <- do.call(rbind, lapply(split(det_df, list(det_df$prep_type, det_df$threshold), drop = TRUE), function(d) {
  data.frame(
    prep_type = d$prep_type[1],
    threshold = d$threshold[1],
    mean_detected = mean(d$n_detected),
    sd_detected = sd(d$n_detected),
    n_samples = length(unique(d$sample_id)),
    stringsAsFactors = FALSE
  )
}))
write_tsv(summary_df, file.path(out_tab, "smartseq3_vs_nebnext_detectability_summary.tsv"))

cat("Wrote:\n")
cat(" - ", file.path(out_fig, "smartseq3_vs_nebnext_detectability_boxplot.pdf"), "\n", sep = "")
cat(" - ", file.path(out_fig, "smartseq3_vs_nebnext_detectability_curves.pdf"), "\n", sep = "")
cat(" - ", file.path(out_fig, "smartseq3_vs_nebnext_expression_floor_boxplot.pdf"), "\n", sep = "")
cat(" - ", file.path(out_tab, "smartseq3_vs_nebnext_detectability_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "smartseq3_vs_nebnext_detectability_by_sample.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "smartseq3_vs_nebnext_detectability_curve_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "smartseq3_vs_nebnext_expression_floor.tsv"), "\n", sep = "")
