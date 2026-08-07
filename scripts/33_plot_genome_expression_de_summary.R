#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
fig_dir <- file.path(paths$output_root, "figures")
tab_dir <- file.path(paths$output_root, "tables")
ensure_dir(fig_dir)

profile_path <- file.path(tab_dir, "gene_de_dynamic_profile.tsv")
tpm_path <- file.path(tab_dir, "group_tpm_means.tsv")

if (!file.exists(profile_path)) stop(sprintf("Missing: %s", profile_path), call. = FALSE)
if (!file.exists(tpm_path)) stop(sprintf("Missing: %s", tpm_path), call. = FALSE)

profile <- read.delim(profile_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tpm <- read.delim(tpm_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("gene_id", "tested_any", "de_any", "class", "log2_range", "constitutive_tpm1_all_stages", "expressed_any_tpm1") %in% colnames(profile))) {
  stop("gene_de_dynamic_profile.tsv missing required columns.", call. = FALSE)
}

# ------------------------------------------------------------------
# 1) Violin/box: dynamic range DE vs stable (tested only)
# ------------------------------------------------------------------
d_test <- profile[profile$tested_any, c("gene_id", "class", "log2_range"), drop = FALSE]
d_test$class <- factor(d_test$class, levels = c("stable_nonDE", "DE_at_least_one_contrast"))

p_violin <- ggplot2::ggplot(d_test, ggplot2::aes(x = class, y = log2_range, fill = class)) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
  ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = c("stable_nonDE" = "#4d4d4d", "DE_at_least_one_contrast" = "#b2182b")) +
  ggplot2::labs(
    title = "Dynamic range across stages (tested genes only)",
    x = NULL,
    y = "log2(max TPM + 1) - log2(min TPM + 1)"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")

# Omitted non-manuscript PDF export: ggplot2::ggsave(file.path(fig_dir, "de_dynamic_range_violin_tested_only.pdf"), p_violin, width = 6.5, height = 5.0)

# ------------------------------------------------------------------
# 2) Quantile ribbons by stage for constitutive genes
# ------------------------------------------------------------------
tpm_cols <- grep("^mean_tpm_", colnames(tpm), value = TRUE)
if (length(tpm_cols) == 0) stop("No mean_tpm_* columns in group_tpm_means.tsv", call. = FALSE)

stage_names <- sub("^mean_tpm_", "", tpm_cols)
mat <- as.matrix(tpm[, tpm_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- tpm$gene_id

# Keep only stages with values.
keep <- colSums(!is.na(mat)) > 0
mat <- mat[, keep, drop = FALSE]
stage_names <- stage_names[keep]
tpm_cols <- tpm_cols[keep]

pref_stage_order <- c(
  "Bm3235f_vegtissue",
  "Bm3235f_procarps",
  "Bm3235m_vegtissue",
  "Bm3235m_spermatangia",
  "Bm3235m_spermatia",
  "Bm3235f_pericarp",
  "Bm3235f_cystocarps_imm",
  "Bm3235f_cystocarps_mat",
  "Bm3235c_carposporophyte",
  "Bm3235t_carpospores",
  "Bm3235t_vegtissue",
  "Bm3235t_stichidia",
  "Bm3235t_tetraspores"
)
stage_levels <- c(intersect(pref_stage_order, stage_names), setdiff(stage_names, pref_stage_order))
ord <- match(stage_levels, stage_names)
ord <- ord[is.finite(ord)]
mat <- mat[, ord, drop = FALSE]
stage_names <- stage_names[ord]

const_stable <- profile$gene_id[profile$constitutive_tpm1_all_stages & profile$tested_any & !profile$de_any]
const_de <- profile$gene_id[profile$constitutive_tpm1_all_stages & profile$tested_any & profile$de_any]

calc_quantiles <- function(gene_ids, label) {
  gids <- intersect(gene_ids, rownames(mat))
  if (length(gids) == 0) return(NULL)
  subm <- log2(mat[gids, , drop = FALSE] + 1)
  out <- data.frame(
    stage = stage_names,
    q10 = apply(subm, 2, function(x) as.numeric(stats::quantile(x, 0.10, na.rm = TRUE))),
    q50 = apply(subm, 2, function(x) as.numeric(stats::quantile(x, 0.50, na.rm = TRUE))),
    q90 = apply(subm, 2, function(x) as.numeric(stats::quantile(x, 0.90, na.rm = TRUE))),
    n_genes = length(gids),
    group = label,
    stringsAsFactors = FALSE
  )
  out
}

q_stable <- calc_quantiles(const_stable, "constitutive_stable_nonDE")
q_de <- calc_quantiles(const_de, "constitutive_DE_any")
qdf <- rbind(q_stable, q_de)
qdf$stage <- factor(qdf$stage, levels = stage_levels)
qdf$group <- factor(qdf$group, levels = c("constitutive_stable_nonDE", "constitutive_DE_any"))

p_ribbon <- ggplot2::ggplot(qdf, ggplot2::aes(x = stage, y = q50, color = group, fill = group, group = group)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q10, ymax = q90), alpha = 0.18, color = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::scale_color_manual(values = c("constitutive_stable_nonDE" = "#4d4d4d", "constitutive_DE_any" = "#b2182b")) +
  ggplot2::scale_fill_manual(values = c("constitutive_stable_nonDE" = "#4d4d4d", "constitutive_DE_any" = "#b2182b")) +
  ggplot2::labs(
    title = "Constitutive genes: stage-wise quantiles",
    subtitle = "Median line with 10-90% ribbon of log2(TPM+1)",
    x = "Stage",
    y = "log2(TPM + 1)",
    color = NULL,
    fill = NULL
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

# Omitted non-manuscript PDF export: ggplot2::ggsave(file.path(fig_dir, "constitutive_quantile_ribbon_by_stage.pdf"), p_ribbon, width = 10.5, height = 5.5)

write_tsv(qdf, file.path(tab_dir, "constitutive_quantiles_by_stage.tsv"))

# ------------------------------------------------------------------
# 3) Donut summary: genome expression/DE classes
# ------------------------------------------------------------------
donut <- profile[, c("gene_id", "expressed_any_tpm1", "tested_any", "de_any"), drop = FALSE]
donut$segment <- ifelse(
  !donut$expressed_any_tpm1, "Not expressed (TPM<=1)",
  ifelse(donut$tested_any & donut$de_any, "Expressed + DE (>=1 contrast)",
    ifelse(donut$tested_any & !donut$de_any, "Expressed + stable non-DE", "Expressed + not DE-tested")
  )
)

seg_tbl <- as.data.frame(table(donut$segment), stringsAsFactors = FALSE)
colnames(seg_tbl) <- c("segment", "n_genes")
seg_tbl$pct <- 100 * seg_tbl$n_genes / sum(seg_tbl$n_genes)
seg_tbl <- seg_tbl[order(seg_tbl$segment), , drop = FALSE]
seg_tbl$label <- sprintf("%s\nn=%d (%.1f%%)", seg_tbl$segment, seg_tbl$n_genes, seg_tbl$pct)
seg_tbl$ymax <- cumsum(seg_tbl$n_genes)
seg_tbl$ymin <- c(0, head(seg_tbl$ymax, -1))
seg_tbl$ypos <- (seg_tbl$ymin + seg_tbl$ymax) / 2

pal <- c(
  "Expressed + DE (>=1 contrast)" = "#b2182b",
  "Expressed + stable non-DE" = "#4d4d4d",
  "Expressed + not DE-tested" = "#9ecae1",
  "Not expressed (TPM<=1)" = "#d9d9d9"
)
seg_tbl$segment <- factor(seg_tbl$segment, levels = names(pal))

p_donut <- ggplot2::ggplot(seg_tbl, ggplot2::aes(x = 2, y = n_genes, fill = segment)) +
  ggplot2::geom_col(color = "white", width = 1) +
  ggplot2::geom_text(
    ggplot2::aes(y = ypos, label = sprintf("n=%d\n%.1f%%", n_genes, pct)),
    x = 2.05,
    size = 3.0,
    lineheight = 0.95
  ) +
  ggplot2::coord_polar(theta = "y") +
  ggplot2::xlim(0.5, 2.5) +
  ggplot2::scale_fill_manual(values = pal) +
  ggplot2::labs(
    title = "Genome-wide expression and DE summary",
    subtitle = "Partition of genes by expression and DE status",
    fill = NULL
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(legend.position = "right")

# Omitted non-manuscript PDF export: ggplot2::ggsave(file.path(fig_dir, "genome_expression_de_donut.pdf"), p_donut, width = 8.0, height = 5.8)

p_bar <- ggplot2::ggplot(seg_tbl, ggplot2::aes(x = "Genome", y = pct, fill = segment)) +
  ggplot2::geom_col(color = "white", width = 0.48) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%d\n%.1f%%", n_genes, pct)),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.1,
    lineheight = 0.95
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
  ggplot2::scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Genome-wide expression and DE summary",
    subtitle = "Partition of genes by expression and DE status",
    x = NULL,
    y = "Genes (%)",
    fill = NULL
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(file.path(fig_dir, "genome_expression_de_stacked_barplot.pdf"), p_bar, width = 8.0, height = 3.2)
write_tsv(seg_tbl[, c("segment", "n_genes", "pct")], file.path(tab_dir, "genome_expression_de_donut_summary.tsv"))

cat("Wrote:\n")
cat(" - ", file.path(fig_dir, "de_dynamic_range_violin_tested_only.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "constitutive_quantile_ribbon_by_stage.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "genome_expression_de_donut.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "genome_expression_de_stacked_barplot.pdf"), "\n", sep = "")
cat(" - ", file.path(tab_dir, "constitutive_quantiles_by_stage.tsv"), "\n", sep = "")
cat(" - ", file.path(tab_dir, "genome_expression_de_donut_summary.tsv"), "\n", sep = "")
