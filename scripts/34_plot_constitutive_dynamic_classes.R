#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_fig <- file.path(paths$output_root, "figures")
out_tab <- file.path(paths$output_root, "tables")
ensure_dir(out_fig)
ensure_dir(out_tab)

tpm_path <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")
if (!file.exists(tpm_path)) stop(sprintf("Missing: %s", tpm_path), call. = FALSE)

tpm <- read.delim(tpm_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!("gene_id" %in% colnames(tpm))) stop("group_tpm_means.tsv missing gene_id.", call. = FALSE)

stage_cols <- grep("^mean_tpm_", colnames(tpm), value = TRUE)
if (length(stage_cols) == 0) stop("No mean_tpm_* columns found.", call. = FALSE)

mat <- as.matrix(tpm[, stage_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- tpm$gene_id
stage_names <- sub("^mean_tpm_", "", stage_cols)

# core metrics per gene (relaxed thresholds)
det_thr <- 1
n_stages <- ncol(mat)
detected_n <- apply(mat, 1, function(x) sum(x > det_thr, na.rm = TRUE))
prevalence <- detected_n / n_stages
max_tpm <- apply(mat, 1, function(x) max(x, na.rm = TRUE))
min_tpm <- apply(mat, 1, function(x) min(x, na.rm = TRUE))
dyn_range <- log2(max_tpm + 1) - log2(min_tpm + 1)
cv_log <- apply(log2(mat + 1), 1, function(x) {
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / m
})

gdf <- data.frame(
  gene_id = rownames(mat),
  prevalence = prevalence,
  detected_n = detected_n,
  dynamic_range = dyn_range,
  cv_log = cv_log,
  max_tpm = max_tpm,
  min_tpm = min_tpm,
  stringsAsFactors = FALSE
)

# Dominant stage assignment
dom_idx <- apply(mat, 1, which.max)
gdf$dominant_stage <- stage_names[dom_idx]

# practical classes
prevalence_cut <- 0.75
dynamic_cut <- 1.5

gdf$class <- ifelse(
  gdf$detected_n == 0, "Non-expressed",
  ifelse(
    gdf$detected_n == 1, "Stage-restricted",
  ifelse(
    gdf$detected_n >= ceiling(prevalence_cut * n_stages) & gdf$dynamic_range >= dynamic_cut, "Constitutive + dynamic",
    ifelse(
      gdf$detected_n >= ceiling(prevalence_cut * n_stages) & gdf$dynamic_range < dynamic_cut, "Constitutive + stable",
      "Intermediate"
    )
  )
  )
)

write_tsv(gdf, file.path(out_tab, "constitutive_dynamic_gene_metrics.tsv"))

# Panel A: detected stage count vs dynamic range
stage_count_cut <- ceiling(prevalence_cut * n_stages)
pA <- ggplot(gdf, aes(x = detected_n, y = dynamic_range, color = class)) +
  geom_point(
    size = 0.7,
    alpha = 0.35,
    position = position_jitter(width = 0.20, height = 0.03)
  ) +
  geom_vline(xintercept = stage_count_cut, linetype = 2, color = "#8c8c8c") +
  geom_hline(yintercept = dynamic_cut, linetype = 2, color = "#8c8c8c") +
  scale_x_continuous(breaks = seq.int(0, n_stages, by = 1), limits = c(-0.5, n_stages + 0.5)) +
  scale_color_manual(values = c(
    "Non-expressed" = "#bdbdbd",
    "Constitutive + stable" = "#636363",
    "Constitutive + dynamic" = "#af8e67",
    "Intermediate" = "#756bb1",
    "Stage-restricted" = "#238b45"
  )) +
  labs(
    title = "Detected stages vs dynamic range across stages",
    x = sprintf("Number of stages with mean TPM > %.1f", det_thr),
    y = "Dynamic range: log2(max TPM + 1) - log2(min TPM + 1)",
    color = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(out_fig, "constitutive_dynamic_quadrant.pdf"), pA, width = 7.2, height = 5.6)

# Panel B: stacked class decomposition
class_tab <- as.data.frame(table(gdf$class), stringsAsFactors = FALSE)
colnames(class_tab) <- c("class", "n_genes")
class_levels <- c(
  "Constitutive + dynamic",
  "Constitutive + stable",
  "Intermediate",
  "Stage-restricted",
  "Non-expressed"
)
class_tab$class <- factor(class_tab$class, levels = class_levels)
class_tab <- class_tab[order(class_tab$class), , drop = FALSE]
class_tab$pct <- 100 * class_tab$n_genes / sum(class_tab$n_genes)
class_tab$label <- sprintf("%s\nn=%d (%.1f%%)", class_tab$class, class_tab$n_genes, class_tab$pct)

pB <- ggplot(class_tab, aes(x = 1, y = n_genes, fill = class)) +
  geom_col(width = 0.75, color = "white") +
  geom_text(aes(label = sprintf("n=%d\n%.1f%%", n_genes, pct)), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(
    "Non-expressed" = "#bdbdbd",
    "Constitutive + stable" = "#636363",
    "Constitutive + dynamic" = "#af8e67",
    "Intermediate" = "#756bb1",
    "Stage-restricted" = "#238b45"
  ), breaks = class_levels) +
  labs(
    title = "Genome-wide expression behavior classes",
    x = NULL, y = "Number of genes", fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(out_fig, "constitutive_dynamic_class_stackedbar.pdf"), pB, width = 5.5, height = 5.8)

# Panel C: per-stage fold-change against gene baseline
gene_baseline <- rowMeans(log2(mat + 1), na.rm = TRUE)
lfc_mat <- sweep(log2(mat + 1), 1, gene_baseline, "-")

lfc_long <- do.call(rbind, lapply(seq_along(stage_names), function(i) {
  data.frame(
    stage = stage_names[i],
    lfc_vs_gene_mean = lfc_mat[, i],
    stringsAsFactors = FALSE
  )
}))
lfc_long <- lfc_long[is.finite(lfc_long$lfc_vs_gene_mean), , drop = FALSE]

# order by median lfc
med_ord <- sort(tapply(lfc_long$lfc_vs_gene_mean, lfc_long$stage, median, na.rm = TRUE), decreasing = TRUE)
lfc_long$stage <- factor(lfc_long$stage, levels = names(med_ord))

pC <- ggplot(lfc_long, aes(x = stage, y = lfc_vs_gene_mean)) +
  geom_violin(fill = "#af8e67", color = NA, alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = 2, color = "#8c8c8c") +
  labs(
    title = "Stage-wise expression shifts (relative to each gene's mean)",
    x = NULL,
    y = "Centered log-expression: log2(TPM+1) - gene mean"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Fixed y-range for consistent readability across runs.
pC <- pC + coord_cartesian(ylim = c(-5.1, 5.1))

ggsave(file.path(out_fig, "stage_shift_violin_centered_expression.pdf"), pC, width = 9.8, height = 5.8)

write_tsv(class_tab, file.path(out_tab, "constitutive_dynamic_class_summary.tsv"))
write_tsv(
  data.frame(stage = names(med_ord), median_centered_logexpr = as.numeric(med_ord), stringsAsFactors = FALSE),
  file.path(out_tab, "stage_shift_median_centered_expression.tsv")
)
write_tsv(
  data.frame(
    metric = c("detect_threshold_tpm", "prevalence_cut", "dynamic_cut"),
    value = c(det_thr, prevalence_cut, dynamic_cut),
    stringsAsFactors = FALSE
  ),
  file.path(out_tab, "constitutive_dynamic_thresholds_used.tsv")
)

# Panel D: dynamic range by class + statistics
dr_df <- gdf[, c("gene_id", "class", "dynamic_range"), drop = FALSE]
dr_df <- dr_df[is.finite(dr_df$dynamic_range) & !is.na(dr_df$class), , drop = FALSE]
dr_df$class <- factor(dr_df$class, levels = class_levels)
dr_df <- dr_df[!is.na(dr_df$class), , drop = FALSE]

pD <- ggplot(dr_df, aes(x = class, y = dynamic_range, fill = class)) +
  geom_violin(trim = TRUE, alpha = 0.75, color = NA, width = 1.15) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = grDevices::adjustcolor("white", alpha.f = 0.45), color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c(
    "Non-expressed" = "#bdbdbd",
    "Constitutive + stable" = "#636363",
    "Constitutive + dynamic" = "#af8e67",
    "Intermediate" = "#756bb1",
    "Stage-restricted" = "#238b45"
  ), breaks = class_levels) +
  labs(
    title = "Dynamic range by expression behavior class",
    x = NULL,
    y = "Dynamic range: log2(max TPM + 1) - log2(min TPM + 1)",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")

dr_sum <- do.call(rbind, lapply(split(dr_df, dr_df$class), function(d) {
  data.frame(
    class = as.character(d$class[1]),
    n_genes = nrow(d),
    median_dynamic_range = median(d$dynamic_range, na.rm = TRUE),
    q25_dynamic_range = as.numeric(quantile(d$dynamic_range, 0.25, na.rm = TRUE)),
    q75_dynamic_range = as.numeric(quantile(d$dynamic_range, 0.75, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}))
write_tsv(dr_sum, file.path(out_tab, "dynamic_range_by_class_stats.tsv"))

med_map <- setNames(dr_sum$median_dynamic_range, dr_sum$class)
xlabs <- setNames(
  nm = class_levels,
  object = vapply(class_levels, function(cl) {
    if (cl %in% names(med_map)) {
      sprintf("%s\nmed=%.2f", cl, med_map[[cl]])
    } else {
      sprintf("%s\nmed=NA", cl)
    }
  }, character(1))
)

cliffs_delta <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if (nx == 0 || ny == 0) return(NA_real_)
  # Mann-Whitney U relation: delta = 2U/(nx*ny) - 1
  r <- rank(c(x, y), ties.method = "average")
  rx <- sum(r[seq_len(nx)])
  U <- rx - nx * (nx + 1) / 2
  2 * U / (nx * ny) - 1
}

cls <- class_levels[class_levels %in% unique(as.character(dr_df$class))]
pairs <- utils::combn(cls, 2, simplify = FALSE)
pw <- do.call(rbind, lapply(pairs, function(pp) {
  x <- dr_df$dynamic_range[dr_df$class == pp[1]]
  y <- dr_df$dynamic_range[dr_df$class == pp[2]]
  wt <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
  data.frame(
    class_a = pp[1],
    class_b = pp[2],
    n_a = sum(dr_df$class == pp[1]),
    n_b = sum(dr_df$class == pp[2]),
    median_a = median(x, na.rm = TRUE),
    median_b = median(y, na.rm = TRUE),
    wilcox_p = wt$p.value,
    cliffs_delta = cliffs_delta(x, y),
    stringsAsFactors = FALSE
  )
}))
pw$wilcox_fdr_bh <- stats::p.adjust(pw$wilcox_p, method = "BH")
write_tsv(pw, file.path(out_tab, "dynamic_range_pairwise_wilcoxon.tsv"))

get_fdr <- function(a, b, tbl) {
  idx <- which((tbl$class_a == a & tbl$class_b == b) | (tbl$class_a == b & tbl$class_b == a))
  if (length(idx) == 0) return(NA_real_)
  tbl$wilcox_fdr_bh[idx[1]]
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("FDR=NA")
  if (p < 1e-300) return("FDR<1e-300")
  sprintf("FDR=%.2e", p)
}

# Bracket-style pairwise annotations
cls_present <- unique(as.character(dr_df$class))
cmp <- list(
  c("Constitutive + dynamic", "Intermediate"),
  c("Intermediate", "Stage-restricted"),
  c("Stage-restricted", "Non-expressed"),
  c("Constitutive + stable", "Stage-restricted")
)
cmp <- cmp[vapply(cmp, function(z) all(z %in% cls_present), logical(1))]

y_cap <- 15
y_base <- 12.2
y_step <- 0.7
seg_df <- data.frame()
lab_df <- data.frame()
for (i in seq_along(cmp)) {
  a <- cmp[[i]][1]
  b <- cmp[[i]][2]
  xa <- match(a, class_levels)
  xb <- match(b, class_levels)
  if (!is.finite(xa) || !is.finite(xb)) next
  y <- y_base + (i - 1) * y_step
  seg_df <- rbind(
    seg_df,
    data.frame(x = xa, xend = xa, y = y - y_step * 0.25, yend = y, stringsAsFactors = FALSE),
    data.frame(x = xb, xend = xb, y = y - y_step * 0.25, yend = y, stringsAsFactors = FALSE),
    data.frame(x = xa, xend = xb, y = y, yend = y, stringsAsFactors = FALSE)
  )
  lab_df <- rbind(
    lab_df,
    data.frame(
      x = (xa + xb) / 2,
      y = y + y_step * 0.12,
      label = fmt_p(get_fdr(a, b, pw)),
      stringsAsFactors = FALSE
    )
  )
}

pD_annot <- pD +
  geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = lab_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.9, vjust = 0) +
  scale_x_discrete(labels = xlabs) +
  coord_cartesian(ylim = c(min(dr_df$dynamic_range, na.rm = TRUE), y_cap))

ggsave(file.path(out_fig, "dynamic_range_by_class_boxplot.pdf"), pD_annot, width = 4.4, height = 5.8)

cat("Wrote:\n")
cat(" - ", file.path(out_fig, "constitutive_dynamic_quadrant.pdf"), "\n", sep = "")
cat(" - ", file.path(out_fig, "constitutive_dynamic_class_stackedbar.pdf"), "\n", sep = "")
cat(" - ", file.path(out_fig, "stage_shift_violin_centered_expression.pdf"), "\n", sep = "")
cat(" - ", file.path(out_fig, "dynamic_range_by_class_boxplot.pdf"), "\n", sep = "")
cat(" - ", file.path(out_tab, "constitutive_dynamic_gene_metrics.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "constitutive_dynamic_class_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "stage_shift_median_centered_expression.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "constitutive_dynamic_thresholds_used.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "dynamic_range_by_class_stats.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "dynamic_range_pairwise_wilcoxon.tsv"), "\n", sep = "")
