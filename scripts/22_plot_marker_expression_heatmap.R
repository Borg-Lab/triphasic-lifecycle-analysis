#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
fig_dir <- file.path(paths$output_root, "figures")
tab_dir <- file.path(paths$output_root, "tables")
marker_path <- file.path(paths$output_root, "markers", "marker_table.tsv")
means_path <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")
marker_tf_path <- file.path(paths$output_root, "tables", "marker_tfs.tsv")
ensure_dir(fig_dir)
ensure_dir(tab_dir)

if (!file.exists(marker_path)) stop(sprintf("Missing: %s", marker_path), call. = FALSE)
if (!file.exists(means_path)) stop(sprintf("Missing: %s", means_path), call. = FALSE)

mk <- read.delim(marker_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
means <- read.delim(means_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

req_mk <- c("gene_id", "marker_group", "marker_flag")
if (!all(req_mk %in% colnames(mk))) stop("marker_table.tsv missing required columns.", call. = FALSE)
if (!("gene_id" %in% colnames(means))) stop("group_tpm_means.tsv missing gene_id.", call. = FALSE)

stage_cols <- grep("^mean_tpm_", colnames(means), value = TRUE)
if (length(stage_cols) == 0) stop("No mean_tpm_* columns in group_tpm_means.tsv.", call. = FALSE)
stage_cols <- stage_cols[!grepl("thall(us|i)$", stage_cols, ignore.case = TRUE)]
if (length(stage_cols) == 0) stop("No non-thallus mean_tpm_* columns remain.", call. = FALSE)

markers <- mk[mk$marker_flag %in% TRUE, c("gene_id", "marker_group"), drop = FALSE]
markers <- markers[!is.na(markers$marker_group) & nzchar(markers$marker_group), , drop = FALSE]
markers <- markers[!duplicated(markers$gene_id), , drop = FALSE]
if (nrow(markers) == 0) stop("No marker_flag==TRUE genes found.", call. = FALSE)

sub <- means[means$gene_id %in% markers$gene_id, c("gene_id", stage_cols), drop = FALSE]
if (nrow(sub) == 0) stop("No marker genes found in group_tpm_means.tsv.", call. = FALSE)

m <- match(sub$gene_id, markers$gene_id)
sub$marker_group <- markers$marker_group[m]

mat <- as.matrix(sub[, stage_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- sub$gene_id
mat <- log2(mat + 1)

# Row-wise z-score
row_mu <- rowMeans(mat, na.rm = TRUE)
row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
zmat <- (mat - row_mu) / row_sd
zmat[!is.finite(zmat)] <- 0

clip_lim <- 1.5
zmat <- pmax(pmin(zmat, clip_lim), -clip_lim)

stage_names <- sub("^mean_tpm_", "", colnames(zmat))

dev_stage_order <- c(
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
col_order_names <- c(intersect(dev_stage_order, stage_names), setdiff(stage_names, dev_stage_order))
col_idx <- match(col_order_names, stage_names)
col_idx <- col_idx[is.finite(col_idx)]
zmat <- zmat[, col_idx, drop = FALSE]
stage_names <- stage_names[col_idx]

# Order rows by marker group, then by expression at that group's stage.
row_anno <- data.frame(marker_group = sub$marker_group, stringsAsFactors = FALSE)
rownames(row_anno) <- sub$gene_id
row_anno <- row_anno[rownames(zmat), , drop = FALSE]

# Optional TF family annotation for marker genes.
row_anno$TAP_family <- "non_TF"
if (file.exists(marker_tf_path)) {
  tf_tab <- read.delim(marker_tf_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("gene_id", "TAP_family") %in% colnames(tf_tab))) {
    tf_tab <- tf_tab[!is.na(tf_tab$gene_id) & nzchar(tf_tab$gene_id), c("gene_id", "TAP_family"), drop = FALSE]
    tf_tab <- tf_tab[!duplicated(tf_tab$gene_id), , drop = FALSE]
    fam <- tf_tab$TAP_family[match(rownames(row_anno), tf_tab$gene_id)]
    fam[is.na(fam) | !nzchar(fam)] <- "non_TF"
    row_anno$TAP_family <- fam
  }
}

grp_levels <- unique(markers$marker_group)
grp_levels <- c(intersect(dev_stage_order, grp_levels), setdiff(sort(unique(row_anno$marker_group)), dev_stage_order))
row_anno$marker_group <- factor(row_anno$marker_group, levels = grp_levels)
tf_levels <- c("non_TF", sort(setdiff(unique(row_anno$TAP_family), "non_TF")))
row_anno$TAP_family <- factor(row_anno$TAP_family, levels = tf_levels)

ord_rows <- unlist(lapply(levels(row_anno$marker_group), function(g) {
  ids <- rownames(row_anno)[row_anno$marker_group == g]
  if (length(ids) == 0) return(character())
  st_col <- paste0("mean_tpm_", g)
  if (st_col %in% colnames(sub)) {
    vals <- sub[match(ids, sub$gene_id), st_col]
    ids[order(vals, decreasing = TRUE, na.last = TRUE)]
  } else {
    ids
  }
}), use.names = FALSE)
ord_rows <- ord_rows[ord_rows %in% rownames(zmat)]
zmat_ord <- zmat[ord_rows, , drop = FALSE]
row_anno_ord <- row_anno[ord_rows, , drop = FALSE]

# Row gaps by marker group
grp_sizes <- table(row_anno_ord$marker_group)
grp_sizes <- grp_sizes[grp_sizes > 0]
gaps_row <- cumsum(as.integer(grp_sizes))
gaps_row <- gaps_row[gaps_row < nrow(zmat_ord)]

# Color palettes
heat_pal <- grDevices::colorRampPalette(c("#cccccc", "#ffffff", "#af8e67"))(201)
brks <- seq(-clip_lim, clip_lim, length.out = length(heat_pal) + 1)

stage_cols_palette <- c(
  "Bm3235f_cystocarps_imm" = "#B2ABD2",
  "Bm3235f_cystocarps_mat" = "#8073AC",
  "Bm3235f_pericarp"       = "#C77C7C",
  "Bm3235c_carposporophyte" = "#3F1F78",
  "Bm3235t_carpospores"    = "#5E3C99",
  "Bm3235f_vegtissue"      = "#B55454",
  "Bm3235f_procarps"       = "#7F1D1D",
  "Bm3235t_vegtissue"      = "#C7DFA3",
  "Bm3235t_stichidia"      = "#A8C77A",
  "Bm3235t_tetraspores"    = "#7FA35A",
  "Bm3235m_vegtissue"      = "#7FB8B2",
  "Bm3235m_spermatangia"   = "#3F9188",
  "Bm3235m_spermatia"      = "#0F766E"
)
col_anno <- data.frame(stage = factor(stage_names, levels = stage_names), stringsAsFactors = FALSE)
rownames(col_anno) <- stage_names
ann_colors <- list(
  stage = stage_cols_palette[intersect(names(stage_cols_palette), stage_names)],
  marker_group = stage_cols_palette[intersect(names(stage_cols_palette), as.character(levels(row_anno_ord$marker_group)))]
)

# Ensure all marker groups shown in annotation colors.
mg_levels <- as.character(levels(row_anno_ord$marker_group))
missing_mg <- setdiff(mg_levels, names(ann_colors$marker_group))
if (length(missing_mg) > 0) {
  extra <- grDevices::hcl.colors(length(missing_mg), palette = "Dark 3")
  names(extra) <- missing_mg
  ann_colors$marker_group <- c(ann_colors$marker_group, extra)
}

# TAP family colors
tf_levels_ord <- as.character(levels(row_anno_ord$TAP_family))
tf_pal <- setNames(rep("#d9d9d9", length(tf_levels_ord)), tf_levels_ord)
if (length(tf_levels_ord) > 1) {
  tf_fams <- setdiff(tf_levels_ord, "non_TF")
  fam_cols <- grDevices::hcl.colors(length(tf_fams), palette = "Set 2")
  names(fam_cols) <- tf_fams
  tf_pal[names(fam_cols)] <- fam_cols
}
ann_colors$TAP_family <- tf_pal

out_pdf <- file.path(fig_dir, "marker_genes_expression_heatmap.pdf")

stage_col_use <- unname(stage_cols_palette[stage_names])
stage_col_use[is.na(stage_col_use)] <- "#bdbdbd"
names(stage_col_use) <- stage_names

ha_top <- HeatmapAnnotation(
  stage = stage_names,
  col = list(stage = stage_col_use),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)

ha_left <- rowAnnotation(
  marker_group = row_anno_ord$marker_group,
  col = list(marker_group = ann_colors$marker_group),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)

tf_idx <- which(as.character(row_anno_ord$TAP_family) != "non_TF")
if (length(tf_idx) > 0) {
  ha_right <- rowAnnotation(
    TAP_family = anno_mark(
      at = tf_idx,
      labels = as.character(row_anno_ord$TAP_family[tf_idx]),
      side = "right",
      labels_gp = gpar(fontsize = 6),
      link_width = unit(4, "mm"),
      padding = unit(1, "mm")
    ),
    width = unit(35, "mm")
  )
} else {
  ha_right <- NULL
}

col_fun <- circlize::colorRamp2(c(-clip_lim, 0, clip_lim), c("#cccccc", "#ffffff", "#663333"))

ht <- Heatmap(
  zmat_ord,
  name = "z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = stage_names,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  top_annotation = ha_top,
  left_annotation = ha_left,
  right_annotation = ha_right,
  row_split = row_anno_ord$marker_group,
  row_title = NULL,
  border = FALSE,
  column_title = sprintf("Marker genes across stages (n=%d)", nrow(zmat_ord)),
  heatmap_legend_param = list(
    at = c(-clip_lim, -1, 0, 1, clip_lim),
    labels = c(sprintf("<=%.1f", -clip_lim), "-1", "0", "1", sprintf(">=%.1f", clip_lim))
  )
)

grDevices::pdf(out_pdf, width = 10.5, height = 11, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

summary_tbl <- data.frame(
  metric = c("n_marker_genes", "n_stages", "z_clip"),
  value = c(nrow(zmat_ord), ncol(zmat_ord), clip_lim),
  stringsAsFactors = FALSE
)
write_tsv(summary_tbl, file.path(tab_dir, "marker_genes_expression_heatmap_summary.tsv"))
write_tsv(data.frame(
  gene_id = rownames(zmat_ord),
  marker_group = as.character(row_anno_ord$marker_group),
  TAP_family = as.character(row_anno_ord$TAP_family),
  stringsAsFactors = FALSE
),
          file.path(tab_dir, "marker_genes_heatmap_row_order.tsv"))

cat("Wrote:\n")
cat(" - ", out_pdf, "\n", sep = "")
cat(" - ", file.path(tab_dir, "marker_genes_expression_heatmap_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(tab_dir, "marker_genes_heatmap_row_order.tsv"), "\n", sep = "")
