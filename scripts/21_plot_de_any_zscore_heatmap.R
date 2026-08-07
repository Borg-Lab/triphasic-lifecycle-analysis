#!/usr/bin/env Rscript

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

req_profile <- c("gene_id", "tested_any", "de_any")
if (!all(req_profile %in% colnames(profile))) {
  stop("gene_de_dynamic_profile.tsv missing required columns.", call. = FALSE)
}
if (!("gene_id" %in% colnames(tpm))) {
  stop("group_tpm_means.tsv missing gene_id.", call. = FALSE)
}

stage_cols <- grep("^mean_tpm_", colnames(tpm), value = TRUE)
if (length(stage_cols) == 0) {
  stop("No mean_tpm_* columns found in group_tpm_means.tsv.", call. = FALSE)
}
# Drop thallus/thalli stages for this view.
stage_cols <- stage_cols[!grepl("thall(us|i)$", stage_cols, ignore.case = TRUE)]
if (length(stage_cols) == 0) {
  stop("No non-thallus mean_tpm_* columns remain after filtering.", call. = FALSE)
}

de_genes <- profile$gene_id[profile$de_any %in% TRUE]
de_genes <- unique(as.character(de_genes))

sub <- tpm[tpm$gene_id %in% de_genes, c("gene_id", stage_cols), drop = FALSE]
if (nrow(sub) == 0) stop("No DE-at-least-one-contrast genes found in group_tpm_means.tsv.", call. = FALSE)

mat <- as.matrix(sub[, stage_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- sub$gene_id
mat <- log2(mat + 1)

row_mu <- rowMeans(mat, na.rm = TRUE)
row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
zmat <- (mat - row_mu) / row_sd
zmat[!is.finite(zmat)] <- 0

# Clip extremes for robust color mapping.
clip_lim <- 1.5
zmat <- pmax(pmin(zmat, clip_lim), -clip_lim)

stage_names <- sub("^mean_tpm_", "", colnames(zmat))

pal <- grDevices::colorRampPalette(c("#cccccc", "#ffffff", "#af8e67"))(201)
brks <- seq(-clip_lim, clip_lim, length.out = length(pal) + 1)

# Hierarchical clustering for genes and stages
row_dist <- stats::dist(zmat, method = "euclidean")
row_hc <- stats::hclust(row_dist, method = "ward.D2")
col_dist <- stats::dist(t(zmat), method = "euclidean")
col_hc <- stats::hclust(col_dist, method = "average")

# Optional cluster labels for interpretation
k_row <- 8L
row_cluster <- stats::cutree(row_hc, k = k_row)
if (ncol(zmat) >= 6) {
  k_col <- 4L
} else if (ncol(zmat) >= 4) {
  k_col <- 3L
} else {
  k_col <- 2L
}
row_cluster_df <- data.frame(
  gene_id = names(row_cluster),
  de_heatmap_cluster = as.integer(row_cluster),
  stringsAsFactors = FALSE
)
write_tsv(row_cluster_df, file.path(tab_dir, "de_any_zscore_gene_clusters.tsv"))

out_pdf <- file.path(fig_dir, "de_any_zscore_heatmap.pdf")
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Package 'pheatmap' is required for heatmap plotting with visible scale and cluster gaps.", call. = FALSE)
}

# Add row cutree breaks as visual separators.
row_cluster_ord <- row_cluster[row_hc$labels[row_hc$order]]
gaps_row <- which(diff(row_cluster_ord) != 0)
col_cluster <- stats::cutree(col_hc, k = k_col)
col_cluster_ord <- col_cluster[col_hc$labels[col_hc$order]]
gaps_col <- which(diff(col_cluster_ord) != 0)

grDevices::pdf(out_pdf, width = 8.5, height = 10.5, useDingbats = FALSE)
pheatmap::pheatmap(
  mat = zmat,
  color = pal,
  breaks = brks,
  cluster_rows = row_hc,
  cluster_cols = col_hc,
  cutree_rows = k_row,
  cutree_cols = k_col,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  labels_col = stage_names,
  border_color = NA,
  fontsize_col = 9,
  main = sprintf("DE genes in >=1 contrast (n=%d; clustered)", nrow(zmat)),
  legend = TRUE,
  legend_breaks = c(-clip_lim, -1, 0, 1, clip_lim),
  legend_labels = c(sprintf("<=%.1f", -clip_lim), "-1", "0", "1", sprintf(">=%.1f", clip_lim))
)
dev.off()

# ------------------------------------------------------------------
# Second view: developmental trajectory ordering
# ------------------------------------------------------------------
traj_pdf <- file.path(fig_dir, "de_any_zscore_heatmap_trajectory.pdf")

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

traj_stage_order <- c(intersect(dev_stage_order, stage_names), setdiff(stage_names, dev_stage_order))
traj_col_idx <- match(traj_stage_order, stage_names)
traj_col_idx <- traj_col_idx[is.finite(traj_col_idx)]
zmat_traj <- zmat[, traj_col_idx, drop = FALSE]
traj_stage_names <- stage_names[traj_col_idx]

# Order row clusters by where their centroid peaks along developmental flow.
cluster_levels <- sort(unique(row_cluster))
centroids <- vapply(
  cluster_levels,
  function(k) {
    idx <- which(row_cluster == k)
    if (length(idx) == 1) return(zmat_traj[idx, ])
    colMeans(zmat_traj[idx, , drop = FALSE], na.rm = TRUE)
  },
  numeric(ncol(zmat_traj))
)
if (!is.matrix(centroids)) centroids <- matrix(centroids, nrow = ncol(zmat_traj))
centroids <- t(centroids)
peak_idx <- apply(centroids, 1, which.max)
peak_val <- apply(centroids, 1, max)
ord_clusters <- cluster_levels[order(peak_idx, -peak_val)]

# Within each cluster, retain local structure from row dendrogram order.
row_ord_base <- row_hc$labels[row_hc$order]
row_blocks <- lapply(ord_clusters, function(k) {
  genes_k <- names(row_cluster)[row_cluster == k]
  row_ord_base[row_ord_base %in% genes_k]
})
row_order_traj <- unlist(row_blocks, use.names = FALSE)
zmat_traj <- zmat_traj[row_order_traj, , drop = FALSE]

block_sizes <- vapply(row_blocks, length, integer(1))
gaps_row_traj <- cumsum(block_sizes)
gaps_row_traj <- gaps_row_traj[gaps_row_traj < nrow(zmat_traj)]

# Add stage-group gaps by top-level class (f/m/t prefix) to emphasize phase shifts.
stage_group <- sub("^Bm3235([fmt]).*$", "\\1", traj_stage_names)
gaps_col_traj <- which(stage_group[-1] != stage_group[-length(stage_group)])

# Omitted non-manuscript PDF export: de_any_zscore_heatmap_trajectory.pdf

traj_cluster_tbl <- data.frame(
  trajectory_cluster_rank = seq_along(ord_clusters),
  de_heatmap_cluster = as.integer(ord_clusters),
  peak_stage = traj_stage_names[peak_idx[match(ord_clusters, cluster_levels)]],
  peak_centroid_z = as.numeric(peak_val[match(ord_clusters, cluster_levels)]),
  n_genes = as.integer(block_sizes),
  stringsAsFactors = FALSE
)
write_tsv(traj_cluster_tbl, file.path(tab_dir, "de_any_zscore_trajectory_cluster_order.tsv"))

summary_tbl <- data.frame(
  metric = c("n_de_genes_in_heatmap", "n_stages", "z_clip", "row_clusters_k", "col_clusters_k"),
  value = c(nrow(zmat), ncol(zmat), clip_lim, k_row, k_col),
  stringsAsFactors = FALSE
)
write_tsv(summary_tbl, file.path(tab_dir, "de_any_zscore_heatmap_summary.tsv"))

cat("Wrote:\n")
cat(" - ", out_pdf, "\n", sep = "")
cat(" - ", file.path(tab_dir, "de_any_zscore_heatmap_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(tab_dir, "de_any_zscore_gene_clusters.tsv"), "\n", sep = "")
cat(" - ", file.path(tab_dir, "de_any_zscore_trajectory_cluster_order.tsv"), "\n", sep = "")
