#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))

design <- design[as.logical(design$include), , drop = FALSE]
expr <- prepare_expression_inputs(paths, design, analysis)
mods <- safe_read_table(file.path(paths$output_root, "modules", "module_membership.tsv"))
label_map <- safe_read_table(file.path(paths$output_root, "modules", "module_label_map.tsv"))

if (nrow(mods) == 0) stop("module_membership.tsv is empty", call. = FALSE)
if (!all(c("gene_id", "module_id") %in% colnames(mods))) {
  stop("module_membership.tsv must contain gene_id and module_id columns", call. = FALSE)
}
counts <- expr$counts
if (nrow(counts) == 0) stop("No counts available", call. = FALSE)

samples <- intersect(design$sample_id, colnames(counts))
if (length(samples) < 4) stop("Need >=4 samples for module-stage enrichment", call. = FALSE)

counts <- counts[, c("gene_id", samples), drop = FALSE]
genes_use <- intersect(counts$gene_id, mods$gene_id)
counts <- counts[counts$gene_id %in% genes_use, , drop = FALSE]
mods <- mods[mods$gene_id %in% genes_use, , drop = FALSE]

out_fig <- file.path(paths$output_root, "figures")
out_mod <- file.path(paths$output_root, "modules")
ensure_dir(out_fig)
ensure_dir(out_mod)

if (nrow(label_map) > 0 && all(c("module_id", "module_label") %in% colnames(label_map))) {
  module_label_lookup <- setNames(label_map$module_label, label_map$module_id)
} else {
  module_label_lookup <- setNames(character(0), character(0))
}
module_label_safe <- function(x) {
  lab <- unname(module_label_lookup[x])
  lab[is.na(lab) | !nzchar(lab)] <- x[is.na(lab) | !nzchar(lab)]
  lab
}

mat <- as.matrix(counts[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- counts$gene_id
mat <- round(mat)
mat <- log2(mat + 1)

# samples x genes
DatExpr <- t(mat)

# Align module colors to DatExpr columns (genes)
mods_u <- mods[!duplicated(mods$gene_id), c("gene_id", "module_id")]
mods_u <- mods_u[match(colnames(DatExpr), mods_u$gene_id), , drop = FALSE]
colors <- sub("^module_", "", mods_u$module_id)

if (requireNamespace("WGCNA", quietly = TRUE)) {
  MEs <- WGCNA::moduleEigengenes(DatExpr, colors = colors)$eigengenes
  MEs <- WGCNA::orderMEs(MEs)
  me_names <- sub("^ME", "module_", colnames(MEs))
  colnames(MEs) <- me_names
} else {
  # Fallback: mean expression per module
  mod_levels <- sort(unique(mods_u$module_id))
  MEs <- sapply(mod_levels, function(m) {
    g <- mods_u$gene_id[mods_u$module_id == m]
    if (length(g) == 1) as.numeric(DatExpr[, g]) else rowMeans(DatExpr[, g, drop = FALSE], na.rm = TRUE)
  })
  if (is.vector(MEs)) MEs <- matrix(MEs, ncol = 1, dimnames = list(rownames(DatExpr), mod_levels))
  colnames(MEs) <- mod_levels
}

ME_df <- data.frame(sample_id = rownames(MEs), MEs, check.names = FALSE, stringsAsFactors = FALSE)
ME_df <- merge(ME_df, design[, c("sample_id", "stage")], by = "sample_id", all.x = TRUE)
ME_df <- ME_df[!is.na(ME_df$stage), , drop = FALSE]

stage_order <- c(
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

present_stages <- unique(ME_df$stage)
stage_levels <- c(intersect(stage_order, present_stages), setdiff(sort(present_stages), stage_order))

module_cols <- setdiff(colnames(ME_df), c("sample_id", "stage"))

# Stage means per module
agg <- aggregate(ME_df[, module_cols, drop = FALSE], by = list(stage = ME_df$stage), FUN = mean, na.rm = TRUE)
agg$stage <- factor(agg$stage, levels = stage_levels)
agg <- agg[order(agg$stage), , drop = FALSE]

# z-score enrichment by module across stages
zmat <- as.matrix(agg[, module_cols, drop = FALSE])
rownames(zmat) <- as.character(agg$stage)
z_by_mod <- apply(zmat, 2, function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
})
if (is.vector(z_by_mod)) z_by_mod <- matrix(z_by_mod, ncol = 1, dimnames = list(rownames(zmat), module_cols))

# matrix for plotting: modules x stages
z_mod_stage <- t(z_by_mod)

write_tsv(
  data.frame(module_id = rownames(z_mod_stage), z_mod_stage, check.names = FALSE, stringsAsFactors = FALSE),
  file.path(out_mod, "module_stage_enrichment_zscore_matrix.tsv")
)
write_tsv(
  data.frame(stage = rownames(zmat), zmat, check.names = FALSE, stringsAsFactors = FALSE),
  file.path(out_mod, "module_stage_mean_eigengenes.tsv")
)

# Module dendrogram
hc <- stats::hclust(stats::dist(z_mod_stage, method = "euclidean"), method = "ward.D2")
hc_plot <- hc

# Optional clade-flip rotation to the readable module rank order while
# preserving tree topology.
if (nrow(label_map) > 0 && all(c("module_id", "module_rank") %in% colnames(label_map))) {
  rank_num <- suppressWarnings(as.integer(sub("^M", "", label_map$module_rank)))
  desired_module_id_order <- label_map$module_id[order(rank_num, label_map$module_id)]
} else {
  desired_module_id_order <- rownames(z_mod_stage)
}

module_rank_lookup <- setNames(rep(NA_character_, nrow(label_map)), character(0))
if (nrow(label_map) > 0 && all(c("module_id", "module_rank") %in% colnames(label_map))) {
  module_rank_lookup <- setNames(as.character(label_map$module_rank), label_map$module_id)
} else {
  # Fallback from module_label prefix if module_rank is unavailable.
  all_mod_ids <- rownames(z_mod_stage)
  inferred_rank <- sub("^((M[0-9]{2})).*$", "\\1", module_label_safe(all_mod_ids))
  inferred_rank[!grepl("^M[0-9]{2}$", inferred_rank)] <- NA_character_
  module_rank_lookup <- setNames(inferred_rank, all_mod_ids)
}

present_ids <- rownames(z_mod_stage)
desired_ids <- intersect(desired_module_id_order, present_ids)
desired_ids <- unique(desired_ids)

if (length(desired_ids) >= 2) {
  rotated_ok <- FALSE
  if (requireNamespace("ape", quietly = TRUE)) {
    phy <- ape::as.phylo(hc)
    phy_rot <- try(ape::rotateConstr(phy, desired_ids), silent = TRUE)
    if (!inherits(phy_rot, "try-error")) {
      hc_plot <- stats::as.hclust(phy_rot)
      rotated_ok <- TRUE
    }
  }
  if (!rotated_ok) {
    dend <- as.dendrogram(hc)
    if (requireNamespace("dendextend", quietly = TRUE)) {
      dend <- dendextend::rotate(dend, order = desired_ids)
    } else {
      # Best-effort fallback: reorder leaves by target rank weights.
      w <- setNames(seq_along(desired_ids), desired_ids)
      dend <- stats::reorder(dend, wts = w, agglo.FUN = mean)
    }
    hc_plot <- as.hclust(dend)
  }
}

# Keep module_id order from rotated tree before relabeling.
row_order_ids <- hc_plot$labels[hc_plot$order]
row_order_ranks <- unname(module_rank_lookup[row_order_ids])
desired_pos_lookup <- setNames(seq_along(desired_ids), desired_ids)
desired_position <- unname(desired_pos_lookup[row_order_ids])

write_tsv(
  data.frame(
    module_id = row_order_ids,
    module_rank = row_order_ranks,
    desired_position = desired_position,
    matches_desired = row_order_ids == desired_ids[seq_along(row_order_ids)],
    module_label = module_label_safe(row_order_ids),
    stringsAsFactors = FALSE
  ),
  file.path(out_mod, "module_stage_enrichment_row_order.tsv")
)

plot_labels <- module_label_safe(hc_plot$labels)
hc_plot$labels <- plot_labels

# Manuscript workflow uses the z-score matrix and row order downstream; exploratory
# dendrogram/heatmap/lineplot outputs from this step are intentionally omitted.

cat("Wrote:\n")
cat(" - ", file.path(out_mod, "module_stage_enrichment_zscore_matrix.tsv"), "\n", sep = "")
cat(" - ", file.path(out_mod, "module_stage_mean_eigengenes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_mod, "module_stage_enrichment_row_order.tsv"), "\n", sep = "")
