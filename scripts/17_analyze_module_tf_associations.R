#!/usr/bin/env Rscript

source("R/io.R")

is_tf <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(trimws(x)) & trimws(x) != "-" & trimws(x) != "NA"
}

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_dir <- file.path(paths$output_root, "tables")
ensure_dir(out_dir)

marker_path <- file.path(paths$output_root, "markers", "marker_table.tsv")
module_path <- file.path(paths$output_root, "modules", "module_membership.tsv")
label_map_path <- file.path(paths$output_root, "modules", "module_label_map.tsv")
cluster_marker_path <- file.path(paths$output_root, "enrichment", "cluster_marker_enrichment.tsv")

if (!file.exists(marker_path)) stop(sprintf("Missing marker table: %s", marker_path), call. = FALSE)
if (!file.exists(module_path)) stop(sprintf("Missing module membership table: %s", module_path), call. = FALSE)
if (!file.exists(label_map_path)) stop(sprintf("Missing module label map: %s", label_map_path), call. = FALSE)
if (!file.exists(cluster_marker_path)) stop(sprintf("Missing cluster marker enrichment table: %s", cluster_marker_path), call. = FALSE)

marker_df <- read.delim(marker_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
module_df <- read.delim(module_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
label_map <- read.delim(label_map_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cluster_marker_df <- read.delim(cluster_marker_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!("TAP_family" %in% colnames(marker_df))) {
  stop("marker_table.tsv has no TAP_family column.", call. = FALSE)
}
if (!("module_id" %in% colnames(module_df))) {
  stop("module_membership.tsv has no module_id column.", call. = FALSE)
}
if (!("module_id" %in% colnames(label_map))) {
  stop("module_label_map.tsv has no module_id column.", call. = FALSE)
}

label_cols <- c(
  "module_id", "module_rank", "module_label",
  "module_label_generation", "module_label_context",
  "dominant_stage", "dominant_generation",
  "second_stage", "second_generation"
)
label_cols <- label_cols[label_cols %in% colnames(label_map)]
label_map <- unique(label_map[, label_cols, drop = FALSE])
existing_label_cols <- intersect(setdiff(label_cols, "module_id"), colnames(module_df))
if (length(existing_label_cols) > 0) {
  module_df <- module_df[, setdiff(colnames(module_df), existing_label_cols), drop = FALSE]
}
module_df <- merge(module_df, label_map, by = "module_id", all.x = TRUE, sort = FALSE)

tf_all <- marker_df[is_tf(marker_df$TAP_family), , drop = FALSE]
tf_all <- merge(tf_all, module_df, by = "gene_id", all.x = TRUE, sort = FALSE)

# (1) Marker genes that are TFs (strict markers by marker_flag).
marker_tf <- tf_all[tf_all$marker_flag == TRUE, , drop = FALSE]
marker_tf <- marker_tf[order(marker_tf$marker_group, marker_tf$TAP_family, -abs(marker_tf$best_log2fc)), , drop = FALSE]
write_tsv(marker_tf, file.path(out_dir, "marker_tfs.tsv"))

if (nrow(marker_tf) > 0) {
  marker_tf_group <- aggregate(
    gene_id ~ marker_group + TAP_family,
    data = marker_tf,
    FUN = length
  )
  colnames(marker_tf_group)[colnames(marker_tf_group) == "gene_id"] <- "n_marker_tfs"
  marker_tf_group <- marker_tf_group[order(marker_tf_group$marker_group, -marker_tf_group$n_marker_tfs), , drop = FALSE]
} else {
  marker_tf_group <- data.frame(marker_group = character(), TAP_family = character(), n_marker_tfs = integer(), stringsAsFactors = FALSE)
}
write_tsv(marker_tf_group, file.path(out_dir, "marker_tfs_by_group_summary.tsv"))

# (2) WGCNA module correlation proxy for TFs: kME per TF.
tf_module_assoc <- tf_all[, c(
  "gene_id", "TAP_family", "Preferred_name", "Description",
  "module_id", "module_rank", "module_label", "module_label_generation", "module_label_context",
  "dominant_stage", "dominant_generation", "second_stage", "second_generation",
  "kME", "hub_flag",
  "marker_group", "marker_class", "marker_flag"
), drop = FALSE]
tf_module_assoc <- tf_module_assoc[order(tf_module_assoc$module_id, -abs(tf_module_assoc$kME)), , drop = FALSE]
write_tsv(tf_module_assoc, file.path(out_dir, "tf_module_associations.tsv"))

# Top TF candidates per module by |kME|.
if (nrow(tf_module_assoc) > 0) {
  split_mod <- split(tf_module_assoc, tf_module_assoc$module_id)
  top_mod <- do.call(rbind, lapply(split_mod, function(d) head(d[order(-abs(d$kME)), , drop = FALSE], 15)))
  rownames(top_mod) <- NULL
} else {
  top_mod <- tf_module_assoc
}
write_tsv(top_mod, file.path(out_dir, "module_top_tfs.tsv"))

# Optional: module-level TF overrepresentation test.
universe <- unique(marker_df$gene_id)
tf_universe <- unique(tf_all$gene_id)
mods <- split(module_df$gene_id, module_df$module_id)
rows <- list()
idx <- 0L
for (mid in names(mods)) {
  mg <- intersect(unique(mods[[mid]]), universe)
  if (length(mg) == 0) next
  a <- length(intersect(mg, tf_universe))
  b <- length(setdiff(mg, tf_universe))
  c <- length(setdiff(tf_universe, mg))
  d <- length(setdiff(universe, union(mg, tf_universe)))
  tab <- matrix(c(a, b, c, d), nrow = 2)
  ft <- stats::fisher.test(tab, alternative = "greater")
  idx <- idx + 1L
  rows[[idx]] <- data.frame(
    module_id = mid,
    module_rank = if ("module_rank" %in% colnames(module_df)) module_df$module_rank[match(mid, module_df$module_id)][1] else NA_character_,
    module_label = if ("module_label" %in% colnames(module_df)) module_df$module_label[match(mid, module_df$module_id)][1] else NA_character_,
    module_label_generation = if ("module_label_generation" %in% colnames(module_df)) module_df$module_label_generation[match(mid, module_df$module_id)][1] else NA_character_,
    module_label_context = if ("module_label_context" %in% colnames(module_df)) module_df$module_label_context[match(mid, module_df$module_id)][1] else NA_character_,
    dominant_stage = if ("dominant_stage" %in% colnames(module_df)) module_df$dominant_stage[match(mid, module_df$module_id)][1] else NA_character_,
    dominant_generation = if ("dominant_generation" %in% colnames(module_df)) module_df$dominant_generation[match(mid, module_df$module_id)][1] else NA_character_,
    second_stage = if ("second_stage" %in% colnames(module_df)) module_df$second_stage[match(mid, module_df$module_id)][1] else NA_character_,
    second_generation = if ("second_generation" %in% colnames(module_df)) module_df$second_generation[match(mid, module_df$module_id)][1] else NA_character_,
    module_size = length(mg),
    tf_in_module = a,
    tf_fraction = if (length(mg) > 0) a / length(mg) else NA_real_,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}
if (length(rows) > 0) {
  tf_enr <- do.call(rbind, rows)
  tf_enr$fdr <- stats::p.adjust(tf_enr$p_value, method = "BH")
  tf_enr <- tf_enr[order(tf_enr$fdr, tf_enr$p_value), , drop = FALSE]
} else {
  tf_enr <- data.frame(
    module_id = character(), module_rank = character(), module_label = character(),
    module_label_generation = character(), module_label_context = character(),
    dominant_stage = character(), dominant_generation = character(),
    second_stage = character(), second_generation = character(),
    module_size = integer(), tf_in_module = integer(), tf_fraction = numeric(),
    odds_ratio = numeric(), p_value = numeric(), fdr = numeric(),
    stringsAsFactors = FALSE
  )
}

# Add module stage/tissue interpretation column from module-marker enrichment results.
tf_enr$enriched_marker_groups <- NA_character_
tf_enr$enriched_marker_groups_fdr05 <- NA_character_
tf_enr$top_marker_groups <- NA_character_

if (nrow(cluster_marker_df) > 0 && all(c("cluster_id", "marker_group", "fdr") %in% colnames(cluster_marker_df))) {
  cluster_marker_df <- cluster_marker_df[is.finite(cluster_marker_df$fdr), , drop = FALSE]
  by_mod <- split(cluster_marker_df, cluster_marker_df$cluster_id)
  for (i in seq_len(nrow(tf_enr))) {
    mid <- tf_enr$module_id[i]
    if (!mid %in% names(by_mod)) next
    d <- by_mod[[mid]]
    d <- d[order(d$fdr), , drop = FALSE]
    if (nrow(d) == 0) next

    sig <- unique(d$marker_group[d$fdr <= 0.05])
    if (length(sig) > 0) {
      tf_enr$enriched_marker_groups_fdr05[i] <- paste(sig, collapse = ",")
    }

    all_ord <- unique(d$marker_group)
    tf_enr$enriched_marker_groups[i] <- paste(all_ord, collapse = ",")
    tf_enr$top_marker_groups[i] <- paste(head(all_ord, 3), collapse = ",")
  }
}

write_tsv(tf_enr, file.path(out_dir, "module_tf_enrichment.tsv"))

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "marker_tfs.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "marker_tfs_by_group_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "tf_module_associations.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "module_top_tfs.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "module_tf_enrichment.tsv"), "\n", sep = "")
