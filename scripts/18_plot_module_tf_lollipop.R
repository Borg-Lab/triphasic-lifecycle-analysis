#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
})

source("R/io.R")

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path), call. = FALSE)
  read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_root <- paths$output_root
out_dir <- file.path(out_root, "modules")
ensure_dir(out_dir)

marker_tf <- read_tsv(file.path(out_root, "tables", "marker_tfs.tsv"))
tf_assoc <- read_tsv(file.path(out_root, "tables", "tf_module_associations.tsv"))
module_top <- read_tsv(file.path(out_root, "tables", "module_top_tfs.tsv"))
module_tf_enr <- read_tsv(file.path(out_root, "tables", "module_tf_enrichment.tsv"))
module_label_map_path <- file.path(out_root, "modules", "module_label_map.tsv")
module_label_map <- if (file.exists(module_label_map_path)) read_tsv(module_label_map_path) else data.frame()
cluster_enr <- read_tsv(file.path(out_root, "enrichment", "cluster_marker_enrichment.tsv"))

if (nrow(tf_assoc) == 0) stop("No rows in tf_module_associations.tsv", call. = FALSE)

module_display <- function(df, id_col = "module_id", label_col = "module_label") {
  out <- as.character(df[[id_col]])
  if (label_col %in% colnames(df)) {
    lab <- as.character(df[[label_col]])
    keep <- !is.na(lab) & nzchar(lab)
    out[keep] <- lab[keep]
  }
  out
}

# Non-manuscript TF network/summary plots are omitted; the manuscript panel is the lollipop below.

# ---------------------------
# 2) Module lollipop of top TF candidates
# ---------------------------
if (nrow(module_top) > 0) {
  d <- module_top
  d$module_display <- module_display(d)
  d$module_facet <- if ("module_rank" %in% colnames(d)) as.character(d$module_rank) else d$module_display
  if (nrow(module_label_map) > 0 && "module_rank" %in% colnames(module_label_map)) {
    module_levels <- as.character(module_label_map$module_rank[order(as.integer(sub("^M", "", module_label_map$module_rank)))])
  } else {
    module_levels <- unique(d$module_facet[order(as.integer(sub("^M", "", d$module_facet)))])
  }
  d$module_facet <- factor(d$module_facet, levels = module_levels)
  d$abs_kME <- abs(d$kME)
  d <- d[order(d$module_facet, -d$abs_kME), , drop = FALSE]
  d$label <- paste0(d$gene_id, " (", d$TAP_family, ")")
  d$label <- factor(d$label, levels = rev(unique(d$label)))

  p_lol <- ggplot2::ggplot(d, ggplot2::aes(x = abs_kME, y = label, color = TAP_family)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = abs_kME, yend = label), color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~module_facet, scales = "free_y", ncol = 4, drop = FALSE) +
    ggplot2::labs(
      title = "Top TF candidates per module",
      subtitle = "Ranked by |kME|",
      x = "|kME|", y = "TF gene (family)"
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing = grid::unit(0.8, "lines"),
      aspect.ratio = 1
    )
  ggplot2::ggsave(file.path(out_dir, "module_top_tfs_lollipop.pdf"), p_lol, width = 12, height = 16)
}

cat("Wrote manuscript TF lollipop plot to: ", file.path(out_dir, "module_top_tfs_lollipop.pdf"), "\n", sep = "")
cat("Wrote TF association plots to: ", out_dir, "\n", sep = "")
