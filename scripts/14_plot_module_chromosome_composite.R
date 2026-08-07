#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))

in_order <- file.path(paths$output_root, "modules", "module_stage_enrichment_row_order.tsv")
in_map <- file.path(paths$output_root, "modules", "module_label_map.tsv")
in_counts <- file.path(paths$output_root, "modules", "module_gene_counts.tsv")
in_blocks <- file.path(paths$output_root, "chrom_cluster", "chrom_proximity_blocks.tsv")
in_cluster <- file.path(paths$output_root, "chrom_cluster", "chrom_cluster_summary_labeled.tsv")

out_fig <- file.path(paths$output_root, "figures", "module_chrom_composite_bars.pdf")
out_tab <- file.path(paths$output_root, "chrom_cluster", "module_chrom_composite_values.tsv")

for (f in c(in_order, in_map, in_counts, in_blocks, in_cluster)) {
  if (!file.exists(f)) stop(sprintf("Missing input: %s", f), call. = FALSE)
}

ensure_dir(dirname(out_fig))
ensure_dir(dirname(out_tab))

ord <- read.delim(in_order, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
map <- read.delim(in_map, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cnt <- read.delim(in_counts, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
blk <- read.delim(in_blocks, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cls <- read.delim(in_cluster, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("module_id") %in% colnames(ord))) stop("module_stage_enrichment_row_order.tsv needs module_id", call. = FALSE)
if (!all(c("module_id", "module_label") %in% colnames(map))) stop("module_label_map.tsv needs module_id,module_label", call. = FALSE)
if (!all(c("module_id", "n_genes") %in% colnames(cnt))) stop("module_gene_counts.tsv needs module_id,n_genes", call. = FALSE)
if (!all(c("set_id", "n_genes") %in% colnames(blk))) stop("chrom_proximity_blocks.tsv needs set_id,n_genes", call. = FALSE)
if (!all(c("set_id", "clustering_strength_log2ratio") %in% colnames(cls))) {
  stop("chrom_cluster_summary_labeled.tsv needs set_id,clustering_strength_log2ratio", call. = FALSE)
}

module_ids <- ord$module_id
module_ids <- module_ids[module_ids %in% map$module_id]
module_ids <- unique(module_ids)

lab_lookup <- setNames(map$module_label, map$module_id)
rank_lookup <- if ("module_rank" %in% colnames(map)) setNames(map$module_rank, map$module_id) else setNames(character(0), character(0))
module_ranks <- unname(rank_lookup[module_ids])
module_labels <- unname(lab_lookup[module_ids])
bad <- is.na(module_labels) | !nzchar(module_labels)
module_labels[bad] <- module_ranks[bad]
bad <- is.na(module_labels) | !nzchar(module_labels)
module_labels[bad] <- module_ids[bad]
module_factor_levels <- rev(module_labels)
display_lookup <- setNames(module_labels, module_ids)
rank_display_lookup <- setNames(module_ranks, module_ids)

## Panel 1: total genes in module
p1_df <- data.frame(
  module_id = module_ids,
  module_label = module_labels,
  n_genes_module = cnt$n_genes[match(module_ids, cnt$module_id)],
  stringsAsFactors = FALSE
)
p1_df$n_genes_module[is.na(p1_df$n_genes_module)] <- 0
p1_df$module_label <- factor(p1_df$module_label, levels = module_factor_levels)

## Panel 2: genes in proximity blocks by block-size class (2, 3, 4+)
b <- blk
b$module_id <- b$set_id
b <- b[b$module_id %in% module_ids, , drop = FALSE]
b$block_class <- ifelse(b$n_genes <= 2, "2", ifelse(b$n_genes == 3, "3", "4+"))

if (nrow(b) > 0) {
  agg <- aggregate(
    n_genes ~ module_id + block_class,
    data = b,
    FUN = sum,
    na.rm = TRUE
  )
} else {
  agg <- data.frame(module_id = character(0), block_class = character(0), n_genes = numeric(0), stringsAsFactors = FALSE)
}

# complete module x class grid
grid <- expand.grid(module_id = module_ids, block_class = c("2", "3", "4+"), stringsAsFactors = FALSE)
p2_df <- merge(grid, agg, by = c("module_id", "block_class"), all.x = TRUE, sort = FALSE)
p2_df$n_genes[is.na(p2_df$n_genes)] <- 0
p2_df$module_label <- unname(display_lookup[p2_df$module_id])
p2_df$module_label[is.na(p2_df$module_label) | !nzchar(p2_df$module_label)] <- p2_df$module_id[is.na(p2_df$module_label) | !nzchar(p2_df$module_label)]
p2_df$module_label <- factor(p2_df$module_label, levels = module_factor_levels)
p2_df$block_class <- factor(p2_df$block_class, levels = c("2", "3", "4+"))

## Panel 3: chromosomal clustering strength
p3_df <- data.frame(
  module_id = module_ids,
  module_label = module_labels,
  clustering_strength_log2ratio = cls$clustering_strength_log2ratio[match(module_ids, cls$set_id)],
  stringsAsFactors = FALSE
)
p3_df$clustering_strength_log2ratio[!is.finite(p3_df$clustering_strength_log2ratio)] <- 0
p3_df$significant <- FALSE
if ("cluster_flag" %in% colnames(cls)) {
  p3_df$significant <- as.logical(cls$cluster_flag[match(module_ids, cls$set_id)])
} else if ("fdr" %in% colnames(cls)) {
  p3_df$significant <- as.numeric(cls$fdr[match(module_ids, cls$set_id)]) <= 0.05
}
p3_df$significant[is.na(p3_df$significant)] <- FALSE
p3_df$sig_label <- ifelse(p3_df$significant, "significant", "not_significant")
p3_df$module_label <- factor(p3_df$module_label, levels = module_factor_levels)

block_totals <- aggregate(n_genes ~ module_id, data = p2_df, FUN = sum, na.rm = TRUE)
composite_out <- data.frame(
  module_id = module_ids,
  module_rank = module_ranks,
  module_label = module_labels,
  n_genes_module = p1_df$n_genes_module,
  clustering_strength_log2ratio = p3_df$clustering_strength_log2ratio,
  n_genes = block_totals$n_genes[match(module_ids, block_totals$module_id)],
  stringsAsFactors = FALSE
)
composite_out$n_genes[is.na(composite_out$n_genes)] <- 0
write_tsv(composite_out, out_tab)

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank()
  )

p1 <- ggplot(p1_df, aes(y = module_label, x = n_genes_module)) +
  geom_col(fill = "#4d4d4d", width = 0.8) +
  geom_text(aes(label = n_genes_module), hjust = -0.1, size = 2.7) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Genes in module", title = "Module size") +
  base_theme

p2 <- ggplot(p2_df, aes(y = module_label, x = n_genes, fill = block_class)) +
  geom_col(width = 0.8, position = "stack") +
  scale_fill_manual(values = c("2" = "#d9d9d9", "3" = "#9ecae1", "4+" = "#3182bd"), name = "Block size") +
  labs(x = "Genes in proximity blocks", title = "Block composition") +
  base_theme

# Normalize block composition to percentages to avoid domination by very large modules.
p2_tot <- aggregate(n_genes ~ module_label, data = p2_df, FUN = sum, na.rm = TRUE)
colnames(p2_tot)[2] <- "total_block_genes"
p2_df <- merge(p2_df, p2_tot, by = "module_label", all.x = TRUE, sort = FALSE)
p2_df$prop <- ifelse(p2_df$total_block_genes > 0, p2_df$n_genes / p2_df$total_block_genes, 0)
p2_df$module_label <- factor(p2_df$module_label, levels = module_factor_levels)

p2_lab <- unique(p2_df[, c("module_label", "total_block_genes"), drop = FALSE])

p3 <- ggplot(p3_df, aes(y = module_label, x = clustering_strength_log2ratio, fill = sig_label)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("significant" = "#a63603", "not_significant" = "#bdbdbd"), name = "Clustering") +
  labs(x = "Chromosomal clustering strength\nlog2(obs median NN / perm median NN)", title = "Clustering strength") +
  base_theme

if (requireNamespace("gridExtra", quietly = TRUE)) {
  p2 <- ggplot(p2_df, aes(y = module_label, x = prop, fill = block_class)) +
    geom_col(width = 0.8, position = "stack") +
    geom_text(
      data = p2_lab,
      aes(y = module_label, x = 1.02, label = total_block_genes),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.6
    ) +
    scale_fill_manual(values = c("2" = "#d9d9d9", "3" = "#9ecae1", "4+" = "#3182bd"), name = "Block size") +
    scale_x_continuous(
      limits = c(0, 1.18),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "25", "50", "75", "100")
    ) +
    labs(
      x = "Block-size composition (%)   [labels: total genes in blocks]",
      title = "Block composition"
    ) +
    base_theme

  grDevices::pdf(out_fig, width = 15, height = 10)
  gridExtra::grid.arrange(p1, p2, p3, ncol = 3, widths = c(1.0, 1.4, 1.1))
  grDevices::dev.off()
} else {
  # Fallback single panel if composition packages are unavailable.
  ggsave(out_fig, p2, width = 8, height = 10)
}

cat("Wrote:\n")
cat(" - ", out_fig, "\n", sep = "")
cat(" - ", out_tab, "\n", sep = "")
