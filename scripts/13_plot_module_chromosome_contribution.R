#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_chr <- file.path(paths$output_root, "chrom_cluster")
out_fig <- file.path(paths$output_root, "figures")
ensure_dir(out_chr)
ensure_dir(out_fig)

members <- safe_read_table(file.path(out_chr, "chrom_cluster_members.tsv"))
blocks <- safe_read_table(file.path(out_chr, "chrom_proximity_blocks.tsv"))
labels <- safe_read_table(file.path(paths$output_root, "modules", "module_label_map.tsv"))
coords <- read_coordinates(paths)

if (nrow(members) == 0 || nrow(blocks) == 0 || nrow(coords) == 0) {
  stop("Missing chromosomal clustering inputs (members/blocks/coords).", call. = FALSE)
}

coords <- coords[, c("gene_id", "chr", "start", "end"), drop = FALSE]
coords$gene_id <- as.character(coords$gene_id)
coords$chr <- as.character(coords$chr)
coords <- coords[!duplicated(coords$gene_id), , drop = FALSE]

members <- merge(members, coords[, c("gene_id", "chr"), drop = FALSE], by = "gene_id", all.x = TRUE, sort = FALSE)
members <- members[!is.na(members$chr) & nzchar(members$chr), , drop = FALSE]

if (nrow(labels) > 0 && all(c("module_id", "module_label") %in% colnames(labels))) {
  keep_cols <- c("module_id", "module_label", "module_rank")
  keep_cols <- keep_cols[keep_cols %in% colnames(labels)]
  keep <- unique(labels[, keep_cols, drop = FALSE])
  members <- merge(members, keep, by.x = "set_id", by.y = "module_id", all.x = TRUE, sort = FALSE)
  blocks <- merge(blocks, keep, by.x = "set_id", by.y = "module_id", all.x = TRUE, sort = FALSE)
  members$module_id_raw <- members$set_id
  blocks$module_id_raw <- blocks$set_id
  members$module <- ifelse(!is.na(members$module_label) & nzchar(members$module_label), members$module_label, members$set_id)
  blocks$module <- ifelse(!is.na(blocks$module_label) & nzchar(blocks$module_label), blocks$module_label, blocks$set_id)
} else {
  members$module_id_raw <- members$set_id
  blocks$module_id_raw <- blocks$set_id
  members$module_rank <- NA_character_
  blocks$module_rank <- NA_character_
  members$module_label <- members$set_id
  blocks$module_label <- blocks$set_id
  members$module <- members$set_id
  blocks$module <- blocks$set_id
}

# Build module numeric order (M01..M24) for plotting.
module_rank_map <- setNames(character(0), character(0))
if (nrow(labels) > 0 && all(c("module_label", "module_rank") %in% colnames(labels))) {
  module_rank_map <- setNames(as.character(labels$module_rank), as.character(labels$module_label))
}

module_tot <- aggregate(gene_id ~ module, data = members, FUN = function(x) length(unique(x)))
colnames(module_tot)[2] <- "module_genes_clustered"

chr_mod <- aggregate(gene_id ~ module + chr, data = members, FUN = function(x) length(unique(x)))
colnames(chr_mod)[3] <- "clustered_genes_on_chr"

blk_n <- aggregate(block_id ~ module + chr, data = blocks, FUN = function(x) length(unique(x)))
colnames(blk_n)[3] <- "n_blocks"
blk_span <- aggregate(span_bp ~ module + chr, data = blocks, FUN = median)
colnames(blk_span)[3] <- "median_block_span_bp"

contrib <- merge(chr_mod, module_tot, by = "module", all.x = TRUE, sort = FALSE)
contrib <- merge(contrib, blk_n, by = c("module", "chr"), all.x = TRUE, sort = FALSE)
contrib <- merge(contrib, blk_span, by = c("module", "chr"), all.x = TRUE, sort = FALSE)
contrib$n_blocks[is.na(contrib$n_blocks)] <- 0L
contrib$median_block_span_bp[is.na(contrib$median_block_span_bp)] <- 0
contrib$pct_module_clustered_on_chr <- 100 * contrib$clustered_genes_on_chr / pmax(contrib$module_genes_clustered, 1)
contrib$module_rank <- unname(module_rank_map[as.character(contrib$module)])
contrib$module_label <- as.character(contrib$module)
contrib$module_rank_num <- suppressWarnings(as.integer(sub("^M", "", contrib$module_rank)))
contrib <- contrib[order(contrib$module_rank_num, contrib$module_label, -contrib$pct_module_clustered_on_chr, contrib$chr, na.last = TRUE), , drop = FALSE]
contrib$module_rank_num <- NULL
write_tsv(contrib[, c("module_rank", "module_label", "module", "chr", "clustered_genes_on_chr", "module_genes_clustered", "n_blocks", "median_block_span_bp", "pct_module_clustered_on_chr"), drop = FALSE], file.path(out_chr, "module_chromosome_contribution.tsv"))

# Heatmap-style overview for all modules.
top_mod <- unique(contrib$module)
heat <- contrib[contrib$module %in% top_mod, , drop = FALSE]

all_chr <- sort(unique(as.character(coords$chr)))
auto_chr <- all_chr[grepl("^Chr[0-9]+$", all_chr)]
if (length(auto_chr) > 0) {
  auto_num <- as.integer(sub("^Chr", "", auto_chr))
  auto_chr <- auto_chr[order(auto_num)]
}
tail_chr <- intersect(c("ChrUf", "ChrVm"), all_chr)
mid_chr <- setdiff(all_chr, c(auto_chr, tail_chr))
chr_order <- c(auto_chr, mid_chr, tail_chr)
heat$chr <- factor(heat$chr, levels = chr_order)

# Prefer numeric module-rank ordering on y-axis.
if (length(module_rank_map) > 0) {
  top_rank <- unname(module_rank_map[top_mod])
  top_rank_num <- suppressWarnings(as.integer(sub("^M", "", top_rank)))
  ord_idx <- order(top_rank_num, na.last = TRUE)
  top_mod <- top_mod[ord_idx]
}
heat$module <- factor(heat$module, levels = rev(top_mod))

# Capped continuous scale: most values are <10%, with rare outliers >40%.
# Compressing at 12% preserves the high-signal outliers while making 1-10%
# differences visible in the heatmap.
fill_cap <- 12
heat$fill_pct_capped <- pmin(heat$pct_module_clustered_on_chr, fill_cap)

p_heat <- ggplot(heat, aes(x = chr, y = module, fill = fill_pct_capped)) +
  geom_tile(color = "white", linewidth = 0.2, height = 0.8) +
  scale_fill_gradientn(
    colours = c("#f7efe5", "#ead8c0", "#d3ad82", "#aa7557", "#663333"),
    values = scales::rescale(c(0, 1, 3, 6, fill_cap)),
    limits = c(0, fill_cap),
    breaks = c(0, 1, 3, 6, 12),
    labels = c("0", "1", "3", "6", ">=12"),
    oob = scales::squish
  ) +
  labs(
    title = "Module chromosome contribution (clustered genes)",
    subtitle = "All modules",
    x = "Chromosome",
    y = "Module",
    fill = "% of module\nclustered genes"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(out_fig, "module_chromosome_contribution_heatmap.pdf"), p_heat, width = 10.8, height = 10)

# Exploratory chromosome block-position plots are omitted from the publication repo.

cat("Wrote:\n")
cat(" - ", file.path(out_chr, "module_chromosome_contribution.tsv"), "\n", sep = "")
cat(" - ", file.path(out_fig, "module_chromosome_contribution_heatmap.pdf"), "\n", sep = "")
