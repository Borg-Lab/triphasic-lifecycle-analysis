#!/usr/bin/env Rscript

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_root <- paths$output_root
chrom_dir <- file.path(out_root, "chrom_cluster")
ensure_dir(chrom_dir)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

write_table <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

value_or_empty <- function(x, name) {
  if (name %in% names(x)) as.character(x[[name]]) else rep("", nrow(x))
}

# -----------------------------------------------------------------------------
# 1) Module proximity blocks with renamed Mxx labels.
# -----------------------------------------------------------------------------
blocks_path <- file.path(chrom_dir, "chrom_proximity_blocks.tsv")
label_path <- file.path(out_root, "modules", "module_label_map.tsv")
blocks <- read_tsv(blocks_path)
labels <- read_tsv(label_path)

label_keep <- intersect(c("module_id", "module_rank", "module_label"), colnames(labels))
labels <- labels[, label_keep, drop = FALSE]
blocks_labeled <- merge(blocks, labels, by.x = "set_id", by.y = "module_id", all.x = TRUE, sort = FALSE)
blocks_labeled$module_rank[is.na(blocks_labeled$module_rank) | !nzchar(blocks_labeled$module_rank)] <- blocks_labeled$set_id[is.na(blocks_labeled$module_rank) | !nzchar(blocks_labeled$module_rank)]
blocks_labeled$module_label[is.na(blocks_labeled$module_label) | !nzchar(blocks_labeled$module_label)] <- blocks_labeled$set_id[is.na(blocks_labeled$module_label) | !nzchar(blocks_labeled$module_label)]
blocks_labeled$cluster_number <- blocks_labeled$block_id
blocks_labeled$cluster_id <- ifelse(
  suppressWarnings(!is.na(as.integer(blocks_labeled$block_id))),
  paste0(blocks_labeled$module_rank, "_cluster_", sprintf("%03d", as.integer(blocks_labeled$block_id))),
  paste0(blocks_labeled$module_rank, "_cluster_", blocks_labeled$block_id)
)
module_blocks_out <- data.frame(
  set_type = "module",
  set_id_raw = blocks_labeled$set_id,
  module_rank = blocks_labeled$module_rank,
  set_label = blocks_labeled$module_label,
  cluster_number = blocks_labeled$cluster_number,
  cluster_id = blocks_labeled$cluster_id,
  chr = blocks_labeled$chr,
  cluster_start = blocks_labeled$block_start,
  cluster_end = blocks_labeled$block_end,
  span_bp = blocks_labeled$span_bp,
  n_genes = blocks_labeled$n_genes,
  gene_ids = blocks_labeled$gene_ids,
  stringsAsFactors = FALSE
)
module_blocks_path <- file.path(chrom_dir, "chrom_proximity_blocks_labeled.tsv")
write_table(module_blocks_out, module_blocks_path)

# -----------------------------------------------------------------------------
# 2) Combined proximity-block lookup across modules, triphasy sets, and SBGs.
# -----------------------------------------------------------------------------
standardize_blocks <- function(x, set_type) {
  if (set_type == "module") return(x)
  set_id <- value_or_empty(x, "set_id")
  set_label <- value_or_empty(x, "set_label")
  set_alt <- value_or_empty(x, "set")
  set_label[!nzchar(set_label)] <- set_alt[!nzchar(set_label)]
  set_label[!nzchar(set_label)] <- set_id[!nzchar(set_label)]
  block_id <- value_or_empty(x, "block_id")
  block_int <- suppressWarnings(as.integer(block_id))
  data.frame(
    set_type = set_type,
    set_id_raw = set_id,
    module_rank = "",
    set_label = set_label,
    cluster_number = block_id,
    cluster_id = ifelse(
      !is.na(block_int),
      paste0(set_id, "_cluster_", sprintf("%03d", block_int)),
      paste0(set_id, "_cluster_", block_id)
    ),
    chr = x$chr,
    cluster_start = x$block_start,
    cluster_end = x$block_end,
    span_bp = x$span_bp,
    n_genes = x$n_genes,
    gene_ids = x$gene_ids,
    stringsAsFactors = FALSE
  )
}

trip_blocks <- read_tsv(file.path(out_root, "generation", "triphasy_up_down_chrom_proximity_blocks.tsv"))
sex_blocks <- read_tsv(file.path(out_root, "sex_bias", "sex_biased_chrom_proximity_blocks.tsv"))
all_blocks <- rbind(
  module_blocks_out,
  standardize_blocks(trip_blocks, "triphasy"),
  standardize_blocks(sex_blocks, "sex_bias")
)
all_blocks_path <- file.path(chrom_dir, "all_chrom_proximity_blocks_labeled.tsv")
write_table(all_blocks, all_blocks_path)

# -----------------------------------------------------------------------------
# 3) Combined chromosomal clustering enrichment/statistics summary.
# -----------------------------------------------------------------------------
standardize_stats <- function(x, set_type) {
  n <- nrow(x)
  set_id_raw <- if (set_type == "module") value_or_empty(x, "module_id") else value_or_empty(x, "set_id")
  set_id_raw[!nzchar(set_id_raw)] <- value_or_empty(x, "set")[!nzchar(set_id_raw)]
  module_rank <- if (set_type == "module") value_or_empty(x, "module_rank") else rep("", n)
  set_label <- if (set_type == "module") value_or_empty(x, "module_label") else value_or_empty(x, "set_label")
  set_label[!nzchar(set_label)] <- value_or_empty(x, "set")[!nzchar(set_label)]
  set_label[!nzchar(set_label)] <- set_id_raw[!nzchar(set_label)]
  set_label_generation <- if (set_type == "module") value_or_empty(x, "set_label_generation") else rep("", n)
  if (set_type == "module" && !any(nzchar(set_label_generation))) set_label_generation <- value_or_empty(x, "module_label_generation")
  log2_oe <- value_or_empty(x, "log2_observed_expected_nn")
  if (!any(nzchar(log2_oe)) && "clustering_strength_log2ratio" %in% colnames(x)) {
    log2_oe <- as.character(-as.numeric(x$clustering_strength_log2ratio))
  }
  data.frame(
    set_type = set_type,
    set_id_raw = set_id_raw,
    module_rank = module_rank,
    set_label = set_label,
    set_label_generation = set_label_generation,
    n_genes = value_or_empty(x, "n_genes"),
    obs_median_nn_dist = value_or_empty(x, "obs_median_nn_dist"),
    perm_median = value_or_empty(x, "perm_median"),
    z_score = value_or_empty(x, "z_score"),
    empirical_p = value_or_empty(x, "empirical_p"),
    fdr = value_or_empty(x, "fdr"),
    cluster_flag = value_or_empty(x, "cluster_flag"),
    n_proximity_blocks = value_or_empty(x, "n_proximity_blocks"),
    n_genes_in_proximity_blocks = value_or_empty(x, "n_genes_in_proximity_blocks"),
    pct_genes_in_proximity_blocks = value_or_empty(x, "pct_genes_in_proximity_blocks"),
    avg_block_span_bp = value_or_empty(x, "avg_block_span_bp"),
    median_block_span_bp = value_or_empty(x, "median_block_span_bp"),
    clustering_strength_log2ratio = value_or_empty(x, "clustering_strength_log2ratio"),
    log2_observed_expected_nn = log2_oe,
    chr_scope = value_or_empty(x, "chr_scope"),
    stringsAsFactors = FALSE
  )
}

module_stats <- read_tsv(file.path(chrom_dir, "chrom_cluster_summary_labeled.tsv"))
trip_stats <- read_tsv(file.path(out_root, "generation", "triphasy_up_down_chrom_cluster_summary.tsv"))
sex_stats <- read_tsv(file.path(out_root, "sex_bias", "sex_biased_chrom_cluster_summary.tsv"))
combined_stats <- rbind(
  standardize_stats(module_stats, "module"),
  standardize_stats(trip_stats, "generation_identity"),
  standardize_stats(sex_stats, "sex_bias")
)
module_num <- suppressWarnings(as.integer(sub("^M", "", combined_stats$module_rank)))
combined_stats <- combined_stats[order(
  match(combined_stats$set_type, c("module", "generation_identity", "sex_bias")),
  ifelse(is.na(module_num), 999L, module_num),
  combined_stats$set_label
), , drop = FALSE]
stats_path <- file.path(chrom_dir, "combined_chrom_cluster_enrichment_stats.tsv")
write_table(combined_stats, stats_path)

cat("Wrote:\n")
cat(" - ", module_blocks_path, "\n", sep = "")
cat(" - ", all_blocks_path, "\n", sep = "")
cat(" - ", stats_path, "\n", sep = "")
