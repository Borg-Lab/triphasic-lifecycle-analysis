#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source("R/io.R")
source("R/qc.R")
source("R/chrom_cluster.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
out_dir <- file.path(paths$output_root, "sex_bias")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

sets_file <- file.path(out_dir, "sex_biased_gene_sets.tsv")
if (!file.exists(sets_file)) stop("Missing ", sets_file, ". Run 79_sex_biased_expression_uv_analysis.R first.", call. = FALSE)

set_levels <- c("female_biased", "male_biased")
set_labels <- c(female_biased = "Female-biased genes", male_biased = "Male-biased genes")
set_cols <- c(female_biased = "#8E2F3F", male_biased = "#B8860B")

sex_sets <- read.delim(sets_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(gene_id = normalize_gene_id(gene_id)) %>%
  filter(set %in% set_levels) %>%
  distinct(gene_id, set)

coords <- read_coordinates(paths) %>%
  filter_coordinates(analysis) %>%
  filter(!chr %in% c("ChrMt", "ChrPt"))
coords$gene_id <- normalize_gene_id(coords$gene_id)

extra_sets <- split(sex_sets$gene_id, sex_sets$set)
cluster_obj <- run_chromosomal_clustering(
  module_df = data.frame(gene_id = character(), module_id = character()),
  coords_df = coords,
  analysis_cfg = analysis,
  extra_sets = extra_sets
)

summary_df <- cluster_obj$summary %>%
  filter(set_id %in% set_levels) %>%
  mutate(
    set = set_id,
    set_label = unname(set_labels[set_id]),
    set_id = factor(set_id, levels = set_levels),
    set_label = factor(set_label, levels = rev(unname(set_labels[set_levels]))),
    clustering_strength_log2ratio = log2(perm_median / obs_median_nn_dist),
    log2_observed_expected_nn = log2(obs_median_nn_dist / perm_median),
    fdr_label = case_when(
      fdr <= 0.001 ~ "FDR <= 0.001",
      fdr <= 0.01 ~ "FDR <= 0.01",
      fdr <= 0.05 ~ "FDR <= 0.05",
      TRUE ~ "ns"
    ),
    point_label = sprintf("%.1f%% blocks\n%s", pct_genes_in_proximity_blocks, fdr_label)
  )

blocks <- cluster_obj$blocks %>%
  filter(set_id %in% set_levels) %>%
  mutate(
    set = set_id,
    set_label = unname(set_labels[set_id]),
    set_id = factor(set_id, levels = set_levels),
    set_label = factor(set_label, levels = unname(set_labels[set_levels])),
    block_size_class = case_when(
      n_genes <= 2 ~ "2-gene blocks",
      n_genes == 3 ~ "3-gene blocks",
      TRUE ~ ">=4-gene blocks"
    )
  )

write_tsv(summary_df, file.path(out_dir, "sex_biased_chrom_cluster_summary.tsv"))
write_tsv(blocks, file.path(out_dir, "sex_biased_chrom_proximity_blocks.tsv"))
write_tsv(cluster_obj$members %>% filter(set_id %in% set_levels), file.path(out_dir, "sex_biased_chrom_cluster_members.tsv"))

block_comp <- blocks %>%
  count(set_id, set_label, block_size_class, wt = n_genes, name = "n_genes") %>%
  bind_rows(
    summary_df %>%
      transmute(set_id, set_label, block_size_class = "not in blocks", n_genes = n_genes - n_genes_in_proximity_blocks)
  ) %>%
  mutate(block_size_class = factor(block_size_class, levels = c("not in blocks", "2-gene blocks", "3-gene blocks", ">=4-gene blocks"))) %>%
  group_by(set_id, set_label) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup()
write_tsv(block_comp, file.path(out_dir, "sex_biased_chrom_block_composition.tsv"))

chrom_order <- function(x) {
  u <- unique(as.character(x))
  autosomes <- grep("^Chr[0-9]+$", u, value = TRUE)
  autosomes <- autosomes[order(as.integer(sub("^Chr", "", autosomes)))]
  c(autosomes, intersect(c("ChrUf", "ChrVm"), u), setdiff(u, c(autosomes, "ChrUf", "ChrVm")))
}
all_chr_levels <- chrom_order(coords$chr)

block_gene_chr <- blocks %>%
  separate_rows(gene_ids, sep = ";") %>%
  distinct(set_id, set_label, chr, gene_id = gene_ids)
chrom_template <- expand.grid(
  set_id = factor(set_levels, levels = set_levels),
  chr = factor(all_chr_levels, levels = all_chr_levels),
  stringsAsFactors = FALSE
) %>%
  mutate(set_label = factor(unname(set_labels[as.character(set_id)]), levels = unname(set_labels[set_levels])))

chrom_contrib <- chrom_template %>%
  left_join(block_gene_chr %>% count(set_id, set_label, chr, name = "clustered_genes_on_chr"), by = c("set_id", "set_label", "chr")) %>%
  left_join(summary_df %>% select(set_id, set_label, n_genes, n_genes_in_proximity_blocks), by = c("set_id", "set_label")) %>%
  mutate(clustered_genes_on_chr = replace_na(clustered_genes_on_chr, 0L)) %>%
  mutate(
    pct_set_genes_clustered_on_chr = 100 * clustered_genes_on_chr / n_genes,
    pct_clustered_genes_on_chr = 100 * clustered_genes_on_chr / n_genes_in_proximity_blocks
  ) %>%
  arrange(set_id, chr)
write_tsv(chrom_contrib, file.path(out_dir, "sex_biased_chromosome_contribution.tsv"))

# Sex-biased chromosomal clustering figures are not manuscript panels; this script
# exports the tables used by the combined supplemental chromosomal clustering table.

cat("Wrote sex-biased chromosomal clustering outputs to:\n")
cat(" - ", out_dir, "\n", sep = "")
cat(" - ", fig_dir, "\n", sep = "")
