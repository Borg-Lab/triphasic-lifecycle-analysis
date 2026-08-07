#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source("R/io.R")
source("R/qc.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))

out_dir <- file.path(paths$output_root, "generation")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

summary_file <- file.path(out_dir, "triphasy_up_down_chrom_cluster_summary.tsv")
blocks_file <- file.path(out_dir, "triphasy_up_down_chrom_proximity_blocks.tsv")
if (!file.exists(summary_file) || !file.exists(blocks_file)) {
  stop("Missing triphasy chromosomal clustering outputs. Run 74_triphasy_biased_molecular_features.R first.", call. = FALSE)
}

set_levels <- c(
  "GAM_core_biased", "GAM_core_downregulated",
  "CSP_biased", "CSP_downregulated",
  "TETRA_biased", "TETRA_downregulated"
)
set_labels <- c(
  GAM_core_biased = "GAM biased",
  GAM_core_downregulated = "GAM repressed",
  CSP_biased = "Carposporophyte biased",
  CSP_downregulated = "Carposporophyte repressed",
  TETRA_biased = "Tetrasporophyte biased",
  TETRA_downregulated = "Tetrasporophyte repressed"
)
set_pal <- c(
  GAM_core_biased = "#8E2F3F",
  GAM_core_downregulated = "#CFA4A8",
  CSP_biased = "#6A3D9A",
  CSP_downregulated = "#B8A6D9",
  TETRA_biased = "#8FA33E",
  TETRA_downregulated = "#C8D48C"
)

chrom_order <- function(x) {
  autosomes <- grep("^Chr[0-9]+$", unique(x), value = TRUE)
  autosomes <- autosomes[order(as.integer(sub("^Chr", "", autosomes)))]
  c(autosomes, intersect(c("ChrUf", "ChrVm"), unique(x)))
}

coords <- filter_coordinates(read_coordinates(paths), analysis)
all_chr_levels <- chrom_order(coords$chr)

summary_df <- read.delim(summary_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(set_id %in% set_levels) %>%
  mutate(
    set_id = factor(set_id, levels = set_levels),
    set_label = factor(unname(set_labels[as.character(set_id)]), levels = rev(unname(set_labels[set_levels]))),
    generation = case_when(
      grepl("^GAM_core", as.character(set_id)) ~ "GAM",
      grepl("^CSP", as.character(set_id)) ~ "Carposporophyte",
      grepl("^TETRA", as.character(set_id)) ~ "Tetrasporophyte",
      TRUE ~ as.character(set_id)
    ),
    direction = if_else(grepl("_downregulated$", as.character(set_id)), "repressed", "biased"),
    generation = factor(generation, levels = c("GAM", "Carposporophyte", "Tetrasporophyte")),
    direction = factor(direction, levels = c("repressed", "biased")),
    log2_observed_expected_nn = log2(obs_median_nn_dist / perm_median),
    fdr_label = case_when(
      fdr <= 0.001 ~ "FDR <= 0.001",
      fdr <= 0.01 ~ "FDR <= 0.01",
      fdr <= 0.05 ~ "FDR <= 0.05",
      TRUE ~ "ns"
    ),
    point_label = sprintf("%.1f%% blocks\n%s", pct_genes_in_proximity_blocks, fdr_label),
    text_vjust = if_else(log2(obs_median_nn_dist / perm_median) < 0, 1.25, -0.25)
  )

write_tsv(
  summary_df %>%
    select(
      set_id, set_label, n_genes, obs_median_nn_dist, perm_median,
      log2_observed_expected_nn, empirical_p, fdr, cluster_flag,
      n_proximity_blocks, n_genes_in_proximity_blocks,
      pct_genes_in_proximity_blocks, avg_block_span_bp, median_block_span_bp
    ),
  file.path(out_dir, "triphasy_up_down_chrom_clustering_plot_summary.tsv")
)

blocks <- read.delim(blocks_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(set_id %in% set_levels) %>%
  mutate(
    set_id = factor(set_id, levels = set_levels),
    set_label = factor(unname(set_labels[as.character(set_id)]), levels = unname(set_labels[set_levels])),
    block_size_class = case_when(
      n_genes <= 2 ~ "2-gene blocks",
      n_genes == 3 ~ "3-gene blocks",
      TRUE ~ ">=4-gene blocks"
    ),
    chr = factor(chr, levels = all_chr_levels)
  )

block_comp <- blocks %>%
  count(set_id, set_label, block_size_class, wt = n_genes, name = "n_genes") %>%
  bind_rows(
    summary_df %>%
      transmute(
        set_id, set_label,
        block_size_class = "not in blocks",
        n_genes = n_genes - n_genes_in_proximity_blocks
      )
  ) %>%
  mutate(
    block_size_class = factor(
      block_size_class,
      levels = c("not in blocks", "2-gene blocks", "3-gene blocks", ">=4-gene blocks")
    )
  ) %>%
  group_by(set_id, set_label) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup()
write_tsv(block_comp, file.path(out_dir, "triphasy_up_down_chrom_block_composition.tsv"))

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
  left_join(
    block_gene_chr %>%
      count(set_id, set_label, chr, name = "clustered_genes_on_chr"),
    by = c("set_id", "set_label", "chr")
  ) %>%
  left_join(
    summary_df %>% select(set_id, set_label, n_genes, n_genes_in_proximity_blocks),
    by = c("set_id", "set_label")
  ) %>%
  mutate(clustered_genes_on_chr = replace_na(clustered_genes_on_chr, 0L)) %>%
  mutate(
    pct_set_genes_clustered_on_chr = 100 * clustered_genes_on_chr / n_genes,
    pct_clustered_genes_on_chr = 100 * clustered_genes_on_chr / n_genes_in_proximity_blocks
  ) %>%
  arrange(set_id, chr)
write_tsv(chrom_contrib, file.path(out_dir, "triphasy_up_down_chromosome_contribution.tsv"))

p_strength_facet <- ggplot(summary_df, aes(x = log2_observed_expected_nn, y = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey35") +
  geom_segment(aes(x = 0, xend = log2_observed_expected_nn, yend = direction), linewidth = 0.65, color = "grey60") +
  geom_point(aes(size = pct_genes_in_proximity_blocks, fill = set_id), shape = 21, color = "black", stroke = 0.25, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.1f%%\n%s", pct_genes_in_proximity_blocks, fdr_label)),
    hjust = 1.08,
    size = 2.35,
    lineheight = 0.86
  ) +
  facet_wrap(~ generation, ncol = 1) +
  scale_fill_manual(values = set_pal) +
  scale_size_continuous(range = c(3.0, 6.2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.28, 0.08))) +
  labs(
    x = "log2(observed / expected nearest-neighbour distance)",
    y = NULL,
    title = "Triphasy Chromosomal Clustering by Generation",
    subtitle = "Values below zero indicate genes are closer than expected"
  ) +
  theme_bw(base_size = 10) +
  theme(
    aspect.ratio = 0.42,
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(fig_dir, "triphasy_up_down_chrom_clustering_strength_faceted_dotplot.pdf"), p_strength_facet, width = 5.8, height = 7.2)

chr_levels <- all_chr_levels
chr_levels <- all_chr_levels
chrom_tile <- chrom_contrib %>%
  mutate(
    chr = factor(as.character(chr), levels = chr_levels),
    set_label = factor(set_label, levels = rev(unname(set_labels[set_levels])))
  )

chr_fill_cap <- 10
p_chr <- ggplot(chrom_tile, aes(x = chr, y = set_label, fill = pct_clustered_genes_on_chr)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradientn(
    colors = c(
      "#F5F5F5", "#EFE4D4", "#E7D5BC", "#D3B089",
      "#BE946C", "#A56E55", "#8F3A3A", "#663333"
    ),
    values = rescale(c(0, 0.5, 1, 2, 3.5, 5, 7.5, chr_fill_cap)),
    limits = c(0, chr_fill_cap),
    oob = scales::squish,
    breaks = c(0, 1, 2.5, 5, 7.5, 10),
    labels = c("0", "1", "2.5", "5", "7.5", ">=10"),
    name = "% clustered\ngenes on chr"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Chromosome Contribution to Triphasy Proximity Blocks"
  ) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave(file.path(fig_dir, "triphasy_up_down_chromosome_contribution_heatmap.pdf"), p_chr, width = 10.5, height = 4.5)

# Schematic gene-neighbourhood view: largest proximity block for each triphasy set.
features_file <- file.path(out_dir, "triphasy_up_down_molecular_features.tsv")
if (file.exists(features_file)) {
  molecular_features <- read.delim(features_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(set %in% set_levels) %>%
    mutate(set_id = set)

  gene_master_file <- file.path(out_dir, "..", "tables", "gene_master_table.tsv")
  if (!file.exists(gene_master_file)) {
    gene_master_file <- file.path(paths$output_root, "tables", "gene_master_table.tsv")
  }
  gene_master <- read.delim(gene_master_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(gene_id = normalize_gene_id(gene_id))

  stage_sets <- list(
    GAM_core = paste0("mean_tpm_", c("Bm3235f_vegtissue", "Bm3235f_procarps", "Bm3235m_vegtissue", "Bm3235m_spermatangia")),
    CSP = paste0("mean_tpm_", c("Bm3235c_carposporophyte")),
    TETRA = paste0("mean_tpm_", c("Bm3235t_vegtissue", "Bm3235t_stichidia"))
  )

  score_block <- function(block_row) {
    genes <- strsplit(block_row$gene_ids, ";", fixed = TRUE)[[1]]
    d <- gene_master %>% filter(gene_id %in% genes)
    if (nrow(d) == 0) return(NULL)
    set_id_chr <- as.character(block_row$set_id)
    generation <- case_when(
      grepl("^GAM_core", set_id_chr) ~ "GAM_core",
      grepl("^CSP", set_id_chr) ~ "CSP",
      grepl("^TETRA", set_id_chr) ~ "TETRA",
      TRUE ~ NA_character_
    )
    target_cols <- intersect(stage_sets[[generation]], names(d))
    opposing_cols <- setdiff(intersect(unlist(stage_sets), names(d)), target_cols)
    target_mean <- mean(as.matrix(log2(d[, target_cols, drop = FALSE] + 1)), na.rm = TRUE)
    opposing_mean <- mean(as.matrix(log2(d[, opposing_cols, drop = FALSE] + 1)), na.rm = TRUE)
    direction <- if_else(grepl("_downregulated$", set_id_chr), "repressed", "biased")
    contrast <- if (direction == "biased") target_mean - opposing_mean else opposing_mean - target_mean
    data.frame(
      block_row,
      generation = generation,
      direction = direction,
      target_log2_tpm = target_mean,
      opposing_log2_tpm = opposing_mean,
      visual_contrast_log2 = contrast,
      stringsAsFactors = FALSE
    )
  }

  block_scores <- bind_rows(lapply(seq_len(nrow(blocks)), function(i) {
    score_block(blocks[i, ])
  })) %>%
    filter(is.finite(visual_contrast_log2))

  write_tsv(
    block_scores %>%
      mutate(set_label = unname(set_labels[as.character(set_id)])) %>%
      select(set_id, set_label, chr, block_id, block_start, block_end, span_bp, n_genes,
             generation, direction, target_log2_tpm, opposing_log2_tpm, visual_contrast_log2, gene_ids) %>%
      arrange(set_id, desc(visual_contrast_log2), desc(n_genes), span_bp),
    file.path(out_dir, "triphasy_up_down_proximity_block_expression_scores.tsv")
  )

  candidate_blocks <- block_scores %>%
    filter(n_genes >= 3) %>%
    arrange(desc(visual_contrast_log2), desc(n_genes), span_bp) %>%
    mutate(block_key = paste(chr, block_start, block_end, sep = ":"))

  # Curated browser-view examples: compact blocks with cleaner target/off-target
  # expression than the purely highest-contrast automatic choices.
  curated_browser_blocks <- data.frame(
    set_id = c(
      "GAM_core_biased", "GAM_core_downregulated",
      "CSP_biased", "CSP_downregulated"
    ),
    chr = c("Chr02", "Chr17", "Chr06", "Chr10"),
    block_start = c(8537473L, 3808799L, 27344315L, 3030031L),
    block_end = c(8544187L, 3813409L, 27345282L, 3068224L),
    stringsAsFactors = FALSE
  ) %>%
    mutate(block_key = paste(chr, block_start, block_end, sep = ":"))

  used_block_keys <- character()
  selected_rows <- list()
  for (sid in set_levels) {
    curated_key <- curated_browser_blocks$block_key[curated_browser_blocks$set_id == sid]
    d <- candidate_blocks %>%
      filter(as.character(set_id) == sid, block_key %in% curated_key)
    if (nrow(d) == 0) {
      d <- candidate_blocks %>%
        filter(as.character(set_id) == sid, !block_key %in% used_block_keys) %>%
        arrange(desc(visual_contrast_log2), desc(n_genes), span_bp)
    }
    if (nrow(d) == 0) {
      d <- candidate_blocks %>%
        filter(as.character(set_id) == sid) %>%
        arrange(desc(visual_contrast_log2), desc(n_genes), span_bp)
    }
    if (nrow(d) > 0) {
      selected_rows[[sid]] <- d[1, , drop = FALSE]
      used_block_keys <- c(used_block_keys, d$block_key[1])
    }
  }

  extreme_blocks <- bind_rows(selected_rows) %>%
    select(-block_key) %>%
    mutate(
      set_label = unname(set_labels[set_id]),
      panel_label = sprintf("%s | %s:%s-%s | %s genes, %.1f kb | contrast %.2f",
                            set_label, chr, block_start, block_end, n_genes, span_bp / 1000, visual_contrast_log2),
      panel_label = factor(panel_label, levels = rev(panel_label))
    )

  extreme_gene_rows <- extreme_blocks %>%
    select(set_id, set_label, panel_label, chr, block_id, block_start, block_end, span_bp, n_genes, gene_ids) %>%
    separate_rows(gene_ids, sep = ";") %>%
    rename(gene_id = gene_ids) %>%
    left_join(
      molecular_features %>%
        select(gene_id, set_id, start, end, strand, gene_length_bp, exon_count, cds_length_bp),
      by = c("gene_id", "set_id")
    ) %>%
    mutate(
      plot_start = if_else(strand == "-", end, start),
      plot_end = if_else(strand == "-", start, end),
      x = (plot_start - block_start) / 1000,
      xend = (plot_end - block_start) / 1000,
      midpoint = ((start + end) / 2 - block_start) / 1000,
      label_y = if_else(row_number() %% 2 == 0, 0.28, -0.28),
      set_id = factor(set_id, levels = set_levels)
    )

  write_tsv(
    extreme_blocks %>%
      select(set_id, set_label, chr, block_id, block_start, block_end, span_bp, n_genes,
             generation, direction, target_log2_tpm, opposing_log2_tpm, visual_contrast_log2, gene_ids),
    file.path(out_dir, "triphasy_up_down_extreme_cluster_schematic_blocks.tsv")
  )
  write_tsv(
    extreme_gene_rows,
    file.path(out_dir, "triphasy_up_down_extreme_cluster_schematic_genes.tsv")
  )

  bed_dir <- file.path(out_dir, "igv_tracks")
  ensure_dir(bed_dir)
  write.table(
    extreme_blocks %>%
      transmute(
        chr,
        start0 = pmax(0, block_start - 1),
        block_end,
        name = paste0("triphasy_", set_id, "_", n_genes, "genes")
      ),
    file = file.path(bed_dir, "triphasy_up_down_extreme_cluster_schematic_blocks.bed"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  p_cluster <- ggplot(extreme_gene_rows) +
    geom_segment(
      aes(x = 0, xend = span_bp / 1000, y = 0, yend = 0),
      color = "grey75", linewidth = 0.35
    ) +
    geom_segment(
      aes(x = x, xend = xend, y = 0, yend = 0, color = set_id),
      arrow = grid::arrow(type = "closed", length = grid::unit(0.08, "inches")),
      linewidth = 2.8,
      lineend = "round"
    ) +
    geom_text(
      aes(x = midpoint, y = label_y, label = gene_id),
      size = 2.2,
      angle = 35,
      hjust = 0.5,
      vjust = 0.5
    ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_x") +
    scale_color_manual(values = set_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.7, 0.7), breaks = NULL) +
    labs(
      x = "Position within selected proximity block (kb)",
      y = NULL,
      title = "Extreme Triphasy Proximity Blocks",
      subtitle = "Largest proximity block selected for each biased/repressed generation gene set"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      axis.text.y = element_blank()
    )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_extreme_cluster_schematic.pdf"), p_cluster, width = 11, height = 8.4)
}

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "triphasy_up_down_chrom_clustering_plot_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_chrom_block_composition.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_chromosome_contribution.tsv"), "\n", sep = "")
if (file.exists(file.path(out_dir, "triphasy_up_down_extreme_cluster_schematic_blocks.tsv"))) {
  cat(" - ", file.path(out_dir, "triphasy_up_down_extreme_cluster_schematic_blocks.tsv"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "triphasy_up_down_extreme_cluster_schematic_genes.tsv"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "igv_tracks", "triphasy_up_down_extreme_cluster_schematic_blocks.bed"), "\n", sep = "")
  cat(" - ", file.path(fig_dir, "triphasy_up_down_extreme_cluster_schematic.pdf"), "\n", sep = "")
}
cat(" - ", file.path(fig_dir, "triphasy_up_down_chrom_clustering_strength_faceted_dotplot.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_chromosome_contribution_heatmap.pdf"), "\n", sep = "")
