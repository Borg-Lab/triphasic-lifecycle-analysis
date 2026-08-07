#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_dir <- file.path(paths$output_root, "generation")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

biased_file <- file.path(out_dir, "triphasy_biased_gene_set_membership.tsv")
down_file <- file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv")
if (!file.exists(biased_file)) stop("Missing ", biased_file, call. = FALSE)
if (!file.exists(down_file)) stop("Missing ", down_file, call. = FALSE)

set_levels <- c(
  "GAM_core_biased", "GAM_core_downregulated",
  "CSP_biased", "CSP_downregulated",
  "TETRA_biased", "TETRA_downregulated"
)
set_labels <- c(
  GAM_core_biased = "GAM biased",
  CSP_biased = "Carposporophyte biased",
  TETRA_biased = "Tetrasporophyte biased",
  GAM_core_downregulated = "GAM repressed",
  CSP_downregulated = "Carposporophyte repressed",
  TETRA_downregulated = "Tetrasporophyte repressed",
  genome_background = "All protein-coding genes"
)
set_cols <- c(
  GAM_core_biased = "#8E2F3F",
  CSP_biased = "#6A3D9A",
  TETRA_biased = "#8FA33E",
  GAM_core_downregulated = "#C99AA3",
  CSP_downregulated = "#B79AC8",
  TETRA_downregulated = "#C8D38A",
  genome_background = "#8A8A8A"
)

membership <- bind_rows(
  read.delim(biased_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
  read.delim(down_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
) %>%
  mutate(gene_id = normalize_gene_id(gene_id)) %>%
  filter(set %in% set_levels) %>%
  distinct(gene_id, set)

anno_raw <- read_functional_annotation(paths)
if (!all(c("gene_id", "TAP_family") %in% colnames(anno_raw))) {
  stop("Functional annotation must contain gene_id and TAP_family columns.", call. = FALSE)
}

# Use the protein-coding gene universe for TAP enrichment. Genes without a TAP
# annotation are treated as non-TAP, rather than being dropped from the test.
universe <- get_protein_coding_gene_ids(paths)
anno_calls <- anno_raw %>%
  transmute(
    gene_id = normalize_gene_id(gene_id),
    TAP_family = as.character(TAP_family)
  ) %>%
  filter(!is.na(gene_id), nzchar(gene_id), gene_id %in% universe) %>%
  distinct(gene_id, .keep_all = TRUE)
anno <- data.frame(gene_id = universe, stringsAsFactors = FALSE) %>%
  left_join(anno_calls, by = "gene_id") %>%
  mutate(
    TAP_family = ifelse(is.na(TAP_family) | !nzchar(trimws(TAP_family)) | trimws(TAP_family) %in% c("-", "NA"), "non_TAP", TAP_family),
    tap_status = ifelse(TAP_family == "non_TAP", "non-TAP", "TAP")
  )

membership <- membership %>% filter(gene_id %in% universe)

summary_rows <- bind_rows(
  lapply(set_levels, function(s) {
    gids <- unique(membership$gene_id[membership$set == s])
    tab <- anno %>%
      filter(gene_id %in% gids) %>%
      count(tap_status, name = "n_genes")
    data.frame(set = s, tab, stringsAsFactors = FALSE)
  }),
  anno %>% count(tap_status, name = "n_genes") %>% mutate(set = "genome_background") %>% select(set, tap_status, n_genes)
) %>%
  tidyr::complete(set = c(set_levels, "genome_background"), tap_status = c("TAP", "non-TAP"), fill = list(n_genes = 0L)) %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup() %>%
  mutate(
    set = factor(set, levels = c(set_levels, "genome_background")),
    set_label = factor(unname(set_labels[as.character(set)]), levels = unname(set_labels[c(set_levels, "genome_background")]))
  )
write_tsv(summary_rows, file.path(out_dir, "triphasy_up_down_tap_proportion_summary.tsv"))

stat_rows <- lapply(set_levels, function(s) {
  gids <- unique(membership$gene_id[membership$set == s])
  in_set <- universe %in% gids
  is_tap <- anno$tap_status[match(universe, anno$gene_id)] == "TAP"
  mat <- matrix(
    c(
      sum(in_set & is_tap),
      sum(in_set & !is_tap),
      sum(!in_set & is_tap),
      sum(!in_set & !is_tap)
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("set", "background_excluding_set"), c("TAP", "non_TAP"))
  )
  ft <- fisher.test(mat, alternative = "two.sided")
  data.frame(
    set = s,
    set_label = unname(set_labels[s]),
    n_set_genes = sum(in_set),
    n_set_taps = sum(in_set & is_tap),
    prop_set_taps = sum(in_set & is_tap) / max(1, sum(in_set)),
    n_background_genes_excluding_set = sum(!in_set),
    n_background_taps_excluding_set = sum(!in_set & is_tap),
    prop_background_taps_excluding_set = sum(!in_set & is_tap) / max(1, sum(!in_set)),
    odds_ratio = unname(ft$estimate),
    log2_odds_ratio = log2(unname(ft$estimate)),
    direction = ifelse(unname(ft$estimate) >= 1, "enriched", "depleted"),
    fisher_p = ft$p.value,
    stringsAsFactors = FALSE
  )
})
stats <- bind_rows(stat_rows) %>%
  mutate(
    fdr = p.adjust(fisher_p, method = "BH"),
    direction = ifelse(log2_odds_ratio >= 0, "enriched", "depleted")
  ) %>%
  arrange(fdr, fisher_p)
write_tsv(stats, file.path(out_dir, "triphasy_up_down_tap_enrichment_fisher.tsv"))

family_summary <- membership %>%
  left_join(anno[, c("gene_id", "TAP_family", "tap_status")], by = "gene_id") %>%
  mutate(TAP_family = ifelse(is.na(TAP_family), "non_TAP", TAP_family)) %>%
  filter(tap_status == "TAP") %>%
  count(set, TAP_family, name = "n_genes") %>%
  arrange(set, desc(n_genes), TAP_family)
write_tsv(family_summary, file.path(out_dir, "triphasy_up_down_tap_family_summary.tsv"))

plot_df <- summary_rows %>%
  filter(tap_status == "TAP") %>%
  left_join(stats[, c("set", "fdr", "log2_odds_ratio", "direction")], by = "set") %>%
  mutate(
    fill_col = unname(set_cols[as.character(set)]),
    fdr_label = case_when(
      is.na(fdr) ~ "",
      fdr <= 0.001 ~ "***",
      fdr <= 0.01 ~ "**",
      fdr <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = sprintf("%.1f%%\n%s/%s", 100 * proportion, n_genes, total_genes),
    stat_label = ifelse(is.na(fdr), "", sprintf("%s", fdr_label))
  )

p <- ggplot(plot_df, aes(x = set_label, y = proportion, fill = set)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 2.85, lineheight = 0.9, color = "grey15") +
  geom_text(aes(y = proportion + max(plot_df$proportion, na.rm = TRUE) * 0.085, label = stat_label), size = 3.2, color = "grey10", na.rm = TRUE) +
  scale_fill_manual(values = set_cols, drop = FALSE) +
  scale_y_continuous(labels = function(x) sprintf("%d%%", round(100 * x)), expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(ylim = c(0, max(plot_df$proportion, na.rm = TRUE) * 1.28)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "TAP genes (%)", caption = "Labels show TAP percentage and n TAP / total genes; * FDR <= 0.05, ** FDR <= 0.01, *** FDR <= 0.001 (two-sided Fisher test).")

# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_tap_proportion_barplot.pdf"), p, width = 7.4, height = 4.4, useDingbats = FALSE)

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "triphasy_up_down_tap_proportion_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_tap_enrichment_fisher.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_tap_family_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_tap_proportion_barplot.pdf"), "\n", sep = "")
