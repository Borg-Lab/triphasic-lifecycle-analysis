suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
out_dir <- file.path(paths$output_root, "tai")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_tpm_path <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")
gene_age_path <- paths$gene_age_tsv

stopifnot(file.exists(group_tpm_path), file.exists(gene_age_path))

stage_levels <- c(
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

stage_tbl <- design %>%
  distinct(stage, generation, tissue, sex) %>%
  filter(stage %in% stage_levels) %>%
  mutate(
    stage_order = match(stage, stage_levels),
    stage_id = stage
  ) %>%
  arrange(stage_order)

missing_stages <- setdiff(stage_levels, stage_tbl$stage)
if (length(missing_stages) > 0) {
  stop("Missing stages in design: ", paste(missing_stages, collapse = ", "), call. = FALSE)
}

# Wide stage-mean TPM matrix
gtpm <- read_tsv(group_tpm_path, show_col_types = FALSE, progress = FALSE)

stage_cols <- paste0("mean_tpm_", stage_levels)
stage_cols <- intersect(stage_cols, names(gtpm))
if (length(stage_cols) == 0) stop("No mean_tpm_* columns found for current stage order.", call. = FALSE)

# Gene age table (no header): transcript_id, clade, rank, support
ages <- read_tsv(gene_age_path, col_names = FALSE, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    gene_id = sub("\\.t\\d+$", "", as.character(X1)),
    phylorank = suppressWarnings(as.integer(X3))
  ) %>%
  filter(!is.na(gene_id), gene_id != "", !is.na(phylorank), phylorank >= 1, phylorank <= 8) %>%
  distinct(gene_id, .keep_all = TRUE)

rank_labels <- c(
  "PS1: cellular organisms",
  "PS2: Eukaryota",
  "PS3: Rhodophyta",
  "PS4: Florideophyceae",
  "PS5: Rhodymeniophycidae",
  "PS6: Ceramiales",
  "PS7: Rhodomelaceae",
  "PS8: Bostrychia moritziana"
)

# Long expression table joined to phylorank
expr_long <- gtpm %>%
  select(gene_id, all_of(stage_cols)) %>%
  pivot_longer(cols = all_of(stage_cols), names_to = "stage", values_to = "mean_tpm") %>%
  mutate(
    stage = sub("^mean_tpm_", "", as.character(stage)),
    log2_tpm1 = log2(mean_tpm + 1)
  ) %>%
  inner_join(ages, by = "gene_id") %>%
  inner_join(stage_tbl, by = "stage")

# Summarise per phylorank x stage
sum_tbl <- expr_long %>%
  group_by(phylorank, stage_order, stage, stage_id, generation, tissue, sex) %>%
  summarise(
    n_genes = n(),
    mean_log2_tpm1 = mean(log2_tpm1, na.rm = TRUE),
    sd_log2_tpm1 = sd(log2_tpm1, na.rm = TRUE),
    se_log2_tpm1 = sd_log2_tpm1 / sqrt(pmax(n_genes, 1)),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2_tpm1 - se_log2_tpm1,
    ymax = mean_log2_tpm1 + se_log2_tpm1,
    phylorank = factor(phylorank, levels = 1:8, labels = rank_labels),
    stage_id = factor(stage_id, levels = stage_levels)
  )

write_tsv(sum_tbl, file.path(out_dir, "phylorank_lifecycle_expression_summary.tsv"))

stage_cols_plot <- c(
  "Bm3235f_vegtissue"       = "#AC5858",
  "Bm3235f_procarps"        = "#78211C",
  "Bm3235m_vegtissue"       = "#8EB8B2",
  "Bm3235m_spermatangia"    = "#5B9388",
  "Bm3235m_spermatia"       = "#387568",
  "Bm3235f_pericarp"        = "#B87B7E",
  "Bm3235f_cystocarps_imm"  = "#ACA4C7",
  "Bm3235f_cystocarps_mat"  = "#8172AA",
  "Bm3235c_carposporophyte" = "#5E3895",
  "Bm3235t_carpospores"     = "#442076",
  "Bm3235t_vegtissue"       = "#C7DA9F",
  "Bm3235t_stichidia"       = "#A2BE72",
  "Bm3235t_tetraspores"     = "#7D9E59"
)

p <- ggplot(sum_tbl, aes(x = stage_id, y = mean_log2_tpm1, group = 1)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#efefef", alpha = 0.7, linewidth = 0) +
  geom_line(color = "black", linewidth = 0.9) +
  geom_point(aes(color = stage_id), size = 2.6) +
  facet_wrap(~ phylorank, ncol = 2, nrow = 4, scales = "free_y") +
  labs(
    title = "Lifecycle Expression Trajectories by Phylostratum",
    subtitle = "Mean expression +/- SE across ordered stages (log2(TPM+1))",
    x = NULL,
    y = "Mean log2(TPM + 1)"
  ) +
  scale_color_manual(values = stage_cols_plot, drop = FALSE) +
  theme_bw(base_size = 10) +
  theme(
    aspect.ratio = 0.85,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "phylorank_lifecycle_expression_faceted.pdf"), p, width = 8.5, height = 14)
message("Wrote: ", file.path(out_dir, "phylorank_lifecycle_expression_faceted.pdf"))
message("Wrote: ", file.path(out_dir, "phylorank_lifecycle_expression_summary.tsv"))
