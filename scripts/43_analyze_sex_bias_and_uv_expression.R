#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("R/io.R")
source("R/qc.R")
source("R/de.R")
source("R/enrichment.R")
source("R/tai.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
expr <- prepare_expression_inputs(paths, design, analysis)

out_dir <- file.path(paths$output_root, "sex_bias")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

lfc_thr <- 1
fdr_thr <- 0.05
female_stages <- c("Bm3235f_vegtissue", "Bm3235f_procarps")
male_stages <- c("Bm3235m_vegtissue", "Bm3235m_spermatangia")
excluded_stages <- c("Bm3235m_spermatia")

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
stage_cols <- c(
  Bm3235f_vegtissue = "#8E2F3F",
  Bm3235f_procarps = "#5C1F2A",
  Bm3235m_vegtissue = "#D6AF3D",
  Bm3235m_spermatangia = "#B8860B",
  Bm3235m_spermatia = "#8C6A00",
  Bm3235f_pericarp = "#C77C7C",
  Bm3235f_cystocarps_imm = "#B79AC8",
  Bm3235f_cystocarps_mat = "#6E4B9E",
  Bm3235c_carposporophyte = "#3F1F78",
  Bm3235t_carpospores = "#4C2A85",
  Bm3235t_vegtissue = "#B5C46A",
  Bm3235t_stichidia = "#8FA33E",
  Bm3235t_tetraspores = "#6E8629"
)
sex_cols <- c(female_biased = "#8E2F3F", male_biased = "#B8860B", other_protein_coding = "#D9D2C8")
region_cols <- c(autosome = "#C8B78A", ChrUf = "#8E2F3F", ChrVm = "#5C1F2A", PAR = "#6A3D9A", other_unplaced = "#9E9E9E")

# --------- DE definition: pooled core female vs core male, excluding spermatia ---------
core_design <- design %>%
  mutate(
    sex_core_group = case_when(
      stage %in% female_stages ~ "female_core",
      stage %in% male_stages ~ "male_core",
      TRUE ~ NA_character_
    ),
    include = include & !stage %in% excluded_stages & !is.na(sex_core_group)
  )

core_samples <- core_design$sample_id[core_design$include]
counts_core <- expr$counts[, c("gene_id", core_samples), drop = FALSE]
contrast <- data.frame(
  contrast_id = "female_core_vs_male_core",
  factor = "sex_core_group",
  numerator = "female_core",
  denominator = "male_core",
  lfc_threshold = lfc_thr,
  fdr_threshold = fdr_thr,
  stringsAsFactors = FALSE
)

de <- run_differential_expression(counts_core, core_design, contrast)
de$gene_id <- normalize_gene_id(de$gene_id)
de$log2fc <- suppressWarnings(as.numeric(de$log2fc))
de$fdr <- suppressWarnings(as.numeric(de$fdr))

pc_genes <- get_protein_coding_gene_ids(paths)
de <- de %>%
  filter(gene_id %in% pc_genes) %>%
  mutate(
    sex_bias = case_when(
      is.finite(fdr) & is.finite(log2fc) & fdr <= fdr_thr & log2fc >= lfc_thr ~ "female_biased",
      is.finite(fdr) & is.finite(log2fc) & fdr <= fdr_thr & log2fc <= -lfc_thr ~ "male_biased",
      TRUE ~ "not_sex_biased"
    )
  )
write_tsv(de, file.path(out_dir, "sex_biased_core_de_results.tsv"))

sex_sets <- de %>%
  filter(sex_bias %in% c("female_biased", "male_biased")) %>%
  transmute(gene_id, set = sex_bias, log2fc, fdr, base_mean)
write_tsv(sex_sets, file.path(out_dir, "sex_biased_gene_sets.tsv"))

# --------- Coordinates, region labels, optional PAR ---------
coords_all <- read_coordinates(paths) %>% filter(gene_id %in% pc_genes)
organellar_chrs <- c("ChrMt", "ChrPt")
organellar_genes_excluded <- unique(coords_all$gene_id[coords_all$chr %in% organellar_chrs])
coords <- coords_all %>% filter(!chr %in% organellar_chrs)
region_df <- coords %>%
  transmute(gene_id, chr, start, end, region = case_when(
    chr == "ChrUf" ~ "ChrUf",
    chr == "ChrVm" ~ "ChrVm",
    grepl("^Chr[0-9]+$", chr) ~ "autosome",
    TRUE ~ "other_unplaced"
  ))

par_files <- c("config/par_regions.tsv", file.path(out_dir, "par_regions.tsv"))
par_file <- par_files[file.exists(par_files)][1]
par_present <- FALSE
par_n_regions <- 0L
par_n_genes <- 0L
if (!is.na(par_file) && nzchar(par_file)) {
  par <- read.delim(par_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("chr", "start", "end", "region_name") %in% colnames(par))) {
    par$start <- as.integer(par$start)
    par$end <- as.integer(par$end)
    par <- par[!is.na(par$chr) & !is.na(par$start) & !is.na(par$end) & par$end >= par$start, , drop = FALSE]
    par_present <- nrow(par) > 0
    par_n_regions <- nrow(par)
    if (par_present) {
      par_genes <- unique(unlist(lapply(seq_len(nrow(par)), function(i) {
        region_df$gene_id[region_df$chr == par$chr[i] & region_df$start <= par$end[i] & region_df$end >= par$start[i]]
      })))
      par_genes <- par_genes[!is.na(par_genes) & nzchar(par_genes)]
      par_n_genes <- length(par_genes)
      region_df$region[region_df$gene_id %in% par_genes] <- "PAR"
    }
  } else {
    warning("PAR file exists but lacks required columns: chr, start, end, region_name. Skipping PAR assignment.")
  }
}
region_df$region <- factor(region_df$region, levels = c("autosome", "ChrUf", "ChrVm", "PAR", "other_unplaced"))
region_cols_present <- region_cols[names(region_cols) %in% unique(as.character(region_df$region))]
write_tsv(region_df, file.path(out_dir, "sex_biased_gene_regions.tsv"))

universe <- unique(region_df$gene_id)
de_region <- de %>%
  filter(gene_id %in% universe) %>%
  left_join(region_df, by = "gene_id")
sex_sets <- sex_sets %>% filter(gene_id %in% universe)
write_tsv(sex_sets, file.path(out_dir, "sex_biased_gene_sets.tsv"))

all_gene_status <- region_df %>%
  left_join(de[, c("gene_id", "log2fc", "p_value", "fdr", "base_mean", "sex_bias")], by = "gene_id") %>%
  mutate(
    sex_bias = ifelse(is.na(sex_bias), "not_sex_biased", sex_bias),
    set = ifelse(sex_bias %in% c("female_biased", "male_biased"), sex_bias, "other_protein_coding")
  )
write_tsv(all_gene_status, file.path(out_dir, "sex_biased_all_gene_status.tsv"))

# --------- Region summaries and enrichment ---------
region_summary <- all_gene_status %>%
  count(set, region, name = "n_genes") %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup()
write_tsv(region_summary, file.path(out_dir, "sex_biased_region_summary.tsv"))

region_tests <- expand.grid(
  set = c("female_biased", "male_biased"),
  region = levels(region_df$region),
  stringsAsFactors = FALSE
) %>%
  filter(region %in% unique(as.character(region_df$region))) %>%
  rowwise() %>%
  do({
    s <- .$set
    r <- .$region
    in_set <- universe %in% de_region$gene_id[de_region$sex_bias == s]
    in_region <- region_df$region[match(universe, region_df$gene_id)] == r
    mat <- matrix(c(
      sum(in_set & in_region), sum(in_set & !in_region),
      sum(!in_set & in_region), sum(!in_set & !in_region)
    ), nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat, alternative = "two.sided")
    data.frame(
      set = s,
      region = r,
      n_set_genes = sum(in_set),
      n_set_region = sum(in_set & in_region),
      prop_set_region = sum(in_set & in_region) / max(1, sum(in_set)),
      n_background_region = sum(!in_set & in_region),
      prop_background_region = sum(!in_set & in_region) / max(1, sum(!in_set)),
      odds_ratio = unname(ft$estimate),
      log2_odds_ratio = log2(unname(ft$estimate)),
      direction = ifelse(unname(ft$estimate) >= 1, "enriched", "depleted"),
      fisher_p = ft$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(fisher_p, method = "BH")) %>%
  arrange(fdr, set, region)
write_tsv(region_tests, file.path(out_dir, "sex_biased_region_enrichment.tsv"))

expr_means <- compute_group_tpm_means(expr$tpm, design, group_col = "stage")
expr_means$gene_id <- normalize_gene_id(expr_means$gene_id)
mean_cols <- paste0("mean_tpm_", c(female_stages, male_stages))
expr_summary <- expr_means %>%
  filter(gene_id %in% universe) %>%
  select(gene_id, all_of(mean_cols)) %>%
  mutate(
    female_core_mean_tpm = rowMeans(across(all_of(paste0("mean_tpm_", female_stages))), na.rm = TRUE),
    male_core_mean_tpm = rowMeans(across(all_of(paste0("mean_tpm_", male_stages))), na.rm = TRUE)
  ) %>%
  select(gene_id, female_core_mean_tpm, male_core_mean_tpm) %>%
  left_join(all_gene_status[, c("gene_id", "sex_bias", "set", "log2fc", "fdr", "region", "chr")], by = "gene_id")
expr_region_summary <- expr_summary %>%
  group_by(set, region) %>%
  summarise(
    n_genes = n(),
    median_log2fc = median(log2fc, na.rm = TRUE),
    median_abs_log2fc = median(abs(log2fc), na.rm = TRUE),
    median_female_core_log2_tpm = median(log2(female_core_mean_tpm + 1), na.rm = TRUE),
    median_male_core_log2_tpm = median(log2(male_core_mean_tpm + 1), na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(expr_region_summary, file.path(out_dir, "sex_biased_expression_summary_by_region.tsv"))

sex_bias_expr_compare <- expr_summary %>%
  filter(sex_bias %in% c("female_biased", "male_biased")) %>%
  mutate(
    sex_bias_label = factor(
      recode(sex_bias, female_biased = "Female-biased genes", male_biased = "Male-biased genes"),
      levels = c("Female-biased genes", "Male-biased genes")
    ),
    male_core_log2_tpm = log2(male_core_mean_tpm + 1),
    female_core_log2_tpm = log2(female_core_mean_tpm + 1),
    m_vs_f_log2fc = male_core_log2_tpm - female_core_log2_tpm,
    f_vs_m_log2fc_tpm = female_core_log2_tpm - male_core_log2_tpm,
    reciprocal_class = case_when(
      female_core_mean_tpm >= 1 & male_core_mean_tpm < 1 ~ "female only",
      male_core_mean_tpm >= 1 & female_core_mean_tpm < 1 ~ "male only",
      female_core_mean_tpm >= 1 & male_core_mean_tpm >= 1 ~ "detected both",
      TRUE ~ "below TPM floor"
    )
  )
write_tsv(sex_bias_expr_compare, file.path(out_dir, "sex_biased_male_vs_female_expression_values.tsv"))

max_expr_axis <- max(c(sex_bias_expr_compare$male_core_log2_tpm, sex_bias_expr_compare$female_core_log2_tpm), na.rm = TRUE)
if (!is.finite(max_expr_axis) || max_expr_axis <= 0) max_expr_axis <- 1
p_mf_scatter <- ggplot(sex_bias_expr_compare, aes(x = male_core_log2_tpm, y = female_core_log2_tpm, color = sex_bias_label)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.42, size = 1.05, stroke = 0) +
  scale_color_manual(values = c("Female-biased genes" = "#8E2F3F", "Male-biased genes" = "#B8860B"), name = NULL) +
  coord_equal(xlim = c(0, max_expr_axis * 1.03), ylim = c(0, max_expr_axis * 1.03), expand = FALSE) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top") +
  labs(
    x = "Male core expression log2(TPM + 1)",
    y = "Female core expression log2(TPM + 1)",
    title = "Sex-biased genes: male versus female core expression",
    subtitle = "Dashed line indicates equal expression; points near axes indicate target-sex-restricted expression"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_male_vs_female_expression_scatter.pdf"), p_mf_scatter, width = 5.2, height = 5.5, useDingbats = FALSE)

p_mf_lfc <- ggplot(sex_bias_expr_compare, aes(x = sex_bias_label, y = m_vs_f_log2fc, fill = sex_bias_label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_violin(width = 0.72, trim = TRUE, color = NA, alpha = 0.86) +
  geom_boxplot(width = 0.16, outlier.shape = NA, color = "grey20", fill = "white", linewidth = 0.25) +
  scale_fill_manual(values = c("Female-biased genes" = "#8E2F3F", "Male-biased genes" = "#B8860B"), guide = "none") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(
    x = NULL,
    y = "Male - female core log2(TPM + 1)",
    title = "Direction of sex-biased expression",
    subtitle = "Positive values indicate male-skewed expression; negative values indicate female-skewed expression"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_male_vs_female_log2fc_violin.pdf"), p_mf_lfc, width = 4.8, height = 4.2, useDingbats = FALSE)

p_mf_density <- ggplot(sex_bias_expr_compare, aes(x = m_vs_f_log2fc, fill = sex_bias_label, color = sex_bias_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey35") +
  geom_density(alpha = 0.28, linewidth = 0.7, adjust = 1.0) +
  scale_fill_manual(values = c("Female-biased genes" = "#8E2F3F", "Male-biased genes" = "#B8860B"), name = NULL) +
  scale_color_manual(values = c("Female-biased genes" = "#8E2F3F", "Male-biased genes" = "#B8860B"), name = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top") +
  labs(
    x = "Male - female core log2(TPM + 1)",
    y = "Density",
    title = "Reciprocal sex-core expression bias",
    subtitle = "Negative values indicate female-skewed expression; positive values indicate male-skewed expression"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_male_vs_female_log2fc_density.pdf"), p_mf_density, width = 5.6, height = 4.0, useDingbats = FALSE)

p_mf_classes <- sex_bias_expr_compare %>%
  count(sex_bias_label, reciprocal_class, name = "n_genes") %>%
  group_by(sex_bias_label) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup() %>%
  filter(reciprocal_class != "below TPM floor") %>%
  mutate(
    reciprocal_class = factor(reciprocal_class, levels = c("detected both", "female only", "male only")),
    label = sprintf("%s", n_genes)
  ) %>%
  ggplot(aes(x = reciprocal_class, y = n_genes, fill = sex_bias_label)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 3.0, color = "grey15") +
  facet_wrap(~ sex_bias_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c("Female-biased genes" = "#8E2F3F", "Male-biased genes" = "#B8860B"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_blank()
  ) +
  labs(x = NULL, y = "Number of genes", title = "Sex-biased genes by reciprocal expression class")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_reciprocal_expression_class_barplot.pdf"), p_mf_classes, width = 6.2, height = 4.2, useDingbats = FALSE)

sex_bias_fc_sector <- sex_bias_expr_compare %>%
  mutate(
    deseq_m_vs_f_log2fc = -log2fc,
    fc_sector = ifelse(deseq_m_vs_f_log2fc >= 0, "M > F", "F > M"),
    sex_bias_label = factor(sex_bias_label, levels = c("Female-biased genes", "Male-biased genes"))
  )
write_tsv(sex_bias_fc_sector, file.path(out_dir, "sex_biased_m_vs_f_log2fc_values.tsv"))

fc_lim <- max(abs(sex_bias_fc_sector$deseq_m_vs_f_log2fc[is.finite(sex_bias_fc_sector$deseq_m_vs_f_log2fc)]), na.rm = TRUE)
if (!is.finite(fc_lim) || fc_lim <= 0) fc_lim <- 1
fc_lim <- ceiling(fc_lim)

p_fc_sector <- ggplot(sex_bias_fc_sector, aes(x = deseq_m_vs_f_log2fc, fill = fc_sector)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#8E2F3F", alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#B8860B", alpha = 0.06) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey25") +
  geom_histogram(binwidth = 0.25, boundary = 0, color = "white", linewidth = 0.15) +
  facet_wrap(~ sex_bias_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("F > M" = "#8E2F3F", "M > F" = "#B8860B"), name = NULL) +
  scale_x_continuous(limits = c(-fc_lim, fc_lim), breaks = pretty(c(-fc_lim, fc_lim), n = 7)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "M vs F log2 fold change",
    y = "Number of genes",
    title = "Sex-biased genes separated by M-vs-F log2 fold change sector",
    subtitle = "Negative sector = higher in female core; positive sector = higher in male core"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_m_vs_f_log2fc_sector_histogram.pdf"), p_fc_sector, width = 5.4, height = 6.2, useDingbats = FALSE)

p_fc_sector_density <- ggplot(sex_bias_fc_sector, aes(x = deseq_m_vs_f_log2fc, y = after_stat(count), fill = fc_sector)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#8E2F3F", alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#B8860B", alpha = 0.06) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey25") +
  geom_density(alpha = 0.72, linewidth = 0.25, adjust = 0.9, color = "grey20") +
  facet_wrap(~ sex_bias_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("F > M" = "#8E2F3F", "M > F" = "#B8860B"), name = NULL) +
  scale_x_continuous(limits = c(-fc_lim, fc_lim), breaks = pretty(c(-fc_lim, fc_lim), n = 7)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "M vs F log2 fold change",
    y = "Density scaled to gene count",
    title = "Sex-biased genes separated by M-vs-F log2 fold change sector",
    subtitle = "Negative sector = higher in female core; positive sector = higher in male core"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_m_vs_f_log2fc_sector_density.pdf"), p_fc_sector_density, width = 5.4, height = 6.2, useDingbats = FALSE)

# --------- Triphasy overlap ---------
trip_files <- c(
  file.path(paths$output_root, "generation", "triphasy_biased_gene_set_membership.tsv"),
  file.path(paths$output_root, "generation", "triphasy_downregulated_gene_set_membership.tsv")
)
trip <- bind_rows(lapply(trip_files[file.exists(trip_files)], function(f) read.delim(f, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))) %>%
  mutate(gene_id = normalize_gene_id(gene_id)) %>%
  filter(gene_id %in% universe) %>%
  distinct(gene_id, set)
trip_overlap <- expand.grid(sex_set = c("female_biased", "male_biased"), triphasy_set = sort(unique(trip$set)), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  do({
    sex_genes <- unique(de_region$gene_id[de_region$sex_bias == .$sex_set])
    tri_genes <- unique(trip$gene_id[trip$set == .$triphasy_set])
    in_sex <- universe %in% sex_genes
    in_tri <- universe %in% tri_genes
    mat <- matrix(c(sum(in_sex & in_tri), sum(in_sex & !in_tri), sum(!in_sex & in_tri), sum(!in_sex & !in_tri)), nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat, alternative = "two.sided")
    data.frame(
      sex_set = .$sex_set,
      triphasy_set = .$triphasy_set,
      n_sex_set = sum(in_sex),
      n_triphasy_set = sum(in_tri),
      n_overlap = sum(in_sex & in_tri),
      jaccard = sum(in_sex & in_tri) / max(1, sum(in_sex | in_tri)),
      odds_ratio = unname(ft$estimate),
      log2_odds_ratio = log2(unname(ft$estimate)),
      fisher_p = ft$p.value,
      overlap_gene_ids = paste(sort(intersect(sex_genes, tri_genes)), collapse = ","),
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(fisher_p, method = "BH")) %>%
  arrange(fdr, sex_set, triphasy_set)
write_tsv(trip_overlap, file.path(out_dir, "sex_biased_triphasy_overlap.tsv"))

# --------- TAP enrichment ---------
anno <- read_functional_annotation(paths) %>%
  transmute(gene_id = normalize_gene_id(gene_id), TAP_family = as.character(TAP_family), GOs, KEGG_ko, PFAMs, COG_category) %>%
  filter(gene_id %in% universe) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(
    TAP_family = ifelse(is.na(TAP_family) | !nzchar(trimws(TAP_family)) | trimws(TAP_family) %in% c("-", "NA"), "non_TAP", TAP_family),
    tap_status = ifelse(TAP_family == "non_TAP", "non-TAP", "TAP")
  )

tap_summary <- all_gene_status %>%
  select(gene_id, set) %>%
  left_join(anno[, c("gene_id", "tap_status")], by = "gene_id") %>%
  count(set, tap_status, name = "n_genes") %>%
  complete(set = c("female_biased", "male_biased", "other_protein_coding"), tap_status = c("TAP", "non-TAP"), fill = list(n_genes = 0L)) %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup()
write_tsv(tap_summary, file.path(out_dir, "sex_biased_tap_proportion_summary.tsv"))

tap_tests <- lapply(c("female_biased", "male_biased"), function(s) {
  in_set <- universe %in% de_region$gene_id[de_region$sex_bias == s]
  is_tap <- anno$tap_status[match(universe, anno$gene_id)] == "TAP"
  is_tap[is.na(is_tap)] <- FALSE
  mat <- matrix(c(sum(in_set & is_tap), sum(in_set & !is_tap), sum(!in_set & is_tap), sum(!in_set & !is_tap)), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "two.sided")
  data.frame(
    set = s,
    n_set_genes = sum(in_set),
    n_set_taps = sum(in_set & is_tap),
    prop_set_taps = sum(in_set & is_tap) / max(1, sum(in_set)),
    n_background_taps = sum(!in_set & is_tap),
    prop_background_taps = sum(!in_set & is_tap) / max(1, sum(!in_set)),
    odds_ratio = unname(ft$estimate),
    log2_odds_ratio = log2(unname(ft$estimate)),
    direction = ifelse(unname(ft$estimate) >= 1, "enriched", "depleted"),
    fisher_p = ft$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>% mutate(fdr = p.adjust(fisher_p, method = "BH"))
write_tsv(tap_tests, file.path(out_dir, "sex_biased_tap_enrichment.tsv"))

# --------- Phylorank distribution/enrichment ---------
age <- read_gene_age_table(paths$gene_age_tsv) %>%
  mutate(
    gene_id = normalize_gene_id(gene_id),
    rank = suppressWarnings(as.integer(rank)),
    rank_label = ifelse(!is.na(rank), paste0("PS", rank, ": ", taxon), as.character(taxon))
  ) %>%
  filter(gene_id %in% universe, !is.na(rank_label), nzchar(rank_label)) %>%
  distinct(gene_id, .keep_all = TRUE)
rank_levels <- age %>% distinct(rank, rank_label) %>% arrange(rank) %>% pull(rank_label)
rank_dist <- all_gene_status %>%
  select(gene_id, set) %>%
  inner_join(age[, c("gene_id", "rank_label")], by = "gene_id") %>%
  count(set, rank_label, name = "value")
write_tsv(rank_dist, file.path(out_dir, "sex_biased_phylorank_distribution.tsv"))

age_universe <- intersect(universe, age$gene_id)
rank_grid <- expand.grid(set = c("female_biased", "male_biased"), rank_label = rank_levels, stringsAsFactors = FALSE)
rank_rows <- lapply(seq_len(nrow(rank_grid)), function(i) {
  s <- rank_grid$set[i]
  rr <- rank_grid$rank_label[i]
  in_set <- age_universe %in% de_region$gene_id[de_region$sex_bias == s]
  in_rank <- age$rank_label[match(age_universe, age$gene_id)] == rr
  mat <- matrix(c(sum(in_set & in_rank), sum(in_set & !in_rank), sum(!in_set & in_rank), sum(!in_set & !in_rank)), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "two.sided")
  data.frame(set = s, rank_label = rr, observed = sum(in_set & in_rank), set_total = sum(in_set), background_total_rank = sum(in_rank), background_total = length(age_universe), odds_ratio = unname(ft$estimate), log2_odds_ratio = log2(unname(ft$estimate)), direction = ifelse(unname(ft$estimate) >= 1, "enriched", "depleted"), fisher_p = ft$p.value, stringsAsFactors = FALSE)
})
rank_tests <- bind_rows(rank_rows) %>%
  group_by(set) %>% mutate(fdr_within_set = p.adjust(fisher_p, method = "BH")) %>% ungroup() %>%
  mutate(fdr_global = p.adjust(fisher_p, method = "BH"))
write_tsv(rank_tests, file.path(out_dir, "sex_biased_phylorank_enrichment.tsv"))

# --------- Functional enrichment for sex sets ---------
func_df <- data.frame(gene_id = pc_genes, stringsAsFactors = FALSE) %>%
  left_join(anno, by = "gene_id") %>%
  mutate(marker_group = ifelse(gene_id %in% de_region$gene_id[de_region$sex_bias == "female_biased"], "female_biased", ifelse(gene_id %in% de_region$gene_id[de_region$sex_bias == "male_biased"], "male_biased", "background")), marker_flag = marker_group %in% c("female_biased", "male_biased"), marker_class = "sex_biased")
func_enrich <- run_marker_functional_enrichment(func_df, marker_classes = NULL, universe_genes = pc_genes)
if (nrow(func_enrich) > 0) {
  func_enrich <- func_enrich %>% rename(set = marker_group)
  func_enrich <- annotate_go_term_names(func_enrich)
  func_enrich <- drop_obsolete_go_terms(func_enrich)
}
write_tsv(func_enrich, file.path(out_dir, "sex_biased_functional_enrichment.tsv"))
write_tsv(func_enrich %>% filter(fdr <= 0.05) %>% arrange(fdr, p_value) %>% group_by(set, source) %>% slice_head(n = 10) %>% ungroup(), file.path(out_dir, "sex_biased_functional_enrichment_top.tsv"))

# --------- Plots ---------
size_df <- all_gene_status %>%
  count(set, name = "n_genes") %>%
  mutate(set = factor(set, levels = c("female_biased", "male_biased", "other_protein_coding"))) %>%
  arrange(set) %>%
  mutate(ymax = cumsum(n_genes), ymin = lag(ymax, default = 0), label_pos = (ymin + ymax)/2, label = sprintf("%s\n%s", set, n_genes))
write_tsv(size_df, file.path(out_dir, "sex_biased_gene_set_sizes_with_other.tsv"))
p_donut <- ggplot(size_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.6, fill = set)) +
  geom_rect(color = "white", linewidth = 0.25) +
  geom_text(aes(x = 4.35, y = label_pos, label = label), size = 2.8, hjust = 0) +
  coord_polar(theta = "y") +
  xlim(c(0, 5.4)) +
  scale_fill_manual(values = sex_cols, drop = FALSE) +
  theme_void(base_size = 10) +
  theme(legend.position = "none")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_gene_set_donut.pdf"), p_donut, width = 5.4, height = 4.4, useDingbats = FALSE)

p_reg <- ggplot(region_tests, aes(x = region, y = log2_odds_ratio, fill = set)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = NA) +
  geom_text(aes(label = ifelse(fdr <= 0.05, "*", "")), position = position_dodge(width = 0.72), vjust = ifelse(region_tests$log2_odds_ratio >= 0, -0.3, 1.2), size = 4) +
  scale_fill_manual(values = sex_cols[c("female_biased", "male_biased")], drop = FALSE) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank()) +
  labs(x = NULL, y = "log2 odds ratio\n(region enrichment)")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_region_enrichment_barplot.pdf"), p_reg, width = 6.4, height = 4.4, useDingbats = FALSE)

region_heat <- region_tests %>%
  mutate(
    set_label = factor(unname(c(female_biased = "Female-biased genes", male_biased = "Male-biased genes")[set]), levels = rev(c("Female-biased genes", "Male-biased genes"))),
    region = factor(region, levels = c("autosome", "ChrUf", "ChrVm", "PAR")),
    log2_or_label = sprintf("%.2f", log2_odds_ratio),
    fdr_sig = is.finite(fdr) & fdr <= 0.05
  ) %>%
  filter(!is.na(region))
region_vals <- region_heat$log2_odds_ratio[is.finite(region_heat$log2_odds_ratio)]
region_min <- min(region_vals, na.rm = TRUE)
region_max <- max(region_vals, na.rm = TRUE)
if (!is.finite(region_min) || !is.finite(region_max) || region_min == region_max) {
  region_min <- -1
  region_max <- 1
}
region_mid_pos <- if (region_max > 0) region_max * 0.45 else region_max
p_reg_heat <- ggplot(region_heat, aes(x = region, y = set_label, fill = log2_odds_ratio)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = log2_or_label, fontface = if_else(fdr_sig, "bold", "plain")), size = 2.8, color = "black") +
  scale_fill_gradientn(
    colors = c("#6B6B6B", "white", "#CC9966", "#663333"),
    values = scales::rescale(c(region_min, 0, region_mid_pos, region_max), from = c(region_min, region_max)),
    limits = c(region_min, region_max),
    oob = scales::squish,
    na.value = "#D9D9D9",
    name = "log2 odds\nratio"
  ) +
  labs(
    title = "Sex-biased gene regional enrichment",
    subtitle = "Positive = enriched; negative = depleted; bold labels indicate Fisher FDR <= 0.05",
    x = NULL,
    y = NULL
  ) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_region_enrichment_log2_oddsratio_heatmap.pdf"), p_reg_heat, width = 5.8, height = 3.4, useDingbats = FALSE)

p_stack <- ggplot(region_summary %>% filter(set %in% c("female_biased", "male_biased", "other_protein_coding")), aes(x = set, y = proportion, fill = region)) +
  geom_col(position = "fill", width = 0.72, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = region_cols_present, drop = TRUE) +
  scale_y_continuous(labels = function(x) sprintf("%d%%", round(100*x))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank()) +
  labs(x = NULL, y = "Gene distribution")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_region_distribution_stackedbar.pdf"), p_stack, width = 6.4, height = 4.5, useDingbats = FALSE)

p_box <- ggplot(de_region %>% filter(is.finite(log2fc), region %in% c("autosome", "ChrUf", "ChrVm", "PAR")), aes(x = region, y = log2fc, fill = region)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, width = 0.62, color = "grey20", linewidth = 0.25) +
  coord_cartesian(ylim = quantile(de_region$log2fc[is.finite(de_region$log2fc)], c(0.01, 0.99), na.rm = TRUE)) +
  scale_fill_manual(values = region_cols_present, drop = TRUE) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = NULL, y = "Female vs male core log2FC")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_bias_log2fc_by_region_boxplot.pdf"), p_box, width = 5.4, height = 4.3, useDingbats = FALSE)

hm_genes <- de_region %>% filter(sex_bias %in% c("female_biased", "male_biased")) %>% arrange(sex_bias, desc(abs(log2fc))) %>% pull(gene_id)
mean_cols_all <- grep("^mean_tpm_", colnames(expr_means), value = TRUE)
tai_order_file <- file.path(paths$output_root, "tai", "tai_plot_input.tsv")
if (file.exists(tai_order_file)) {
  tai_order <- read.delim(tai_order_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("stage_order", "stage") %in% colnames(tai_order))) {
    tai_order <- tai_order[order(tai_order$stage_order), , drop = FALSE]
    mean_cols <- intersect(as.character(tai_order$stage), mean_cols_all)
    mean_cols <- c(mean_cols, setdiff(mean_cols_all, mean_cols))
  } else {
    mean_cols <- order_stage_columns(mean_cols_all)
  }
} else {
  mean_cols <- order_stage_columns(mean_cols_all)
}
stage_ids <- sub("^mean_tpm_", "", mean_cols)
stage_cols_present <- stage_cols[stage_ids]
missing_stage_cols <- setdiff(stage_ids, names(stage_cols_present))
if (length(missing_stage_cols) > 0) {
  fallback_cols <- rep("#999999", length(missing_stage_cols))
  names(fallback_cols) <- missing_stage_cols
  stage_cols_present <- c(stage_cols_present, fallback_cols)
}
stage_cols_present <- stage_cols_present[stage_ids]

hm_expr <- expr_means %>% filter(gene_id %in% hm_genes) %>% select(gene_id, all_of(mean_cols))
hm_mat <- as.matrix(hm_expr[, mean_cols, drop = FALSE])
storage.mode(hm_mat) <- "numeric"
hm_mat <- log2(hm_mat + 1)
rownames(hm_mat) <- hm_expr$gene_id
hm_mat <- hm_mat[hm_genes[hm_genes %in% rownames(hm_mat)], , drop = FALSE]
hm_z <- t(scale(t(hm_mat)))
hm_z[!is.finite(hm_z)] <- 0
colnames(hm_z) <- stage_ids

row_split <- factor(de_region$sex_bias[match(rownames(hm_z), de_region$gene_id)], levels = c("female_biased", "male_biased"), labels = c("Female-biased", "Male-biased"))
row_region <- factor(de_region$region[match(rownames(hm_z), de_region$gene_id)], levels = names(region_cols_present))

heatmap_values <- data.frame(gene_id = rownames(hm_z), hm_z, check.names = FALSE)
heatmap_values <- cbind(
  heatmap_values,
  sex_bias = as.character(de_region$sex_bias[match(rownames(hm_z), de_region$gene_id)]),
  region = as.character(de_region$region[match(rownames(hm_z), de_region$gene_id)])
)
write_tsv(heatmap_values, file.path(out_dir, "sex_biased_expression_heatmap_values.tsv"))

pdf(file.path(fig_dir, "sex_biased_expression_heatmap.pdf"), width = 8.5, height = 10.5)
ht <- Heatmap(
  hm_z,
  name = "row z",
  col = circlize::colorRamp2(c(-2, 0, 2), c("#cccccc", "white", "#663333")),
  top_annotation = HeatmapAnnotation(
    stage = stage_ids,
    col = list(stage = stage_cols_present),
    show_annotation_name = FALSE,
    annotation_legend_param = list(stage = list(ncol = 1))
  ),
  left_annotation = rowAnnotation(
    sex_bias = row_split,
    region = row_region,
    col = list(sex_bias = c("Female-biased" = sex_cols[["female_biased"]], "Male-biased" = sex_cols[["male_biased"]]), region = region_cols_present),
    show_annotation_name = FALSE
  ),
  row_split = row_split,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = stage_ids,
  column_names_rot = 45,
  use_raster = TRUE,
  raster_quality = 2,
  heatmap_legend_param = list(title = "row z")
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

p_tap <- tap_summary %>% filter(tap_status == "TAP") %>% mutate(set = factor(set, levels = c("female_biased", "male_biased", "other_protein_coding")), label = sprintf("%.1f%%\n%s/%s", 100*proportion, n_genes, total_genes)) %>%
  ggplot(aes(x = set, y = proportion, fill = set)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = sex_cols, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) sprintf("%d%%", round(100*x)),
    breaks = seq(0, 0.15, by = 0.025),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.10))
  ) +
  coord_cartesian(ylim = c(0, 0.12), clip = "off") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = NULL, y = "TAP genes (%)")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_tap_proportion_barplot.pdf"), p_tap, width = 5.4, height = 4.0, useDingbats = FALSE)

p_rank <- rank_dist %>% mutate(rank_label = factor(rank_label, levels = rank_levels), set = factor(set, levels = c("female_biased", "male_biased", "other_protein_coding"))) %>%
  ggplot(aes(x = set, y = value, fill = rank_label)) +
  geom_col(position = "fill", width = 0.72, color = "white", linewidth = 0.1) +
  scale_y_continuous(labels = function(x) sprintf("%d%%", round(100*x))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "right", legend.title = element_blank()) +
  labs(x = NULL, y = "Phylorank distribution")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_phylorank_stacked_bar_relative.pdf"), p_rank, width = 7.6, height = 4.8, useDingbats = FALSE)


# --------- True sex-specific expression among FBGs/MBGs ---------
# Specificity is split by genomic compartment so UV/PAR-linked absence does not get
# conflated with genuinely autosomal sex-specific expression.
all_stage_mean_cols <- intersect(paste0("mean_tpm_", stage_order), colnames(expr_means))
gametophyte_stages <- c(female_stages, male_stages, excluded_stages)
non_gametophyte_stages <- setdiff(stage_order, gametophyte_stages)
non_gametophyte_cols <- intersect(paste0("mean_tpm_", non_gametophyte_stages), colnames(expr_means))
female_outside_core_cols <- setdiff(all_stage_mean_cols, paste0("mean_tpm_", female_stages))
male_outside_core_cols <- setdiff(all_stage_mean_cols, paste0("mean_tpm_", male_stages))

sex_specificity <- sex_sets %>%
  select(gene_id, set) %>%
  left_join(region_df[, c("gene_id", "chr", "region")], by = "gene_id") %>%
  left_join(expr_means %>% select(gene_id, all_of(all_stage_mean_cols)), by = "gene_id") %>%
  rowwise() %>%
  mutate(
    female_core_max_tpm = max(c_across(all_of(paste0("mean_tpm_", female_stages))), na.rm = TRUE),
    male_core_max_tpm = max(c_across(all_of(paste0("mean_tpm_", male_stages))), na.rm = TRUE),
    non_gametophyte_max_tpm = if (length(non_gametophyte_cols) > 0) max(c_across(all_of(non_gametophyte_cols)), na.rm = TRUE) else NA_real_,
    female_outside_core_max_tpm = if (length(female_outside_core_cols) > 0) max(c_across(all_of(female_outside_core_cols)), na.rm = TRUE) else NA_real_,
    male_outside_core_max_tpm = if (length(male_outside_core_cols) > 0) max(c_across(all_of(male_outside_core_cols)), na.rm = TRUE) else NA_real_,
    outside_target_core_max_tpm = ifelse(set == "female_biased", female_outside_core_max_tpm, male_outside_core_max_tpm),
    target_max_tpm = ifelse(set == "female_biased", female_core_max_tpm, male_core_max_tpm),
    opposing_max_tpm_raw = ifelse(set == "female_biased", male_core_max_tpm, female_core_max_tpm),
    sdr_expected_hemizygous = (set == "female_biased" & as.character(region) == "ChrUf") |
      (set == "male_biased" & as.character(region) == "ChrVm"),
    opposing_max_tpm = ifelse(sdr_expected_hemizygous, NA_real_, opposing_max_tpm_raw),
    target_detected = is.finite(target_max_tpm) & target_max_tpm >= 1,
    opposing_detected_raw = is.finite(opposing_max_tpm_raw) & opposing_max_tpm_raw >= 1,
    opposing_detected = is.finite(opposing_max_tpm) & opposing_max_tpm >= 1,
    genomic_compartment = case_when(
      as.character(region) == "autosome" ~ "autosomal",
      as.character(region) %in% c("ChrUf", "ChrVm") ~ "SDR_nonPAR",
      as.character(region) == "PAR" ~ "PAR",
      TRUE ~ "other_unplaced"
    ),
    specificity_class = case_when(
      target_detected & !opposing_detected & genomic_compartment == "autosomal" ~ "autosomal_target_detectable_only",
      target_detected & !opposing_detected & genomic_compartment == "SDR_nonPAR" ~ "SDR_target_detectable_only",
      target_detected & !opposing_detected & genomic_compartment == "PAR" ~ "PAR_target_detectable_only",
      target_detected & !opposing_detected ~ "other_target_detectable_only",
      target_detected & opposing_detected ~ "target_biased_but_detected_elsewhere",
      TRUE ~ "not_detected_at_floor"
    ),
    non_gametophyte_detected = is.finite(non_gametophyte_max_tpm) & non_gametophyte_max_tpm >= 1,
    outside_target_core_detected = is.finite(outside_target_core_max_tpm) & outside_target_core_max_tpm >= 1,
    expression_specificity = case_when(
      target_detected & !opposing_detected ~ "target_detectable_only",
      target_detected & opposing_detected ~ "target_biased_but_detected_elsewhere",
      TRUE ~ "not_detected_at_floor"
    ),
    lifecycle_specificity_class = case_when(
      !target_detected ~ "not_detected_at_floor",
      target_detected & opposing_detected ~ "sex_biased_detected_in_both_sex_cores",
      target_detected & !opposing_detected & !non_gametophyte_detected ~ "sex_core_target_only_gametophyte_restricted",
      target_detected & !opposing_detected & non_gametophyte_detected ~ "sex_core_target_only_but_lifecycle_detected",
      TRUE ~ "unclassified"
    ),
    target_core_lifecycle_specific = target_detected & !outside_target_core_detected
  ) %>%
  ungroup()
write_tsv(sex_specificity, file.path(out_dir, "sex_biased_target_specificity.tsv"))

sex_specificity_target_detected_audit <- sex_specificity %>%
  group_by(set) %>%
  summarise(
    total_de_defined_sbg = n(),
    target_detected_sbg = sum(target_detected, na.rm = TRUE),
    dropped_below_target_tpm_floor = sum(!target_detected, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(sex_specificity_target_detected_audit, file.path(out_dir, "sex_biased_target_detected_audit.tsv"))

# Refined SBGs for downstream summaries: keep only DE-defined SBGs that are
# detectably expressed in their target sex core. This matches the specificity
# plots and avoids carrying below-floor DE calls into TAP/phylorank summaries.
sex_sets_target_detected <- sex_specificity %>%
  filter(target_detected) %>%
  transmute(gene_id, set) %>%
  left_join(sex_sets[, c("gene_id", "log2fc", "fdr", "base_mean")], by = "gene_id") %>%
  distinct(gene_id, set, .keep_all = TRUE)

all_gene_status_target_detected <- all_gene_status %>%
  mutate(set = case_when(
    gene_id %in% sex_sets_target_detected$gene_id[sex_sets_target_detected$set == "female_biased"] ~ "female_biased",
    gene_id %in% sex_sets_target_detected$gene_id[sex_sets_target_detected$set == "male_biased"] ~ "male_biased",
    TRUE ~ "other_protein_coding"
  ))

tap_summary <- all_gene_status_target_detected %>%
  select(gene_id, set) %>%
  left_join(anno[, c("gene_id", "tap_status")], by = "gene_id") %>%
  mutate(tap_status = ifelse(is.na(tap_status), "non-TAP", tap_status)) %>%
  count(set, tap_status, name = "n_genes") %>%
  complete(set = c("female_biased", "male_biased", "other_protein_coding"), tap_status = c("TAP", "non-TAP"), fill = list(n_genes = 0L)) %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup()
write_tsv(tap_summary, file.path(out_dir, "sex_biased_tap_proportion_summary.tsv"))

tap_tests <- lapply(c("female_biased", "male_biased"), function(s) {
  in_set <- universe %in% sex_sets_target_detected$gene_id[sex_sets_target_detected$set == s]
  is_tap <- anno$tap_status[match(universe, anno$gene_id)] == "TAP"
  is_tap[is.na(is_tap)] <- FALSE
  mat <- matrix(c(sum(in_set & is_tap), sum(in_set & !is_tap), sum(!in_set & is_tap), sum(!in_set & !is_tap)), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "two.sided")
  data.frame(
    set = s,
    n_set_genes = sum(in_set),
    n_set_taps = sum(in_set & is_tap),
    prop_set_taps = sum(in_set & is_tap) / max(1, sum(in_set)),
    n_background_taps = sum(!in_set & is_tap),
    prop_background_taps = sum(!in_set & is_tap) / max(1, sum(!in_set)),
    odds_ratio = unname(ft$estimate),
    log2_odds_ratio = log2(unname(ft$estimate)),
    direction = ifelse(unname(ft$estimate) >= 1, "enriched", "depleted"),
    fisher_p = ft$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>% mutate(fdr = p.adjust(fisher_p, method = "BH"))
write_tsv(tap_tests, file.path(out_dir, "sex_biased_tap_enrichment.tsv"))

p_tap <- tap_summary %>%
  filter(tap_status == "TAP") %>%
  mutate(
    set = factor(set, levels = c("female_biased", "male_biased", "other_protein_coding")),
    label = sprintf("%.1f%%\n%s/%s", 100 * proportion, n_genes, total_genes)
  ) %>%
  ggplot(aes(x = set, y = proportion, fill = set)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = sex_cols, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) sprintf("%d%%", round(100*x)),
    breaks = seq(0, 0.10, by = 0.025),
    limits = c(0, 0.10),
    expand = expansion(mult = c(0, 0.06))
  ) +
  coord_cartesian(ylim = c(0, 0.10), clip = "off") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = NULL, y = "TAP genes (%)")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_tap_proportion_barplot.pdf"), p_tap, width = 5.4, height = 4.0, useDingbats = FALSE)

sex_specificity_summary <- sex_specificity %>%
  filter(target_detected) %>%
  count(set, specificity_class, name = "n_genes") %>%
  complete(
    set = c("female_biased", "male_biased"),
    specificity_class = c(
      "target_biased_but_detected_elsewhere",
      "other_target_detectable_only",
      "PAR_target_detectable_only",
      "SDR_target_detectable_only",
      "autosomal_target_detectable_only"
    ),
    fill = list(n_genes = 0L)
  ) %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()
write_tsv(sex_specificity_summary, file.path(out_dir, "sex_biased_target_specificity_summary.tsv"))

sex_specificity_region_summary <- sex_specificity %>%
  filter(target_detected) %>%
  count(set, genomic_compartment, expression_specificity, name = "n_genes") %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()
write_tsv(sex_specificity_region_summary, file.path(out_dir, "sex_biased_target_specificity_by_region_summary.tsv"))

sex_lifecycle_specificity_summary <- sex_specificity %>%
  filter(target_detected) %>%
  count(set, lifecycle_specificity_class, name = "n_genes") %>%
  complete(
    set = c("female_biased", "male_biased"),
    lifecycle_specificity_class = c(
      "sex_biased_detected_in_both_sex_cores",
      "sex_core_target_only_but_lifecycle_detected",
      "sex_core_target_only_gametophyte_restricted"
    ),
    fill = list(n_genes = 0L)
  ) %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()
write_tsv(sex_lifecycle_specificity_summary, file.path(out_dir, "sex_biased_lifecycle_specificity_summary.tsv"))

sex_lifecycle_specificity_by_region <- sex_specificity %>%
  filter(target_detected) %>%
  count(set, genomic_compartment, lifecycle_specificity_class, name = "n_genes") %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()
write_tsv(sex_lifecycle_specificity_by_region, file.path(out_dir, "sex_biased_lifecycle_specificity_by_region.tsv"))

spec_cols <- c(
  autosomal_target_detectable_only = "#663333",
  SDR_target_detectable_only = "#8E2F3F",
  PAR_target_detectable_only = "#6A3D9A",
  other_target_detectable_only = "#6B6B6B",
  target_biased_but_detected_elsewhere = "#D9D2C8",
  not_detected_at_floor = "#9E9E9E"
)
spec_labels <- c(
  female_biased = "Female-biased genes",
  male_biased = "Male-biased genes"
)
spec_order <- c(
  "target_biased_but_detected_elsewhere",
  "other_target_detectable_only",
  "PAR_target_detectable_only",
  "SDR_target_detectable_only",
  "autosomal_target_detectable_only"
)
p_spec <- sex_specificity_summary %>%
  filter(n_genes > 0) %>%
  mutate(
    set_label = factor(unname(spec_labels[set]), levels = unname(spec_labels[c("female_biased", "male_biased")])),
    specificity_class = factor(specificity_class, levels = spec_order),
    class_label = recode(
      as.character(specificity_class),
      autosomal_target_detectable_only = "Autosomal only",
      SDR_target_detectable_only = "SDR only",
      PAR_target_detectable_only = "PAR only",
      other_target_detectable_only = "Other/unplaced only",
      target_biased_but_detected_elsewhere = "Detected in both sex cores",
      not_detected_at_floor = "Below TPM floor in target core"
    ),
    class_label = factor(class_label, levels = c(
      "Autosomal only",
      "SDR only",
      "PAR only",
      "Other/unplaced only",
      "Detected in both sex cores",
      "Below TPM floor in target core"
    )),
    label = sprintf("%s", n_genes)
  ) %>%
  ggplot(aes(x = class_label, y = n_genes, fill = specificity_class)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 3.0, color = "grey15") +
  facet_wrap(~ set_label, nrow = 1) +
  scale_fill_manual(values = spec_cols, drop = TRUE, labels = c(
    autosomal_target_detectable_only = "Autosomal, detectable only in target sex core",
    SDR_target_detectable_only = "SDR, assumed sex-limited; opposite-sex TPM ignored",
    PAR_target_detectable_only = "PAR, detectable only in target sex core",
    other_target_detectable_only = "Other/unplaced, detectable only in target sex core",
    target_biased_but_detected_elsewhere = "Biased but detected in both sex cores",
    not_detected_at_floor = "Below TPM floor in target sex core"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Number of sex-biased genes",
    caption = "TPM floor is 1. Non-PAR ChrUf/ChrVm genes are treated as SDR-linked; opposite-sex TPM for the non-carried SDR is ignored for specificity. Genes below the target-sex TPM floor are excluded from this plot."
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_target_specificity_barplot.pdf"), p_spec, width = 7.2, height = 4.6, useDingbats = FALSE)

lifecycle_cols <- c(
  sex_biased_detected_in_both_sex_cores = "#D9D2C8",
  sex_core_target_only_but_lifecycle_detected = "#CC9966",
  sex_core_target_only_gametophyte_restricted = "#663333",
  not_detected_at_floor = "#9E9E9E"
)
p_lifecycle_spec <- sex_lifecycle_specificity_summary %>%
  mutate(
    set_label = factor(unname(spec_labels[set]), levels = unname(spec_labels[c("female_biased", "male_biased")])),
    lifecycle_specificity_class = factor(lifecycle_specificity_class, levels = names(lifecycle_cols)),
    class_label = recode(
      as.character(lifecycle_specificity_class),
      sex_biased_detected_in_both_sex_cores = "Detected in both sex cores",
      sex_core_target_only_but_lifecycle_detected = "Target-sex core only; also non-gametophyte",
      sex_core_target_only_gametophyte_restricted = "Target-sex core only; gametophyte-restricted",
      not_detected_at_floor = "Below target TPM floor"
    ),
    class_label = factor(class_label, levels = c(
      "Detected in both sex cores",
      "Target-sex core only; also non-gametophyte",
      "Target-sex core only; gametophyte-restricted",
      "Below target TPM floor"
    )),
    label = sprintf("%s", n_genes)
  ) %>%
  ggplot(aes(x = class_label, y = n_genes, fill = lifecycle_specificity_class)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 3.0, color = "grey15") +
  facet_wrap(~ set_label, nrow = 1) +
  scale_fill_manual(values = lifecycle_cols, drop = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Number of sex-biased genes",
    caption = "Gametophyte-restricted means detectable in the target sex core, not detectable in the opposing sex core, and TPM < 1 across cystocarp, carpospore, tetrasporophyte and tetraspore stages. Spermatia are treated as gametophyte, not non-gametophyte."
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "sex_biased_lifecycle_specificity_barplot.pdf"), p_lifecycle_spec, width = 7.6, height = 4.7, useDingbats = FALSE)

combined_specificity_region_plot_df <- sex_specificity %>%
  filter(target_detected, genomic_compartment %in% c("autosomal", "SDR_nonPAR", "PAR")) %>%
  mutate(
    set_label = factor(unname(spec_labels[set]), levels = unname(spec_labels[c("female_biased", "male_biased")])),
    compartment_label = recode(
      genomic_compartment,
      autosomal = "Autosome",
      SDR_nonPAR = "SDR/non-PAR",
      PAR = "PAR"
    ),
    compartment_label = factor(compartment_label, levels = c("Autosome", "SDR/non-PAR", "PAR")),
    combined_specificity_class = case_when(
      genomic_compartment == "SDR_nonPAR" & non_gametophyte_detected ~ "SDR-linked; detected elsewhere",
      genomic_compartment == "SDR_nonPAR" & !non_gametophyte_detected ~ "SDR-linked; gametophyte-restricted",
      lifecycle_specificity_class == "sex_biased_detected_in_both_sex_cores" ~ "Detected in both sex cores",
      lifecycle_specificity_class == "sex_core_target_only_but_lifecycle_detected" ~ "Target-sex only; detected elsewhere",
      lifecycle_specificity_class == "sex_core_target_only_gametophyte_restricted" ~ "Target-sex only; gametophyte-restricted",
      TRUE ~ "Other"
    ),
    combined_specificity_class = factor(
      combined_specificity_class,
      levels = c(
        "Other",
        "SDR-linked; detected elsewhere",
        "SDR-linked; gametophyte-restricted",
        "Detected in both sex cores",
        "Target-sex only; detected elsewhere",
        "Target-sex only; gametophyte-restricted"
      )
    )
  ) %>%
  count(set, set_label, compartment_label, combined_specificity_class, name = "n_genes") %>%
  group_by(set, compartment_label) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()
write_tsv(combined_specificity_region_plot_df, file.path(out_dir, "sex_biased_specificity_region_combined_summary.tsv"))

combined_specificity_cols <- c(
  "Detected in both sex cores" = "#D9D2C8",
  "Target-sex only; detected elsewhere" = "#CC9966",
  "Target-sex only; gametophyte-restricted" = "#663333",
  "SDR-linked; detected elsewhere" = "#8E2F3F",
  "SDR-linked; gametophyte-restricted" = "#5C1F2A",
  "Other" = "#9E9E9E"
)
combined_specificity_y_max <- max(combined_specificity_region_plot_df$total_genes, na.rm = TRUE)
combined_specificity_y_break_low <- 10
combined_specificity_y_break_high <- 25
combined_specificity_y_gap <- 3
combined_specificity_y_top <- 55
combined_specificity_y_scale <- (combined_specificity_y_top - combined_specificity_y_break_low - combined_specificity_y_gap) /
  max(1, combined_specificity_y_max - combined_specificity_y_break_high)
combined_specificity_y_transform <- function(y) {
  ifelse(
    y <= combined_specificity_y_break_low,
    y,
    combined_specificity_y_break_low + combined_specificity_y_gap +
      (y - combined_specificity_y_break_high) * combined_specificity_y_scale
  )
}
combined_specificity_y_breaks_raw <- sort(unique(c(
  0, 5, 10, 25,
  pretty(c(25, combined_specificity_y_max), n = 5)
)))
combined_specificity_y_breaks_raw <- combined_specificity_y_breaks_raw[
  combined_specificity_y_breaks_raw >= 0 & combined_specificity_y_breaks_raw <= combined_specificity_y_max
]
combined_specificity_rect_df <- combined_specificity_region_plot_df %>%
  arrange(set_label, compartment_label, combined_specificity_class) %>%
  group_by(set_label, compartment_label) %>%
  mutate(ymin_raw = lag(cumsum(n_genes), default = 0), ymax_raw = cumsum(n_genes)) %>%
  ungroup() %>%
  mutate(xnum = as.numeric(compartment_label), xmin = xnum - 0.36, xmax = xnum + 0.36)
combined_specificity_rect_low <- combined_specificity_rect_df %>%
  filter(ymin_raw < combined_specificity_y_break_low) %>%
  mutate(
    ymin_plot = combined_specificity_y_transform(ymin_raw),
    ymax_plot = combined_specificity_y_transform(pmin(ymax_raw, combined_specificity_y_break_low))
  ) %>%
  filter(ymax_plot > ymin_plot)
combined_specificity_rect_high <- combined_specificity_rect_df %>%
  filter(ymax_raw > combined_specificity_y_break_high) %>%
  mutate(
    ymin_plot = combined_specificity_y_transform(pmax(ymin_raw, combined_specificity_y_break_high)),
    ymax_plot = combined_specificity_y_transform(ymax_raw)
  ) %>%
  filter(ymax_plot > ymin_plot)
combined_specificity_rect_plot_df <- bind_rows(combined_specificity_rect_low, combined_specificity_rect_high)
combined_specificity_label_df <- combined_specificity_rect_plot_df %>%
  mutate(
    ymid_plot = (ymin_plot + ymax_plot) / 2,
    segment_height_plot = ymax_plot - ymin_plot,
    label = ifelse(n_genes > 0 & segment_height_plot >= 0.7, as.character(n_genes), ""),
    label_color = ifelse(combined_specificity_class %in% c(
      "Target-sex only; gametophyte-restricted",
      "SDR-linked; detected elsewhere",
      "SDR-linked; gametophyte-restricted"
    ), "white", "grey15")
  )
combined_specificity_break_df <- expand.grid(
  set_label = levels(combined_specificity_region_plot_df$set_label),
  y = c(
    combined_specificity_y_break_low + 0.85,
    combined_specificity_y_break_low + combined_specificity_y_gap - 0.85
  ),
  stringsAsFactors = FALSE
)
combined_specificity_break_df$set_label <- factor(
  combined_specificity_break_df$set_label,
  levels = levels(combined_specificity_region_plot_df$set_label)
)

p_combined_specificity_region <- ggplot(combined_specificity_rect_plot_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin_plot, ymax = ymax_plot, fill = combined_specificity_class),
    color = "white",
    linewidth = 0.2
  ) +
  geom_text(
    data = combined_specificity_label_df,
    aes(x = xnum, y = ymid_plot, label = label, color = label_color),
    size = 2.8,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  geom_hline(
    data = combined_specificity_break_df,
    aes(yintercept = y),
    color = "grey55",
    linewidth = 0.25,
    linetype = "22"
  ) +
  facet_wrap(~ set_label, nrow = 1) +
  scale_x_continuous(
    breaks = seq_along(levels(combined_specificity_region_plot_df$compartment_label)),
    labels = levels(combined_specificity_region_plot_df$compartment_label),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(
    breaks = combined_specificity_y_transform(combined_specificity_y_breaks_raw),
    labels = combined_specificity_y_breaks_raw,
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_fill_manual(values = combined_specificity_cols, drop = FALSE, name = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    legend.position = "top"
  ) +
  labs(
    x = NULL,
    y = "Target-detected sex-biased genes (axis split: 0-10, then 25+)",
    caption = "Facets split FBGs and MBGs; bars are genomic compartments. The y-axis hides 10-25 so low-count PAR classes remain visible while retaining the larger autosomal/SDR bars. Fill classes are stacked from broad/shared expression at the bottom to strongest target/lifecycle restriction at the top. Non-PAR ChrUf/ChrVm genes are treated as structurally SDR-linked, so opposite-sex TPM is ignored rather than interpreted as true expression."
  )
ggsave(file.path(fig_dir, "sex_biased_specificity_region_combined_barplot.pdf"), p_combined_specificity_region, width = 8.2, height = 4.8, useDingbats = FALSE)

# SDR-only expression heatmap across the lifecycle. This uses all non-PAR ChrUf/ChrVm
# protein-coding genes that are sex-biased and target-detected, with stage-mean TPMs.
sdr_expr_gene_df <- sex_specificity %>%
  filter(target_detected, genomic_compartment == "SDR_nonPAR") %>%
  left_join(coords[, c("gene_id", "start", "end", "strand")], by = "gene_id") %>%
  arrange(set, chr, start, end)
sdr_expr_genes <- sdr_expr_gene_df$gene_id
if (length(sdr_expr_genes) > 0) {
  sdr_hm_expr <- expr_means %>% filter(gene_id %in% sdr_expr_genes) %>% select(gene_id, all_of(mean_cols))
  sdr_hm_mat <- as.matrix(sdr_hm_expr[, mean_cols, drop = FALSE])
  storage.mode(sdr_hm_mat) <- "numeric"
  sdr_hm_mat <- log2(sdr_hm_mat + 1)
  rownames(sdr_hm_mat) <- sdr_hm_expr$gene_id
  sdr_hm_mat <- sdr_hm_mat[sdr_expr_genes[sdr_expr_genes %in% rownames(sdr_hm_mat)], , drop = FALSE]
  sdr_hm_z <- t(scale(t(sdr_hm_mat)))
  sdr_hm_z[!is.finite(sdr_hm_z)] <- 0
  colnames(sdr_hm_z) <- stage_ids

  sdr_row_df <- sdr_expr_gene_df %>%
    filter(gene_id %in% rownames(sdr_hm_z)) %>%
    distinct(gene_id, set, chr, region, start, end, lifecycle_specificity_class)
  sdr_row_split <- factor(
    paste0(ifelse(sdr_row_df$set[match(rownames(sdr_hm_z), sdr_row_df$gene_id)] == "female_biased", "FBG", "MBG"), " ", sdr_row_df$chr[match(rownames(sdr_hm_z), sdr_row_df$gene_id)]),
    levels = c("FBG ChrUf", "FBG ChrVm", "MBG ChrUf", "MBG ChrVm")
  )
  sdr_row_chr <- factor(sdr_row_df$chr[match(rownames(sdr_hm_z), sdr_row_df$gene_id)], levels = c("ChrUf", "ChrVm"))
  sdr_row_set <- factor(sdr_row_df$set[match(rownames(sdr_hm_z), sdr_row_df$gene_id)], levels = c("female_biased", "male_biased"), labels = c("Female-biased", "Male-biased"))

  sdr_heatmap_values <- data.frame(gene_id = rownames(sdr_hm_z), sdr_hm_z, check.names = FALSE)
  sdr_heatmap_values <- cbind(
    sdr_heatmap_values,
    sex_bias = as.character(sdr_row_set),
    chr = as.character(sdr_row_chr)
  )
  write_tsv(sdr_heatmap_values, file.path(out_dir, "sex_biased_sdr_expression_heatmap_values.tsv"))

  sdr_cell_mm <- 3.0
  sdr_pdf_width <- max(4.8, length(stage_ids) * sdr_cell_mm / 25.4 + 2.4)
  sdr_pdf_height <- max(5.0, nrow(sdr_hm_z) * sdr_cell_mm / 25.4 + 2.0)
  # Omitted non-manuscript PDF export: sex_biased_sdr_expression_heatmap.pdf
}


# All U/V chromosome gene expression heatmap across the lifecycle. This is not
# restricted to sex-biased genes: rows are all protein-coding genes on ChrUf/ChrVm,
# ordered by chromosome, PAR/SDR compartment, and genomic coordinate.
uv_gene_df <- region_df %>%
  filter(chr %in% c("ChrUf", "ChrVm")) %>%
  mutate(
    uv_compartment = ifelse(as.character(region) == "PAR", "PAR", "SDR"),
    uv_split = paste(chr, uv_compartment),
    chr = factor(chr, levels = c("ChrUf", "ChrVm")),
    uv_compartment = factor(uv_compartment, levels = c("PAR", "SDR")),
    uv_split = factor(uv_split, levels = c("ChrUf PAR", "ChrUf SDR", "ChrVm PAR", "ChrVm SDR"))
  ) %>%
  arrange(chr, uv_compartment, start, end, gene_id)
uv_genes <- uv_gene_df$gene_id
if (length(uv_genes) > 0) {
  uv_hm_expr <- expr_means %>% filter(gene_id %in% uv_genes) %>% select(gene_id, all_of(mean_cols))
  uv_hm_mat <- as.matrix(uv_hm_expr[, mean_cols, drop = FALSE])
  storage.mode(uv_hm_mat) <- "numeric"
  uv_hm_mat <- log2(uv_hm_mat + 1)
  rownames(uv_hm_mat) <- uv_hm_expr$gene_id
  uv_hm_mat <- uv_hm_mat[uv_genes[uv_genes %in% rownames(uv_hm_mat)], , drop = FALSE]
  uv_hm_z <- t(scale(t(uv_hm_mat)))
  uv_hm_z[!is.finite(uv_hm_z)] <- 0
  colnames(uv_hm_z) <- stage_ids

  uv_row_df <- uv_gene_df %>%
    filter(gene_id %in% rownames(uv_hm_z)) %>%
    distinct(gene_id, chr, uv_compartment, uv_split, start, end) %>%
    arrange(chr, uv_compartment, start, end, gene_id)
  uv_row_df <- uv_row_df[match(rownames(uv_hm_z), uv_row_df$gene_id), , drop = FALSE]
  uv_row_chr <- factor(as.character(uv_row_df$chr), levels = c("ChrUf", "ChrVm"))
  uv_row_compartment <- factor(as.character(uv_row_df$uv_compartment), levels = c("PAR", "SDR"))
  uv_row_split <- factor(as.character(uv_row_df$uv_split), levels = c("ChrUf PAR", "ChrUf SDR", "ChrVm PAR", "ChrVm SDR"))

  uv_heatmap_values <- data.frame(
    gene_id = rownames(uv_hm_z),
    chr = as.character(uv_row_chr),
    compartment = as.character(uv_row_compartment),
    start = uv_row_df$start,
    end = uv_row_df$end,
    uv_hm_z,
    check.names = FALSE
  )
  write_tsv(uv_heatmap_values, file.path(out_dir, "uv_all_genes_expression_heatmap_values.tsv"))

  uv_hm_z_masked <- uv_hm_z
  male_haploid_stages <- intersect(c("Bm3235m_vegtissue", "Bm3235m_spermatangia", "Bm3235m_spermatia"), colnames(uv_hm_z_masked))
  female_haploid_stages <- intersect(c("Bm3235f_vegtissue", "Bm3235f_procarps", "Bm3235f_pericarp"), colnames(uv_hm_z_masked))
  u_sdr_rows <- as.character(uv_row_chr) == "ChrUf" & as.character(uv_row_compartment) == "SDR"
  v_sdr_rows <- as.character(uv_row_chr) == "ChrVm" & as.character(uv_row_compartment) == "SDR"
  if (length(male_haploid_stages) > 0) {
    uv_hm_z_masked[u_sdr_rows, male_haploid_stages] <- NA_real_
  }
  if (length(female_haploid_stages) > 0) {
    uv_hm_z_masked[v_sdr_rows, female_haploid_stages] <- NA_real_
  }
  uv_heatmap_masked_values <- data.frame(
    gene_id = rownames(uv_hm_z_masked),
    chr = as.character(uv_row_chr),
    compartment = as.character(uv_row_compartment),
    start = uv_row_df$start,
    end = uv_row_df$end,
    uv_hm_z_masked,
    check.names = FALSE
  )
  write_tsv(uv_heatmap_masked_values, file.path(out_dir, "uv_all_genes_expression_heatmap_presence_masked_values.tsv"))

  uv_cell_w_mm <- 3.0
  uv_cell_h_mm <- 0.55
  uv_pdf_width <- max(5.4, length(stage_ids) * uv_cell_w_mm / 25.4 + 2.8)
  uv_pdf_height <- max(7.0, min(14.0, nrow(uv_hm_z) * uv_cell_h_mm / 25.4 + 2.6))
  # Omitted non-manuscript PDF export: uv_all_genes_expression_heatmap.pdf

  # Omitted non-manuscript PDF export: uv_all_genes_expression_heatmap_presence_masked.pdf

  # Omitted non-manuscript PDF export: uv_all_genes_expression_heatmap_clustered.pdf

  pdf(file.path(fig_dir, "uv_all_genes_expression_heatmap_presence_masked_clustered.pdf"), width = uv_pdf_width, height = uv_pdf_height)
  ht_uv_masked_clustered <- Heatmap(
    uv_hm_z_masked,
    name = "row z",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#6B6B6B", "white", "#663333")),
    na_col = "#D9D9D9",
    top_annotation = HeatmapAnnotation(
      stage = stage_ids,
      col = list(stage = stage_cols_present),
      show_annotation_name = FALSE,
      annotation_legend_param = list(stage = list(ncol = 1))
    ),
    left_annotation = rowAnnotation(
      chr = uv_row_chr,
      region = uv_row_compartment,
      col = list(
        chr = c(ChrUf = "#8E2F3F", ChrVm = "#5C1F2A"),
        region = c(PAR = "#6A3D9A", SDR = "#C8B78A")
      ),
      show_annotation_name = TRUE
    ),
    row_split = uv_row_split,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_labels = stage_ids,
    column_names_rot = 45,
    width = grid::unit(ncol(uv_hm_z_masked) * uv_cell_w_mm, "mm"),
    height = grid::unit(nrow(uv_hm_z_masked) * uv_cell_h_mm, "mm"),
    row_gap = grid::unit(1.2, "mm"),
    rect_gp = grid::gpar(col = NA),
    use_raster = TRUE,
    heatmap_legend_param = list(title = "row z")
  )
  draw(ht_uv_masked_clustered, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# --------- Integrated comparison: triphasy identity sets + sex-biased sets ---------
trip_set_levels <- c(
  "GAM_core_biased", "GAM_core_downregulated",
  "CSP_biased", "CSP_downregulated",
  "TETRA_biased", "TETRA_downregulated"
)
combined_set_levels <- c(trip_set_levels, "female_biased", "male_biased", "genome_background")
combined_set_labels <- c(
  GAM_core_biased = "GAM biased",
  GAM_core_downregulated = "GAM repressed",
  CSP_biased = "Carposporophyte biased",
  CSP_downregulated = "Carposporophyte repressed",
  TETRA_biased = "Tetrasporophyte biased",
  TETRA_downregulated = "Tetrasporophyte repressed",
  female_biased = "Female-biased genes",
  male_biased = "Male-biased genes",
  genome_background = "All protein-coding genes"
)
combined_set_cols <- c(
  GAM_core_biased = "#8E2F3F",
  GAM_core_downregulated = "#C99AA3",
  CSP_biased = "#6A3D9A",
  CSP_downregulated = "#B79AC8",
  TETRA_biased = "#8FA33E",
  TETRA_downregulated = "#C8D38A",
  female_biased = "#5C1F2A",
  male_biased = "#B8860B",
  genome_background = "#8A8A8A"
)

trip_membership_files <- c(
  file.path(paths$output_root, "generation", "triphasy_biased_gene_set_membership.tsv"),
  file.path(paths$output_root, "generation", "triphasy_downregulated_gene_set_membership.tsv")
)
trip_membership <- bind_rows(lapply(trip_membership_files[file.exists(trip_membership_files)], function(f) {
  read.delim(f, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)[, c("gene_id", "set"), drop = FALSE]
})) %>%
  mutate(gene_id = normalize_gene_id(gene_id)) %>%
  filter(set %in% trip_set_levels, gene_id %in% pc_genes) %>%
  distinct(gene_id, set)
sex_membership <- sex_sets_target_detected %>%
  transmute(gene_id, set) %>%
  filter(gene_id %in% pc_genes) %>%
  distinct(gene_id, set)
combined_membership <- bind_rows(trip_membership, sex_membership) %>%
  filter(set %in% combined_set_levels) %>%
  distinct(gene_id, set)
write_tsv(combined_membership, file.path(out_dir, "triphasy_sex_bias_gene_set_membership.tsv"))

# TAP: use all protein-coding genes as the background; missing TAP annotation is treated as non-TAP.
pc_universe <- pc_genes
anno_pc <- read_functional_annotation(paths) %>%
  transmute(gene_id = normalize_gene_id(gene_id), TAP_family = as.character(TAP_family)) %>%
  filter(gene_id %in% pc_universe) %>%
  distinct(gene_id, .keep_all = TRUE)
anno_pc <- data.frame(gene_id = pc_universe, stringsAsFactors = FALSE) %>%
  left_join(anno_pc, by = "gene_id") %>%
  mutate(
    TAP_family = ifelse(is.na(TAP_family) | !nzchar(trimws(TAP_family)) | trimws(TAP_family) %in% c("-", "NA"), "non_TAP", TAP_family),
    tap_status = ifelse(TAP_family == "non_TAP", "non-TAP", "TAP")
  )
combined_tap_summary <- bind_rows(
  lapply(combined_set_levels[combined_set_levels != "genome_background"], function(s) {
    gids <- unique(combined_membership$gene_id[combined_membership$set == s])
    anno_pc %>% filter(gene_id %in% gids) %>% count(tap_status, name = "n_genes") %>% mutate(set = s)
  }),
  anno_pc %>% count(tap_status, name = "n_genes") %>% mutate(set = "genome_background")
) %>%
  complete(set = combined_set_levels, tap_status = c("TAP", "non-TAP"), fill = list(n_genes = 0L)) %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), proportion = n_genes / total_genes) %>%
  ungroup() %>%
  mutate(
    set = factor(set, levels = combined_set_levels),
    set_label = factor(unname(combined_set_labels[as.character(set)]), levels = unname(combined_set_labels[combined_set_levels]))
  )
write_tsv(combined_tap_summary, file.path(out_dir, "triphasy_sex_bias_tap_proportion_summary.tsv"))

combined_tap_tests <- bind_rows(lapply(combined_set_levels[combined_set_levels != "genome_background"], function(s) {
  gids <- unique(combined_membership$gene_id[combined_membership$set == s])
  in_set <- pc_universe %in% gids
  is_tap <- anno_pc$tap_status[match(pc_universe, anno_pc$gene_id)] == "TAP"
  mat <- matrix(c(sum(in_set & is_tap), sum(in_set & !is_tap), sum(!in_set & is_tap), sum(!in_set & !is_tap)), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "two.sided")
  data.frame(
    set = s,
    set_label = unname(combined_set_labels[s]),
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
})) %>%
  mutate(fdr = p.adjust(fisher_p, method = "BH")) %>%
  arrange(fdr, fisher_p)
write_tsv(combined_tap_tests, file.path(out_dir, "triphasy_sex_bias_tap_enrichment_fisher.tsv"))

combined_tap_plot_df <- combined_tap_summary %>%
  filter(tap_status == "TAP") %>%
  left_join(combined_tap_tests[, c("set", "fdr", "log2_odds_ratio")], by = "set") %>%
  mutate(
    label = sprintf("%.1f%%\n%s/%s", 100 * proportion, n_genes, total_genes),
    sig_label = case_when(
      is.na(fdr) ~ "",
      fdr <= 0.001 ~ "***",
      fdr <= 0.01 ~ "**",
      fdr <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
tap_ymax <- max(0.12, max(combined_tap_plot_df$proportion, na.rm = TRUE) * 1.25)

p_combined_tap <- ggplot(combined_tap_plot_df, aes(x = set_label, y = proportion, fill = set)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = label), vjust = -0.25, size = 2.65, lineheight = 0.88) +
  geom_text(aes(y = proportion + max(proportion, na.rm = TRUE) * 0.09, label = sig_label), size = 3.0, color = "grey10", na.rm = TRUE) +
  scale_fill_manual(values = combined_set_cols, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) sprintf("%d%%", round(100*x)),
    breaks = seq(0, 0.20, by = 0.025),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.06))
  ) +
  coord_cartesian(ylim = c(0, tap_ymax), clip = "off") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), legend.position = "none") +
  labs(x = NULL, y = "TAP genes (%)", caption = "Labels show TAP percentage and n TAP / total genes; * FDR <= 0.05, ** FDR <= 0.01, *** FDR <= 0.001 (two-sided Fisher test against protein-coding background excluding the set).")
ggsave(file.path(fig_dir, "triphasy_sex_bias_tap_proportion_barplot.pdf"), p_combined_tap, width = 9.4, height = 4.8, useDingbats = FALSE)

tap_heat <- combined_tap_tests %>%
  mutate(
    set_label = factor(set_label, levels = unname(combined_set_labels[combined_set_levels[combined_set_levels != "genome_background"]])),
    metric = "TAP enrichment",
    log2_or_label = sprintf("%.2f", log2_odds_ratio),
    fdr_sig = is.finite(fdr) & fdr <= 0.05
  )
tap_vals <- tap_heat$log2_odds_ratio[is.finite(tap_heat$log2_odds_ratio)]
tap_min <- min(tap_vals, na.rm = TRUE)
tap_max <- max(tap_vals, na.rm = TRUE)
if (!is.finite(tap_min) || !is.finite(tap_max) || tap_min == tap_max) {
  tap_min <- -1
  tap_max <- 1
}
tap_mid_pos <- if (tap_max > 0) tap_max * 0.45 else tap_max
p_tap_heat <- ggplot(tap_heat, aes(x = set_label, y = metric, fill = log2_odds_ratio)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = log2_or_label, fontface = if_else(fdr_sig, "bold", "plain")), size = 2.55, color = "black") +
  scale_fill_gradientn(
    colors = c("#6B6B6B", "white", "#CC9966", "#663333"),
    values = scales::rescale(c(tap_min, 0, tap_mid_pos, tap_max), from = c(tap_min, tap_max)),
    limits = c(tap_min, tap_max),
    oob = scales::squish,
    na.value = "#D9D9D9",
    name = "log2 odds\nratio"
  ) +
  labs(
    title = "TAP enrichment across triphasy and sex-biased gene sets",
    subtitle = "Positive = TAP enriched; negative = TAP depleted; bold labels indicate Fisher FDR <= 0.05",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_sex_bias_tap_enrichment_log2_oddsratio_heatmap.pdf"), p_tap_heat, width = 9.2, height = 2.5, useDingbats = FALSE)

# Phylorank composition and enrichment for same integrated sets.
age_pc <- read_gene_age_table(paths$gene_age_tsv) %>%
  mutate(
    gene_id = normalize_gene_id(gene_id),
    rank = suppressWarnings(as.integer(rank)),
    rank_label = ifelse(!is.na(rank), paste0("PS", rank, ": ", taxon), as.character(taxon))
  ) %>%
  filter(gene_id %in% pc_universe, !is.na(rank), !is.na(rank_label), nzchar(rank_label)) %>%
  distinct(gene_id, .keep_all = TRUE)
rank_levels_combined <- age_pc %>% distinct(rank, rank_label) %>% arrange(rank) %>% pull(rank_label)
pal_rank_combined <- setNames(rep("#999999", length(rank_levels_combined)), rank_levels_combined)
pal_rank_combined[grep("^PS1:", rank_levels_combined)] <- "#f7efe5"
pal_rank_combined[grep("^PS2:", rank_levels_combined)] <- "#e7d5bc"
pal_rank_combined[grep("^PS3:", rank_levels_combined)] <- "#8CC8C1"
pal_rank_combined[grep("^PS4:", rank_levels_combined)] <- "#4FA79D"
pal_rank_combined[grep("^PS5:", rank_levels_combined)] <- "#2F7F78"
pal_rank_combined[grep("^PS6:", rank_levels_combined)] <- "#B35E5E"
pal_rank_combined[grep("^PS7:", rank_levels_combined)] <- "#8F3A3A"
pal_rank_combined[grep("^PS8:", rank_levels_combined)] <- "#662424"

combined_rank_dist <- bind_rows(
  combined_membership %>% inner_join(age_pc[, c("gene_id", "rank_label")], by = "gene_id") %>% count(set, rank_label, name = "value"),
  age_pc %>% count(rank_label, name = "value") %>% mutate(set = "genome_background")
) %>%
  complete(set = combined_set_levels, rank_label = rank_levels_combined, fill = list(value = 0L)) %>%
  mutate(
    set = factor(set, levels = combined_set_levels),
    set_label = factor(unname(combined_set_labels[as.character(set)]), levels = unname(combined_set_labels[combined_set_levels])),
    rank_label = factor(rank_label, levels = rank_levels_combined)
  )
write_tsv(combined_rank_dist, file.path(out_dir, "triphasy_sex_bias_phylorank_distribution.tsv"))

age_universe_combined <- age_pc$gene_id
rank_test_grid <- expand.grid(set = combined_set_levels[combined_set_levels != "genome_background"], rank_label = rank_levels_combined, stringsAsFactors = FALSE)
combined_rank_tests <- bind_rows(lapply(seq_len(nrow(rank_test_grid)), function(i) {
  s <- rank_test_grid$set[i]
  rr <- rank_test_grid$rank_label[i]
  gids <- unique(combined_membership$gene_id[combined_membership$set == s])
  in_set <- age_universe_combined %in% gids
  in_rank <- age_pc$rank_label[match(age_universe_combined, age_pc$gene_id)] == rr
  mat <- matrix(c(sum(in_set & in_rank), sum(in_set & !in_rank), sum(!in_set & in_rank), sum(!in_set & !in_rank)), nrow = 2, byrow = TRUE)
  ft_enrich <- fisher.test(mat, alternative = "greater")
  ft_deplete <- fisher.test(mat, alternative = "less")
  data.frame(
    set = s,
    set_label = unname(combined_set_labels[s]),
    rank_label = rr,
    observed = sum(in_set & in_rank),
    set_total = sum(in_set),
    background_total_rank = sum(in_rank),
    background_total = length(age_universe_combined),
    set_proportion = sum(in_set & in_rank) / max(1, sum(in_set)),
    background_proportion = sum(in_rank) / max(1, length(age_universe_combined)),
    log2_obs_exp = log2((sum(in_set & in_rank) + 0.5) / ((sum(in_set) * sum(in_rank) / length(age_universe_combined)) + 0.5)),
    odds_ratio_enrich = unname(ft_enrich$estimate),
    fisher_enrich_p = ft_enrich$p.value,
    fisher_deplete_p = ft_deplete$p.value,
    stringsAsFactors = FALSE
  )
})) %>%
  group_by(set) %>%
  mutate(
    fisher_enrich_fdr_within_set = p.adjust(fisher_enrich_p, method = "BH"),
    fisher_deplete_fdr_within_set = p.adjust(fisher_deplete_p, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    fisher_enrich_fdr_global = p.adjust(fisher_enrich_p, method = "BH"),
    fisher_deplete_fdr_global = p.adjust(fisher_deplete_p, method = "BH"),
    significant_enriched_within_set = fisher_enrich_fdr_within_set <= 0.05,
    significant_depleted_within_set = fisher_deplete_fdr_within_set <= 0.05
  ) %>%
  arrange(factor(set, levels = combined_set_levels), factor(rank_label, levels = rank_levels_combined))
write_tsv(combined_rank_tests, file.path(out_dir, "triphasy_sex_bias_phylorank_enrichment.tsv"))

p_combined_rank <- ggplot(combined_rank_dist, aes(x = set_label, y = value, fill = rank_label)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.72, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = pal_rank_combined, breaks = rank_levels_combined, limits = rank_levels_combined, drop = FALSE) +
  scale_y_continuous(labels = function(x) sprintf("%d%%", round(100*x))) +
  coord_flip() +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right", legend.title = element_blank()) +
  labs(x = NULL, y = "Relative composition")
ggsave(file.path(fig_dir, "triphasy_sex_bias_phylorank_stacked_bar_relative.pdf"), p_combined_rank, width = 10.2, height = 7.6, useDingbats = FALSE)

p_combined_rank_heat <- combined_rank_tests %>%
  mutate(
    set_label = factor(set_label, levels = rev(unname(combined_set_labels[combined_set_levels[combined_set_levels != "genome_background"]]))),
    rank_label = factor(rank_label, levels = rank_levels_combined),
    sig_mark = ifelse(fisher_enrich_fdr_within_set <= 0.05, "*", ifelse(fisher_deplete_fdr_within_set <= 0.05, "-", ""))
  ) %>%
  ggplot(aes(x = rank_label, y = set_label, fill = log2_obs_exp)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig_mark), size = 3.4, color = "black") +
  scale_fill_gradient2(low = "#bdbdbd", mid = "white", high = "#663333", midpoint = 0, name = "log2(obs/exp)") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylorank", y = NULL, subtitle = "* enriched, - depleted; within-set Fisher BH FDR <= 0.05")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_sex_bias_phylorank_enrichment_heatmap.pdf"), p_combined_rank_heat, width = 9.8, height = 5.8, useDingbats = FALSE)

# --------- Audit ---------
aud <- data.frame(
  metric = c(
    "female_core_samples", "male_core_samples", "excluded_spermatia_samples",
    "protein_coding_universe_with_coordinates", "tested_protein_coding_genes",
    "female_biased_genes", "male_biased_genes",
    "chruf_genes", "chrvm_genes", "autosome_genes",
    "par_file", "par_present", "par_regions", "par_genes",
    "organellar_chromosomes_excluded", "organellar_protein_coding_genes_excluded"
  ),
  value = c(
    paste(core_design$sample_id[core_design$include & core_design$sex_core_group == "female_core"], collapse = ","),
    paste(core_design$sample_id[core_design$include & core_design$sex_core_group == "male_core"], collapse = ","),
    paste(design$sample_id[design$stage %in% excluded_stages & design$include], collapse = ","),
    length(universe), sum(de$gene_id %in% universe & is.finite(de$fdr)),
    sum(de_region$sex_bias == "female_biased"), sum(de_region$sex_bias == "male_biased"),
    sum(region_df$region == "ChrUf"), sum(region_df$region == "ChrVm"), sum(region_df$region == "autosome"),
    ifelse(is.na(par_file) || !nzchar(par_file), "not_found", par_file), par_present, par_n_regions, par_n_genes,
    paste(organellar_chrs, collapse = ","), length(organellar_genes_excluded)
  ),
  stringsAsFactors = FALSE
)
write_tsv(aud, file.path(out_dir, "sex_biased_audit.tsv"))

cat("Wrote sex-bias analysis outputs to:\n")
cat(" - ", out_dir, "\n", sep = "")
cat(" - ", fig_dir, "\n", sep = "")
