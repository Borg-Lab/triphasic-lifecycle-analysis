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
source("R/tai.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))

out_dir <- file.path(paths$output_root, "generation")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

membership_file <- file.path(out_dir, "triphasy_biased_gene_set_membership.tsv")
if (!file.exists(membership_file)) {
  stop("Missing ", membership_file, ". Run 73_generation_transcriptome_overlap.R first.", call. = FALSE)
}

set_levels <- c("GAM_core_biased", "CSP_biased", "TETRA_biased")
set_labels <- c(
  GAM_core_biased = "GAM biased",
  TETRA_biased = "Tetrasporophyte biased",
  CSP_biased = "Carposporophyte biased"
)
set_cols <- c(
  GAM_core_biased = "#8E2F3F",
  TETRA_biased = "#8FA33E",
  CSP_biased = "#6A3D9A"
)

down_membership_file <- file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv")
down_set_levels <- c("GAM_core_downregulated", "CSP_downregulated", "TETRA_downregulated")
down_set_labels <- c(
  GAM_core_downregulated = "GAM downregulated",
  TETRA_downregulated = "Tetrasporophyte downregulated",
  CSP_downregulated = "Carposporophyte downregulated"
)
down_set_cols <- c(
  GAM_core_downregulated = "#8E2F3F",
  TETRA_downregulated = "#8FA33E",
  CSP_downregulated = "#6A3D9A"
)

make_tap_row_annotation <- function(tap_family, max_families = 10) {
  tap_family <- as.character(tap_family)
  tap_family[is.na(tap_family) | !nzchar(trimws(tap_family)) | trimws(tap_family) %in% c("-", "NA")] <- "non_TAP"
  tap_status <- factor(ifelse(tap_family == "non_TAP", "non-TAP", "TAP"), levels = c("non-TAP", "TAP"))

  fam_counts <- sort(table(tap_family[tap_family != "non_TAP"]), decreasing = TRUE)
  top_fams <- names(head(fam_counts, max_families))
  fam_limited <- ifelse(tap_family == "non_TAP", "non_TAP", ifelse(tap_family %in% top_fams, tap_family, "other_TAP"))
  fam_levels <- unique(c("non_TAP", top_fams, if (any(fam_limited == "other_TAP")) "other_TAP" else character()))
  fam_limited <- factor(fam_limited, levels = fam_levels)

  fam_cols <- setNames(rep("#F0F0F0", length(fam_levels)), fam_levels)
  tap_levels <- setdiff(fam_levels, "non_TAP")
  if (length(tap_levels) > 0) {
    fam_cols[tap_levels] <- grDevices::hcl.colors(length(tap_levels), palette = "Dark 3")
  }
  if ("other_TAP" %in% names(fam_cols)) fam_cols[["other_TAP"]] <- "#8C8C8C"

  rowAnnotation(
    TAP = tap_status,
    TAP_family = fam_limited,
    col = list(
      TAP = c("non-TAP" = "#F0F0F0", TAP = "#2B2B2B"),
      TAP_family = fam_cols
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 7),
    annotation_legend_param = list(
      TAP = list(title = "TAP"),
      TAP_family = list(title = "TAP family", ncol = 1)
    ),
    width = unit(8, "mm")
  )
}

stage_cols <- c(
  "Bm3235f_vegtissue"      = "#8E2F3F",
  "Bm3235f_procarps"       = "#5C1F2A",
  "Bm3235m_vegtissue"      = "#D6AF3D",
  "Bm3235m_spermatangia"   = "#B8860B",
  "Bm3235m_spermatia"      = "#8C6A00",
  "Bm3235f_cystocarps_imm" = "#B79AC8",
  "Bm3235f_cystocarps_mat" = "#6E4B9E",
  "Bm3235t_carpospores"    = "#4C2A85",
  "Bm3235t_vegtissue"      = "#B5C46A",
  "Bm3235t_stichidia"      = "#8FA33E",
  "Bm3235t_tetraspores"    = "#6E8629"
)

expr <- prepare_expression_inputs(paths, design, analysis)
means <- compute_group_tpm_means(expr$tpm, design, group_col = "stage")
means$gene_id <- normalize_gene_id(means$gene_id)

membership <- read.delim(membership_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
membership$gene_id <- normalize_gene_id(membership$gene_id)
membership <- membership %>%
  filter(set %in% set_levels) %>%
  mutate(
    set = factor(set, levels = set_levels),
    set_label = factor(unname(set_labels[as.character(set)]), levels = unname(set_labels[set_levels]))
  )

tai_order_file <- file.path(paths$output_root, "tai", "tai_plot_input.tsv")
mean_cols_all <- grep("^mean_tpm_", colnames(means), value = TRUE)
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
missing_cols <- setdiff(stage_ids, names(stage_cols_present))
if (length(missing_cols) > 0) {
  fallback <- rep("#999999", length(missing_cols))
  names(fallback) <- missing_cols
  stage_cols_present <- c(stage_cols_present, fallback)
}
stage_cols_present <- stage_cols_present[stage_ids]

expr_use <- means %>%
  inner_join(membership, by = "gene_id")

mat_log <- as.matrix(expr_use[, mean_cols, drop = FALSE])
storage.mode(mat_log) <- "numeric"
mat_log <- log2(mat_log + 1)
rownames(mat_log) <- expr_use$gene_id

row_z <- t(scale(t(mat_log)))
row_z[!is.finite(row_z)] <- 0
colnames(row_z) <- stage_ids

row_split <- expr_use$set_label
names(row_split) <- expr_use$gene_id
row_split <- factor(row_split, levels = unname(set_labels[set_levels]))

tap_anno <- data.frame(gene_id = rownames(row_z), TAP_family = "non_TAP", stringsAsFactors = FALSE)
marker_table_path <- file.path(paths$output_root, "markers", "marker_table.tsv")
if (file.exists(marker_table_path)) {
  marker_tab <- read.delim(marker_table_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("gene_id", "TAP_family") %in% colnames(marker_tab))) {
    marker_tab$gene_id <- normalize_gene_id(marker_tab$gene_id)
    marker_tab <- marker_tab[!duplicated(marker_tab$gene_id), c("gene_id", "TAP_family"), drop = FALSE]
    fam <- marker_tab$TAP_family[match(tap_anno$gene_id, marker_tab$gene_id)]
    fam <- as.character(fam)
    fam[is.na(fam) | !nzchar(trimws(fam)) | trimws(fam) %in% c("-", "NA")] <- "non_TAP"
    tap_anno$TAP_family <- fam
  }
} else {
  func_anno <- read_functional_annotation(paths)
  if (all(c("gene_id", "TAP_family") %in% colnames(func_anno))) {
    func_anno <- func_anno[!duplicated(func_anno$gene_id), c("gene_id", "TAP_family"), drop = FALSE]
    fam <- func_anno$TAP_family[match(tap_anno$gene_id, func_anno$gene_id)]
    fam <- as.character(fam)
    fam[is.na(fam) | !nzchar(trimws(fam)) | trimws(fam) %in% c("-", "NA")] <- "non_TAP"
    tap_anno$TAP_family <- fam
  }
}
tap_anno$set <- as.character(expr_use$set)
tap_anno$set_label <- as.character(expr_use$set_label)
tap_hits <- tap_anno[tap_anno$TAP_family != "non_TAP", , drop = FALSE]
tap_summary <- tap_hits %>%
  count(set, set_label, TAP_family, name = "n_genes") %>%
  arrange(set, desc(n_genes), TAP_family)
write_tsv(tap_hits, file.path(out_dir, "triphasy_biased_tap_genes.tsv"))
write_tsv(tap_summary, file.path(out_dir, "triphasy_biased_tap_family_summary.tsv"))

heatmap_values <- data.frame(gene_id = rownames(row_z), row_z, check.names = FALSE)
heatmap_values <- cbind(
  heatmap_values,
  set = as.character(expr_use$set),
  set_label = as.character(expr_use$set_label),
  TAP_family = tap_anno$TAP_family
)
write_tsv(heatmap_values, file.path(out_dir, "triphasy_biased_expression_heatmap_values.tsv"))

profile_long <- expr_use %>%
  select(gene_id, set, set_label, all_of(mean_cols)) %>%
  pivot_longer(cols = all_of(mean_cols), names_to = "stage_col", values_to = "mean_tpm") %>%
  mutate(
    stage = sub("^mean_tpm_", "", stage_col),
    stage = factor(stage, levels = stage_ids),
    stage_order = as.integer(stage),
    log2_tpm = log2(mean_tpm + 1)
  )

profile_summary <- profile_long %>%
  group_by(set, set_label, stage, stage_order) %>%
  summarise(
    n_genes = n(),
    mean_log2_tpm = mean(log2_tpm, na.rm = TRUE),
    se_log2_tpm = sd(log2_tpm, na.rm = TRUE) / sqrt(n()),
    median_log2_tpm = median(log2_tpm, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(profile_long, file.path(out_dir, "triphasy_biased_expression_profile_gene_values.tsv"))
write_tsv(profile_summary, file.path(out_dir, "triphasy_biased_expression_profile_summary.tsv"))

ha_right <- make_tap_row_annotation(tap_anno$TAP_family)

pdf(file.path(fig_dir, "triphasy_biased_expression_heatmap.pdf"), width = 8.5, height = 10.5)
ht <- Heatmap(
  row_z,
  name = "row z",
  col = circlize::colorRamp2(c(-2, 0, 2), c("#cccccc", "white", "#663333")),
  top_annotation = HeatmapAnnotation(
    stage = stage_ids,
    col = list(stage = stage_cols_present),
    show_annotation_name = FALSE,
    annotation_legend_param = list(stage = list(ncol = 1))
  ),
  left_annotation = rowAnnotation(
    category = row_split,
    col = list(category = setNames(unname(set_cols[set_levels]), unname(set_labels[set_levels]))),
    show_annotation_name = FALSE
  ),
  right_annotation = ha_right,
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

p_profile <- ggplot(profile_summary, aes(x = stage, y = mean_log2_tpm, group = set_label)) +
  geom_ribbon(
    aes(ymin = mean_log2_tpm - se_log2_tpm, ymax = mean_log2_tpm + se_log2_tpm),
    fill = "grey80",
    alpha = 0.55,
    color = NA
  ) +
  geom_line(color = "black", linewidth = 1.0) +
  geom_point(aes(color = stage), size = 3.0) +
  facet_wrap(~ set_label, ncol = 1, scales = "free_y") +
  scale_color_manual(values = stage_cols_present, drop = FALSE) +
  labs(
    x = NULL,
    y = "Mean log2(TPM + 1)",
    title = "Expression Profiles of Triphasy-Biased Gene Sets",
    subtitle = "Mean across genes with +/- SE ribbon; points use the TAI/PCA stage palette"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_expression_line_profiles.pdf"), p_profile, width = 8.4, height = 8.8)
# Omitted non-manuscript PDF export: triphasy_biased_expression_line_profiles_with_legend.pdf

p_profile_overlay <- ggplot(profile_summary, aes(x = stage, y = mean_log2_tpm, group = set_label)) +
  geom_ribbon(
    aes(ymin = mean_log2_tpm - se_log2_tpm, ymax = mean_log2_tpm + se_log2_tpm, fill = set_label),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(aes(color = set_label), linewidth = 1.15) +
  geom_point(aes(color = set_label), size = 2.8) +
  scale_color_manual(values = setNames(unname(set_cols[set_levels]), unname(set_labels[set_levels])), name = NULL) +
  scale_fill_manual(values = setNames(unname(set_cols[set_levels]), unname(set_labels[set_levels])), name = NULL) +
  labs(
    x = NULL,
    y = "Mean log2(TPM + 1)",
    title = "Expression Profiles of Triphasy-Biased Gene Sets",
    subtitle = "Mean across genes with +/- SE ribbon"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_expression_line_profiles_overlay.pdf"), p_profile_overlay, width = 8.8, height = 4.8)

if (file.exists(down_membership_file)) {
  down_membership <- read.delim(down_membership_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  down_membership$gene_id <- normalize_gene_id(down_membership$gene_id)
  down_membership <- down_membership %>%
    filter(set %in% down_set_levels) %>%
    mutate(
      set = factor(set, levels = down_set_levels),
      set_label = factor(unname(down_set_labels[as.character(set)]), levels = unname(down_set_labels[down_set_levels]))
    )

  down_expr_use <- means %>%
    inner_join(down_membership, by = "gene_id")

  down_mat_log <- as.matrix(down_expr_use[, mean_cols, drop = FALSE])
  storage.mode(down_mat_log) <- "numeric"
  down_mat_log <- log2(down_mat_log + 1)
  rownames(down_mat_log) <- down_expr_use$gene_id

  down_row_z <- t(scale(t(down_mat_log)))
  down_row_z[!is.finite(down_row_z)] <- 0
  colnames(down_row_z) <- stage_ids

  down_row_split <- down_expr_use$set_label
  names(down_row_split) <- down_expr_use$gene_id
  down_row_split <- factor(down_row_split, levels = unname(down_set_labels[down_set_levels]))

  down_tap_anno <- data.frame(gene_id = rownames(down_row_z), TAP_family = "non_TAP", stringsAsFactors = FALSE)
  if (file.exists(marker_table_path)) {
    marker_tab <- read.delim(marker_table_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    if (all(c("gene_id", "TAP_family") %in% colnames(marker_tab))) {
      marker_tab$gene_id <- normalize_gene_id(marker_tab$gene_id)
      marker_tab <- marker_tab[!duplicated(marker_tab$gene_id), c("gene_id", "TAP_family"), drop = FALSE]
      fam <- marker_tab$TAP_family[match(down_tap_anno$gene_id, marker_tab$gene_id)]
      fam <- as.character(fam)
      fam[is.na(fam) | !nzchar(trimws(fam)) | trimws(fam) %in% c("-", "NA")] <- "non_TAP"
      down_tap_anno$TAP_family <- fam
    }
  } else {
    func_anno <- read_functional_annotation(paths)
    if (all(c("gene_id", "TAP_family") %in% colnames(func_anno))) {
      func_anno <- func_anno[!duplicated(func_anno$gene_id), c("gene_id", "TAP_family"), drop = FALSE]
      fam <- func_anno$TAP_family[match(down_tap_anno$gene_id, func_anno$gene_id)]
      fam <- as.character(fam)
      fam[is.na(fam) | !nzchar(trimws(fam)) | trimws(fam) %in% c("-", "NA")] <- "non_TAP"
      down_tap_anno$TAP_family <- fam
    }
  }
  down_tap_anno$set <- as.character(down_expr_use$set)
  down_tap_anno$set_label <- as.character(down_expr_use$set_label)
  down_tap_hits <- down_tap_anno[down_tap_anno$TAP_family != "non_TAP", , drop = FALSE]
  down_tap_summary <- down_tap_hits %>%
    count(set, set_label, TAP_family, name = "n_genes") %>%
    arrange(set, desc(n_genes), TAP_family)
  write_tsv(down_tap_hits, file.path(out_dir, "triphasy_downregulated_tap_genes.tsv"))
  write_tsv(down_tap_summary, file.path(out_dir, "triphasy_downregulated_tap_family_summary.tsv"))

  down_heatmap_values <- data.frame(gene_id = rownames(down_row_z), down_row_z, check.names = FALSE)
  down_heatmap_values <- cbind(
    down_heatmap_values,
    set = as.character(down_expr_use$set),
    set_label = as.character(down_expr_use$set_label),
    TAP_family = down_tap_anno$TAP_family
  )
  write_tsv(down_heatmap_values, file.path(out_dir, "triphasy_downregulated_expression_heatmap_values.tsv"))

  ha_down_right <- make_tap_row_annotation(down_tap_anno$TAP_family)

  pdf(file.path(fig_dir, "triphasy_downregulated_expression_heatmap.pdf"), width = 8.5, height = 10.5)
  ht_down <- Heatmap(
    down_row_z,
    name = "row z",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#cccccc", "white", "#663333")),
    top_annotation = HeatmapAnnotation(
      stage = stage_ids,
      col = list(stage = stage_cols_present),
      show_annotation_name = FALSE,
      annotation_legend_param = list(stage = list(ncol = 1))
    ),
    left_annotation = rowAnnotation(
      category = down_row_split,
      col = list(category = setNames(unname(down_set_cols[down_set_levels]), unname(down_set_labels[down_set_levels]))),
      show_annotation_name = FALSE
    ),
    right_annotation = ha_down_right,
    row_split = down_row_split,
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
  draw(ht_down, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  down_profile_long <- down_expr_use %>%
    select(gene_id, set, set_label, all_of(mean_cols)) %>%
    pivot_longer(cols = all_of(mean_cols), names_to = "stage_col", values_to = "mean_tpm") %>%
    mutate(
      stage = sub("^mean_tpm_", "", stage_col),
      stage = factor(stage, levels = stage_ids),
      stage_order = as.integer(stage),
      log2_tpm = log2(mean_tpm + 1)
    )

  down_profile_summary <- down_profile_long %>%
    group_by(set, set_label, stage, stage_order) %>%
    summarise(
      n_genes = n(),
      mean_log2_tpm = mean(log2_tpm, na.rm = TRUE),
      se_log2_tpm = sd(log2_tpm, na.rm = TRUE) / sqrt(n()),
      median_log2_tpm = median(log2_tpm, na.rm = TRUE),
      .groups = "drop"
    )

  write_tsv(down_profile_long, file.path(out_dir, "triphasy_downregulated_expression_profile_gene_values.tsv"))
  write_tsv(down_profile_summary, file.path(out_dir, "triphasy_downregulated_expression_profile_summary.tsv"))

  p_down_profile <- ggplot(down_profile_summary, aes(x = stage, y = mean_log2_tpm, group = set_label)) +
    geom_ribbon(
      aes(ymin = mean_log2_tpm - se_log2_tpm, ymax = mean_log2_tpm + se_log2_tpm),
      fill = "grey80",
      alpha = 0.55,
      color = NA
    ) +
    geom_line(color = "black", linewidth = 1.0) +
    geom_point(aes(color = stage), size = 3.0) +
    facet_wrap(~ set_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = stage_cols_present, drop = FALSE) +
    labs(
      x = NULL,
      y = "Mean log2(TPM + 1)",
      title = "Expression Profiles of Triphasy-Downregulated Gene Sets",
      subtitle = "Mean across genes with +/- SE ribbon; points use the TAI/PCA stage palette"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_downregulated_expression_line_profiles.pdf"), p_down_profile, width = 8.4, height = 8.8)
  # Omitted non-manuscript PDF export: triphasy_downregulated_expression_line_profiles_with_legend.pdf
  p_down_profile_overlay <- ggplot(down_profile_summary, aes(x = stage, y = mean_log2_tpm, group = set_label)) +
    geom_ribbon(
      aes(ymin = mean_log2_tpm - se_log2_tpm, ymax = mean_log2_tpm + se_log2_tpm, fill = set_label),
      alpha = 0.16,
      color = NA
    ) +
    geom_line(aes(color = set_label), linewidth = 1.15) +
    geom_point(aes(color = set_label), size = 2.8) +
    scale_color_manual(values = setNames(unname(down_set_cols[down_set_levels]), unname(down_set_labels[down_set_levels])), name = NULL) +
    scale_fill_manual(values = setNames(unname(down_set_cols[down_set_levels]), unname(down_set_labels[down_set_levels])), name = NULL) +
    labs(
      x = NULL,
      y = "Mean log2(TPM + 1)",
      title = "Expression Profiles of Triphasy-Downregulated Gene Sets",
      subtitle = "Mean across genes with +/- SE ribbon"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_downregulated_expression_line_profiles_overlay.pdf"), p_down_profile_overlay, width = 8.8, height = 4.8)

  combined_profile_summary <- bind_rows(
    profile_summary %>%
      mutate(
        direction = "Upregulated",
        generation = case_when(
          set == "GAM_core_biased" ~ "GAM",
          set == "CSP_biased" ~ "CSP",
          set == "TETRA_biased" ~ "TETRA",
          TRUE ~ as.character(set)
        )
      ),
    down_profile_summary %>%
      mutate(
        direction = "Downregulated",
        generation = case_when(
          set == "GAM_core_downregulated" ~ "GAM",
          set == "CSP_downregulated" ~ "CSP",
          set == "TETRA_downregulated" ~ "TETRA",
          TRUE ~ as.character(set)
        )
      )
  ) %>%
    mutate(
      generation = factor(generation, levels = c("GAM", "CSP", "TETRA")),
      direction = factor(direction, levels = c("Upregulated", "Downregulated"))
    ) %>%
    arrange(generation, direction, stage_order)

  write_tsv(combined_profile_summary, file.path(out_dir, "triphasy_up_down_expression_profile_summary.tsv"))

  p_up_down_by_generation <- ggplot(
    combined_profile_summary,
    aes(x = stage, y = mean_log2_tpm, group = direction)
  ) +
    geom_ribbon(
      aes(ymin = mean_log2_tpm - se_log2_tpm, ymax = mean_log2_tpm + se_log2_tpm, fill = direction),
      alpha = 0.14,
      color = NA
    ) +
    geom_line(aes(color = direction), linewidth = 1.15) +
    geom_point(aes(color = stage), size = 2.6) +
    facet_wrap(~ generation, ncol = 1, scales = "free_y") +
    scale_color_manual(
      values = c(
        Upregulated = "#663333",
        Downregulated = "#6A7F8F",
        stage_cols_present
      ),
      breaks = c("Upregulated", "Downregulated"),
      name = NULL
    ) +
    scale_fill_manual(values = c(Upregulated = "#663333", Downregulated = "#6A7F8F"), name = NULL) +
    labs(
      x = NULL,
      y = "Mean log2(TPM + 1)",
      title = "Upregulated and Downregulated Triphasy Gene Sets by Generation",
      subtitle = "Mean across genes with +/- SE ribbon"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave(
    file.path(fig_dir, "triphasy_up_down_expression_line_profiles_by_generation.pdf"),
    p_up_down_by_generation,
    width = 5.9,
    height = 8.8
  )
}

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "triphasy_biased_expression_heatmap_values.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_expression_profile_gene_values.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_expression_profile_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_tap_genes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_tap_family_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_expression_heatmap.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_expression_line_profiles.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_expression_line_profiles_with_legend.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_expression_line_profiles_overlay.pdf"), "\n", sep = "")
if (file.exists(down_membership_file)) {
  cat(" - ", file.path(out_dir, "triphasy_downregulated_expression_heatmap_values.tsv"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "triphasy_downregulated_tap_genes.tsv"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "triphasy_downregulated_tap_family_summary.tsv"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "triphasy_downregulated_expression_profile_gene_values.tsv"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "triphasy_downregulated_expression_profile_summary.tsv"), "\n", sep = "")
  cat(" - ", file.path(fig_dir, "triphasy_downregulated_expression_heatmap.pdf"), "\n", sep = "")
  cat(" - ", file.path(fig_dir, "triphasy_downregulated_expression_line_profiles.pdf"), "\n", sep = "")
  cat(" - ", file.path(fig_dir, "triphasy_downregulated_expression_line_profiles_overlay.pdf"), "\n", sep = "")
  cat(" - ", file.path(out_dir, "triphasy_up_down_expression_profile_summary.tsv"), "\n", sep = "")
  cat(" - ", file.path(fig_dir, "triphasy_up_down_expression_line_profiles_by_generation.pdf"), "\n", sep = "")
}
