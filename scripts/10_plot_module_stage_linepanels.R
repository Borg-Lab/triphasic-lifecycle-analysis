#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))

in_z <- file.path(paths$output_root, "modules", "module_stage_enrichment_zscore_matrix.tsv")
in_order <- file.path(paths$output_root, "modules", "module_stage_enrichment_row_order.tsv")
in_map <- file.path(paths$output_root, "modules", "module_label_map.tsv")
out_fig <- file.path(paths$output_root, "figures", "module_stage_enrichment_barpanels_positive.pdf")
out_fig_line <- file.path(paths$output_root, "figures", "module_stage_enrichment_linepanels_positive.pdf")
out_tab <- file.path(paths$output_root, "modules", "module_stage_enrichment_barpanels_values.tsv")

if (!file.exists(in_z)) stop(sprintf("Missing input: %s", in_z), call. = FALSE)
if (!file.exists(in_order)) stop(sprintf("Missing input: %s", in_order), call. = FALSE)
if (!file.exists(in_map)) stop(sprintf("Missing input: %s", in_map), call. = FALSE)

ensure_dir(dirname(out_fig))
ensure_dir(dirname(out_tab))

z <- read.delim(in_z, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
ord <- read.delim(in_order, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
map <- read.delim(in_map, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"module_id" %in% colnames(z)) stop("module_stage_enrichment_zscore_matrix.tsv needs module_id", call. = FALSE)
if (!"module_id" %in% colnames(ord)) stop("module_stage_enrichment_row_order.tsv needs module_id", call. = FALSE)
if (!all(c("module_id", "module_label") %in% colnames(map))) stop("module_label_map.tsv needs module_id + module_label", call. = FALSE)

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
  "Bm3235f_cystocarps_imm" = "#B2ABD2",
  "Bm3235f_cystocarps_mat" = "#8073AC",
  "Bm3235t_carpospores"    = "#5E3C99",
  "Bm3235f_vegtissue"      = "#B55454",
  "Bm3235f_procarps"       = "#7F1D1D",
  "Bm3235t_vegtissue"      = "#C7DFA3",
  "Bm3235t_stichidia"      = "#A8C77A",
  "Bm3235t_tetraspores"    = "#7FA35A",
  "Bm3235m_vegtissue"      = "#7FB8B2",
  "Bm3235m_spermatangia"   = "#3F9188",
  "Bm3235m_spermatia"      = "#0F766E",
  "Bm3235f_pericarp"       = "#D1BBD7",
  "Bm3235c_carposporophyte" = "#3F237A"
)

module_ids <- intersect(ord$module_id, z$module_id)
z <- z[match(module_ids, z$module_id), , drop = FALSE]

lab_lookup <- setNames(map$module_label, map$module_id)
rank_lookup <- if ("module_rank" %in% colnames(map)) setNames(map$module_rank, map$module_id) else setNames(character(0), character(0))
module_label <- unname(lab_lookup[z$module_id])
bad <- is.na(module_label) | !nzchar(module_label)
module_label[bad] <- unname(rank_lookup[z$module_id])[bad]
bad <- is.na(module_label) | !nzchar(module_label)
module_label[bad] <- z$module_id[bad]

stage_cols_present <- intersect(setdiff(colnames(z), "module_id"), stage_order)
stage_other <- setdiff(setdiff(colnames(z), "module_id"), stage_cols_present)
stage_levels <- c(stage_cols_present, stage_other)

long <- data.frame(
  module_id = rep(z$module_id, each = length(stage_levels)),
  module_label = rep(module_label, each = length(stage_levels)),
  stage = rep(stage_levels, times = nrow(z)),
  z = as.numeric(t(as.matrix(z[, stage_levels, drop = FALSE]))),
  stringsAsFactors = FALSE
)

# Positive-only enrichment for plotting.
long$z_pos <- pmax(long$z, 0)

module_levels <- unique(module_label)
long$module_label <- factor(long$module_label, levels = module_levels)
long$stage <- factor(long$stage, levels = stage_levels)

write_tsv(
  long[, c("module_id", "module_label", "stage", "z", "z_pos"), drop = FALSE],
  out_tab
)

n_mod <- length(unique(long$module_label))

p <- ggplot(long, aes(x = stage, y = z_pos, fill = stage)) +
  geom_col(width = 0.8, color = NA) +
  facet_grid(rows = vars(module_label), scales = "fixed", switch = "y") +
  scale_fill_manual(values = stage_cols, breaks = intersect(stage_levels, names(stage_cols)), drop = FALSE) +
  labs(
    x = "Stage",
    y = "Stage enrichment z-score (positive only)",
    fill = "Stage"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    strip.text.y.left = element_text(size = 8, angle = 0),
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.25, "lines"),
    aspect.ratio = 1,
    legend.position = "right"
  )

# Omitted non-manuscript PDF export: ggsave(out_fig, p, width = 9, height = max(10, n_mod * 1.5))

p_line <- ggplot(long, aes(x = stage, y = z_pos, group = 1)) +
  geom_line(color = "#333333", linewidth = 0.7) +
  geom_point(aes(color = stage), size = 5.4, stroke = 0) +
  facet_grid(rows = vars(module_label), scales = "fixed", switch = "y") +
  scale_color_manual(values = stage_cols, breaks = intersect(stage_levels, names(stage_cols)), drop = FALSE) +
  labs(
    x = "Stage",
    y = "Stage enrichment z-score (positive only)",
    color = "Stage"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    strip.text.y.left = element_text(size = 8, angle = 0),
    strip.background = element_blank(),
    panel.spacing.y = grid::unit(0.25, "lines"),
    aspect.ratio = 0.3,
    legend.position = "right"
  )

ggsave(out_fig_line, p_line, width = 15, height = max(10, n_mod * 1.5))

cat("Wrote:\n")
cat(" - ", out_fig, "\n", sep = "")
cat(" - ", out_fig_line, "\n", sep = "")
cat(" - ", out_tab, "\n", sep = "")
