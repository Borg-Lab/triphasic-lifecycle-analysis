#!/usr/bin/env Rscript
source("R/io.R")
source("R/qc.R")
source("R/tai.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
expr <- prepare_expression_inputs(paths, design, analysis)
means <- compute_group_tpm_means(expr$tpm, design, group_col = "stage")
tai <- run_tai_analysis(means, paths)
write_tsv(tai, file.path(paths$output_root, "tai", "tai_summary.tsv"))
run_tai_diagnostics(means, paths, file.path(paths$output_root, "tai"))

# Plot-ready TAI table with cleaned stage labels and metadata
if (nrow(tai) > 0) {
  tai_plot <- tai
  tai_plot$stage_id <- sub("^mean_tpm_", "", tai_plot$stage)

  ord_cols <- order_stage_columns(tai_plot$stage)
  ord_map <- setNames(seq_along(ord_cols), ord_cols)
  tai_plot$stage_order <- as.integer(ord_map[tai_plot$stage])

  design_meta <- unique(design[, intersect(c("stage", "generation", "tissue", "sex"), colnames(design)), drop = FALSE])
  if ("stage" %in% colnames(design_meta)) {
    names(design_meta)[names(design_meta) == "stage"] <- "stage_id"
    tai_plot <- merge(tai_plot, design_meta, by = "stage_id", all.x = TRUE, sort = FALSE)
  }

  tai_plot <- tai_plot[order(tai_plot$stage_order), , drop = FALSE]
  tai_plot$stage_label <- tai_plot$stage_id
  tai_plot <- tai_plot[, c(
    "stage_order", "stage", "stage_id", "stage_label",
    "generation", "tissue", "sex", "TAI", "n_genes"
  )[c(
    "stage_order", "stage", "stage_id", "stage_label",
    "generation", "tissue", "sex", "TAI", "n_genes"
  ) %in% colnames(tai_plot)], drop = FALSE]

  write_tsv(tai_plot, file.path(paths$output_root, "tai", "tai_plot_input.tsv"))
}

cat("TAI step completed.\n")
