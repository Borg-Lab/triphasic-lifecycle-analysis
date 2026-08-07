#!/usr/bin/env Rscript
source("R/io.R")
source("R/qc.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
expr <- prepare_expression_inputs(paths, design, analysis)
qc <- build_qc_metrics(expr$counts, expr$tpm, design, analysis)
means <- compute_group_tpm_means(expr$tpm, design, group_col = "stage")

# Append functional annotation columns to group means for easier interpretation.
anno <- read_functional_annotation(paths)
if (nrow(anno) > 0) {
  keep_cols <- intersect(
    c(
      "gene_id", "Preferred_name", "Description", "TAP_family",
      "GOs", "KEGG_ko", "PFAMs", "COG_category"
    ),
    colnames(anno)
  )
  anno_use <- anno[, keep_cols, drop = FALSE]
  means <- merge(means, anno_use, by = "gene_id", all.x = TRUE, sort = FALSE)
}

write_tsv(qc, file.path(paths$output_root, "qc", "sample_qc.tsv"))
write_tsv(means, file.path(paths$output_root, "tables", "group_tpm_means.tsv"))
cat("Prepared expression matrices and QC summaries.\n")
