#!/usr/bin/env Rscript
source("R/io.R")
source("R/de.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
contrasts <- read_contrast_table(Sys.getenv("CONTRASTS_CONFIG", unset = "config/contrasts.tsv"))
expr <- prepare_expression_inputs(paths, design, analysis)
de <- run_differential_expression(expr$counts, design, contrasts)
write_tsv(de, file.path(paths$output_root, "de", "de_results.tsv"))
cat("DE step completed.\n")
