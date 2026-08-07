#!/usr/bin/env Rscript
source("R/io.R")
source("R/wgcna_modules.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
expr <- prepare_expression_inputs(paths, design, analysis)
mods <- run_module_inference(expr$counts, analysis)
write_tsv(mods, file.path(paths$output_root, "modules", "module_membership.tsv"))
cat("Module inference step completed.\n")
