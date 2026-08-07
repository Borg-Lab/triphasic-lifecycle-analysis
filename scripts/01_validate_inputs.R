#!/usr/bin/env Rscript
source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
contrasts <- read_contrast_table(Sys.getenv("CONTRASTS_CONFIG", unset = "config/contrasts.tsv"))
validate_configs(paths, analysis, design, contrasts)
cat("Input/config validation passed.\n")
