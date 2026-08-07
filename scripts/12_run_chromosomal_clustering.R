#!/usr/bin/env Rscript
source("R/io.R")
source("R/chrom_cluster.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
mods <- safe_read_table(file.path(paths$output_root, "modules", "module_membership.tsv"))
coords <- read_coordinates(paths)
coords <- filter_coordinates(coords, analysis)
res <- run_chromosomal_clustering(mods, coords, analysis)
write_tsv(res$summary, file.path(paths$output_root, "chrom_cluster", "chrom_cluster_summary.tsv"))
write_tsv(res$members, file.path(paths$output_root, "chrom_cluster", "chrom_cluster_members.tsv"))
write_tsv(res$blocks, file.path(paths$output_root, "chrom_cluster", "chrom_proximity_blocks.tsv"))
cat("Chromosomal clustering step completed.\n")
