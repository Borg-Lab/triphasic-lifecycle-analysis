#!/usr/bin/env Rscript

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))

in_mod <- file.path(paths$output_root, "modules", "module_membership.tsv")
in_map <- file.path(paths$output_root, "modules", "module_label_map.tsv")
out_file <- file.path(paths$output_root, "modules", "module_gene_counts.tsv")

if (!file.exists(in_mod)) stop(sprintf("Missing: %s", in_mod), call. = FALSE)
if (!file.exists(in_map)) stop(sprintf("Missing: %s", in_map), call. = FALSE)

mods <- read.delim(in_mod, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
map <- read.delim(in_map, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("gene_id", "module_id") %in% colnames(mods))) {
  stop("module_membership.tsv must contain gene_id and module_id", call. = FALSE)
}
if (!all(c("module_id", "module_rank", "module_label") %in% colnames(map))) {
  stop("module_label_map.tsv must contain module_id, module_rank, module_label", call. = FALSE)
}

counts <- aggregate(gene_id ~ module_id, data = unique(mods[, c("gene_id", "module_id"), drop = FALSE]), FUN = length)
colnames(counts)[colnames(counts) == "gene_id"] <- "n_genes"

keep <- c(
  "module_id", "module_rank", "module_label", "module_label_generation", "module_label_context",
  "dominant_stage", "dominant_generation", "second_stage", "second_generation", "specificity_gap"
)
keep <- keep[keep %in% colnames(map)]
out <- merge(counts, unique(map[, keep, drop = FALSE]), by = "module_id", all.x = TRUE, sort = FALSE)

if ("module_rank" %in% colnames(out)) {
  out$module_rank_num <- suppressWarnings(as.integer(sub("^M", "", out$module_rank)))
  out <- out[order(out$module_rank_num, out$module_id), , drop = FALSE]
  out$module_rank_num <- NULL
} else {
  out <- out[order(-out$n_genes, out$module_id), , drop = FALSE]
}

write_tsv(out, out_file)
cat("Wrote:\n")
cat(" - ", out_file, "\n", sep = "")
