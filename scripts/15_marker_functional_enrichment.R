#!/usr/bin/env Rscript

source("R/io.R")
source("R/enrichment.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_dir <- file.path(paths$output_root, "enrichment")
ensure_dir(out_dir)

marker_path <- file.path(paths$output_root, "markers", "marker_table.tsv")
if (!file.exists(marker_path)) {
  stop(sprintf("Missing marker table: %s", marker_path), call. = FALSE)
}

marker_df <- read.delim(marker_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

enrich_df <- run_marker_functional_enrichment(
  marker_df = marker_df,
  term_specs = list(
    list(source = "GO", col = "GOs", split = ","),
    list(source = "KEGG_ko", col = "KEGG_ko", split = ","),
    list(source = "PFAM", col = "PFAMs", split = ","),
    list(source = "COG", col = "COG_category", split = ",")
  ),
  min_term_genes = 3L,
  marker_classes = NULL
)

out_file <- file.path(out_dir, "marker_functional_enrichment.tsv")
write_tsv(enrich_df, out_file)

top_file <- file.path(out_dir, "marker_functional_enrichment_top.tsv")
if (nrow(enrich_df) > 0) {
  ord <- order(enrich_df$fdr, enrich_df$p_value)
  top <- enrich_df[ord, , drop = FALSE]
  top <- do.call(rbind, lapply(split(top, top$marker_group), function(d) head(d, 20)))
  write_tsv(top, top_file)
} else {
  write_tsv(enrich_df, top_file)
}

cat("Wrote:\n")
cat(" - ", out_file, "\n", sep = "")
cat(" - ", top_file, "\n", sep = "")
