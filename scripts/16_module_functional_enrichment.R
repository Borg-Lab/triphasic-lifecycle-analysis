#!/usr/bin/env Rscript

source("R/io.R")
source("R/enrichment.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_dir <- file.path(paths$output_root, "enrichment")
ensure_dir(out_dir)

module_path <- file.path(paths$output_root, "modules", "module_membership.tsv")
label_map_path <- file.path(paths$output_root, "modules", "module_label_map.tsv")
anno_path <- file.path(paths$output_root, "markers", "marker_table.tsv")

if (!file.exists(module_path)) stop(sprintf("Missing module table: %s", module_path), call. = FALSE)
if (!file.exists(label_map_path)) stop(sprintf("Missing module label map: %s", label_map_path), call. = FALSE)
if (!file.exists(anno_path)) stop(sprintf("Missing annotation table: %s", anno_path), call. = FALSE)

module_df <- read.delim(module_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
label_map <- read.delim(label_map_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
anno_df <- read.delim(anno_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("module_id", "module_label", "module_rank") %in% colnames(label_map))) {
  stop("module_label_map.tsv must contain module_id, module_label, and module_rank.", call. = FALSE)
}

relabel_modules <- function(df, label_df) {
  if (nrow(df) == 0 || !"module_id" %in% colnames(df)) return(df)
  if ("module_id_raw" %in% colnames(df)) return(df)
  keep <- unique(label_df[, c("module_id", "module_rank", "module_label"), drop = FALSE])
  out <- merge(df, keep, by = "module_id", all.x = TRUE, sort = FALSE)
  out$module_id_raw <- out$module_id
  has_label <- !is.na(out$module_label) & nzchar(out$module_label)
  out$module_id[has_label] <- out$module_label[has_label]
  # Reorder columns for readability.
  first_cols <- c("module_id", "module_rank", "module_id_raw")
  other_cols <- setdiff(colnames(out), c(first_cols, "module_label"))
  out <- out[, c(first_cols, other_cols), drop = FALSE]
  out
}

enrich_df <- run_module_functional_enrichment(
  module_df = module_df,
  annotation_df = anno_df,
  term_specs = list(
    list(source = "GO", col = "GOs", split = ","),
    list(source = "KEGG_ko", col = "KEGG_ko", split = ","),
    list(source = "PFAM", col = "PFAMs", split = ","),
    list(source = "COG", col = "COG_category", split = ",")
  ),
  min_term_genes = 5L,
  min_set_size = 10L
)
enrich_df <- relabel_modules(enrich_df, label_map)

out_file <- file.path(out_dir, "module_functional_enrichment.tsv")
write_tsv(enrich_df, out_file)

top_file <- file.path(out_dir, "module_functional_enrichment_top.tsv")
if (nrow(enrich_df) > 0) {
  ord <- order(enrich_df$fdr, enrich_df$p_value)
  top <- enrich_df[ord, , drop = FALSE]
  top <- do.call(rbind, lapply(split(top, top$module_id), function(d) head(d, 20)))
  write_tsv(top, top_file)
} else {
  write_tsv(enrich_df, top_file)
}

strict_file <- file.path(out_dir, "module_functional_enrichment_strict.tsv")
strict_top_file <- file.path(out_dir, "module_functional_enrichment_top_strict.tsv")

if (nrow(enrich_df) > 0) {
  # Remove very broad/uninformative terms and keep robust enrichments.
  generic_go_ids <- c(
    "GO:0008150", "GO:0003674", "GO:0005575",
    "GO:0009987", "GO:0005622", "GO:0005623",
    "GO:0043226", "GO:0043229", "GO:0044424", "GO:0044464"
  )
  strict_df <- enrich_df[
    enrich_df$fdr <= 0.05 &
      enrich_df$overlap >= 3 &
      enrich_df$odds_ratio >= 2 &
      enrich_df$term_genes >= 10 &
      enrich_df$term_genes <= 500 &
      !(enrich_df$source == "GO" & enrich_df$term %in% generic_go_ids),
    ,
    drop = FALSE
  ]
} else {
  strict_df <- enrich_df
}

write_tsv(strict_df, strict_file)

if (nrow(strict_df) > 0) {
  ord <- order(strict_df$fdr, strict_df$p_value)
  strict_top <- strict_df[ord, , drop = FALSE]
  strict_top <- do.call(rbind, lapply(split(strict_top, strict_top$module_id), function(d) head(d, 20)))
  write_tsv(strict_top, strict_top_file)
} else {
  write_tsv(strict_df, strict_top_file)
}

cat("Wrote:\n")
cat(" - ", out_file, "\n", sep = "")
cat(" - ", top_file, "\n", sep = "")
cat(" - ", strict_file, "\n", sep = "")
cat(" - ", strict_top_file, "\n", sep = "")
