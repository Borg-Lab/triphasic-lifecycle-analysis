#!/usr/bin/env Rscript

source("R/io.R")
source("R/enrichment.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))

means_path <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")
anno_path <- file.path(paths$output_root, "markers", "marker_table.tsv")
out_dir <- file.path(paths$output_root, "enrichment")
ensure_dir(out_dir)

if (!file.exists(means_path)) stop(sprintf("Missing file: %s", means_path), call. = FALSE)
if (!file.exists(anno_path)) stop(sprintf("Missing file: %s", anno_path), call. = FALSE)

means <- read.delim(means_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
anno <- read.delim(anno_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

stage_cols <- grep("^mean_tpm_", colnames(means), value = TRUE)
if (length(stage_cols) == 0) stop("No mean_tpm_ columns in group_tpm_means.tsv", call. = FALSE)

mat <- as.matrix(means[, stage_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- means$gene_id

detect_threshold <- if (!is.null(analysis$tpm_detect_threshold)) as.numeric(analysis$tpm_detect_threshold) else 1
if (!is.finite(detect_threshold)) detect_threshold <- 1

detect <- mat > detect_threshold
n_detect <- rowSums(detect, na.rm = TRUE)
u_idx <- which(n_detect == 1)

if (length(u_idx) == 0) {
  out_empty <- data.frame()
  write_tsv(out_empty, file.path(out_dir, "stage_unique_gene_sets.tsv"))
  write_tsv(out_empty, file.path(out_dir, "stage_unique_functional_enrichment.tsv"))
  write_tsv(out_empty, file.path(out_dir, "stage_unique_functional_enrichment_top.tsv"))
  cat("No stage-unique genes found at threshold > ", detect_threshold, "\n", sep = "")
  quit(save = "no", status = 0)
}

stage_idx <- apply(detect[u_idx, , drop = FALSE], 1, which.max)
stage_group <- gsub("^mean_tpm_", "", stage_cols[stage_idx])
unique_df <- data.frame(
  gene_id = rownames(mat)[u_idx],
  marker_group = stage_group,
  stringsAsFactors = FALSE
)
unique_df <- unique_df[order(unique_df$marker_group, unique_df$gene_id), , drop = FALSE]

keep_cols <- intersect(
  c("gene_id", "Preferred_name", "Description", "TAP_family", "GOs", "KEGG_ko", "PFAMs", "COG_category"),
  colnames(anno)
)
anno_min <- anno[, keep_cols, drop = FALSE]
unique_df_anno <- merge(unique_df, anno_min, by = "gene_id", all.x = TRUE, sort = FALSE)
write_tsv(unique_df_anno, file.path(out_dir, "stage_unique_gene_sets.tsv"))

universe_df <- data.frame(gene_id = rownames(mat), stringsAsFactors = FALSE)
tmp <- merge(universe_df, anno_min, by = "gene_id", all.x = TRUE, sort = FALSE)
tmp$marker_group <- unique_df$marker_group[match(tmp$gene_id, unique_df$gene_id)]
tmp$marker_flag <- !is.na(tmp$marker_group)
tmp$marker_group[!tmp$marker_flag] <- "background"
tmp$marker_class <- "stage_unique"

enrich_df <- run_marker_functional_enrichment(
  marker_df = tmp,
  term_specs = list(
    list(source = "GO", col = "GOs", split = ","),
    list(source = "KEGG_ko", col = "KEGG_ko", split = ","),
    list(source = "PFAM", col = "PFAMs", split = ","),
    list(source = "COG", col = "COG_category", split = ",")
  ),
  min_term_genes = 3L,
  marker_classes = NULL,
  universe_genes = rownames(mat)
)

if (nrow(enrich_df) > 0) {
  enrich_df$global_fdr <- enrich_df$fdr
  split_stage <- split(enrich_df, enrich_df$marker_group)
  enrich_df <- do.call(rbind, lapply(split_stage, function(d) {
    d$fdr <- stats::p.adjust(d$p_value, method = "BH")
    d
  }))
  rownames(enrich_df) <- NULL
}

out_file <- file.path(out_dir, "stage_unique_functional_enrichment.tsv")
write_tsv(enrich_df, out_file)

top_file <- file.path(out_dir, "stage_unique_functional_enrichment_top.tsv")
if (nrow(enrich_df) > 0) {
  ord <- order(enrich_df$fdr, enrich_df$p_value)
  top <- enrich_df[ord, , drop = FALSE]
  top <- do.call(rbind, lapply(split(top, top$marker_group), function(d) head(d, 20)))
  write_tsv(top, top_file)
} else {
  write_tsv(enrich_df, top_file)
}

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "stage_unique_gene_sets.tsv"), "\n", sep = "")
cat(" - ", out_file, "\n", sep = "")
cat(" - ", top_file, "\n", sep = "")
