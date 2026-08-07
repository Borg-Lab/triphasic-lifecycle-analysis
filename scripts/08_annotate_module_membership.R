#!/usr/bin/env Rscript

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))

in_mod <- file.path(paths$output_root, "modules", "module_membership_labeled.tsv")
if (!file.exists(in_mod)) in_mod <- file.path(paths$output_root, "modules", "module_membership.tsv")
in_anno <- file.path(paths$output_root, "markers", "marker_table.tsv")
out_file <- file.path(paths$output_root, "modules", "module_membership_labeled_annotated.tsv")

if (!file.exists(in_mod)) stop(sprintf("Missing module membership table: %s", in_mod), call. = FALSE)
if (!file.exists(in_anno)) stop(sprintf("Missing annotation table: %s", in_anno), call. = FALSE)

mods <- read.delim(in_mod, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
anno <- read.delim(in_anno, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"gene_id" %in% colnames(mods)) stop("Module table must contain gene_id.", call. = FALSE)
if (!"gene_id" %in% colnames(anno)) stop("Annotation table must contain gene_id.", call. = FALSE)

# Keep all module columns, append non-overlapping annotation columns.
anno_keep <- setdiff(colnames(anno), colnames(mods))
anno_sub <- anno[, c("gene_id", anno_keep), drop = FALSE]
anno_sub <- anno_sub[!duplicated(anno_sub$gene_id), , drop = FALSE]

out <- merge(mods, anno_sub, by = "gene_id", all.x = TRUE, sort = FALSE)

# Keep module columns first.
mod_cols <- colnames(mods)
other_cols <- setdiff(colnames(out), mod_cols)
out <- out[, c(mod_cols, other_cols), drop = FALSE]

write_tsv(out, out_file)
cat("Wrote:\n")
cat(" - ", out_file, "\n", sep = "")
