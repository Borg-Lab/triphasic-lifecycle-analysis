#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))

z_file <- file.path(paths$output_root, "modules", "module_stage_enrichment_zscore_matrix.tsv")
mod_file <- file.path(paths$output_root, "modules", "module_membership.tsv")
design_file <- "config/design.tsv"
out_file <- file.path(paths$output_root, "modules", "module_label_map.tsv")
mod_labeled_file <- file.path(paths$output_root, "modules", "module_membership_labeled.tsv")

if (!file.exists(z_file)) stop(sprintf("Missing: %s", z_file), call. = FALSE)
if (!file.exists(mod_file)) stop(sprintf("Missing: %s", mod_file), call. = FALSE)
if (!file.exists(design_file)) stop(sprintf("Missing: %s", design_file), call. = FALSE)

z <- read.delim(z_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
mods <- read.delim(mod_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
design <- read.delim(design_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"module_id" %in% colnames(z)) stop("module_stage_enrichment_zscore_matrix.tsv needs module_id", call. = FALSE)
if (!"module_id" %in% colnames(mods)) stop("module_membership.tsv needs module_id", call. = FALSE)
if (!all(c("stage", "generation") %in% colnames(design))) {
  stop("config/design.tsv must contain stage and generation columns", call. = FALSE)
}

zmat <- as.matrix(z[, setdiff(colnames(z), "module_id"), drop = FALSE])
storage.mode(zmat) <- "numeric"
rownames(zmat) <- z$module_id

best_idx <- apply(zmat, 1, which.max)
best_stage <- colnames(zmat)[best_idx]
best_z <- zmat[cbind(seq_len(nrow(zmat)), best_idx)]

top2_idx <- t(apply(zmat, 1, function(x) order(x, decreasing = TRUE)[1:2]))
second_idx <- top2_idx[, 2]
second_stage <- colnames(zmat)[second_idx]
second_z <- zmat[cbind(seq_len(nrow(zmat)), second_idx)]
gap <- as.numeric(best_z - second_z)

stage_map <- unique(design[, c("stage", "generation")])
stage_gen_map <- setNames(stage_map$generation, stage_map$stage)
dominant_generation <- unname(stage_gen_map[best_stage])
second_generation <- unname(stage_gen_map[second_stage])

mod_sizes <- as.data.frame(table(mods$module_id), stringsAsFactors = FALSE)
colnames(mod_sizes) <- c("module_id", "module_size")

base <- data.frame(
  module_id = rownames(zmat),
  dominant_stage = best_stage,
  dominant_generation = dominant_generation,
  dominant_stage_z = as.numeric(best_z),
  second_stage = second_stage,
  second_generation = second_generation,
  second_stage_z = as.numeric(second_z),
  specificity_gap = gap,
  stringsAsFactors = FALSE
)
base <- merge(base, mod_sizes, by = "module_id", all.x = TRUE)
base$module_size[is.na(base$module_size)] <- 0L

# Canonical readable order: module dendrogram order from the same
# module-by-stage z-score matrix used in the module-stage heatmap.
if (nrow(zmat) >= 2) {
  hc <- stats::hclust(stats::dist(zmat, method = "euclidean"), method = "ward.D2")
  dend_order <- rownames(zmat)[hc$order]
} else {
  dend_order <- rownames(zmat)
}
base <- base[match(dend_order, base$module_id), , drop = FALSE]
base <- base[!is.na(base$module_id), , drop = FALSE]

base$module_rank <- sprintf("M%02d", seq_len(nrow(base)))
stage_short <- sub("^Bm3235", "", base$dominant_stage)
stage_short <- gsub("[^A-Za-z0-9]+", "_", stage_short)
stage2_short <- sub("^Bm3235", "", base$second_stage)
stage2_short <- gsub("[^A-Za-z0-9]+", "_", stage2_short)
gen_short <- ifelse(is.na(base$dominant_generation), "NA", base$dominant_generation)
gen2_short <- ifelse(is.na(base$second_generation), "NA", base$second_generation)
base$module_label <- paste0(base$module_rank, "_", stage_short)
base$module_label_generation <- paste0(base$module_rank, "_", gen_short, "_", stage_short)
base$module_label_context <- paste0(
  base$module_rank, "_", gen_short, "_", stage_short, "__alt__", gen2_short, "_", stage2_short
)

write_tsv(base, out_file)

# Propagate label columns into module membership output for downstream readability.
label_cols <- c(
  "module_id", "module_rank", "module_label",
  "module_label_generation", "module_label_context",
  "dominant_stage", "dominant_generation",
  "second_stage", "second_generation",
  "dominant_stage_z", "second_stage_z", "specificity_gap"
)
label_map <- base[, label_cols, drop = FALSE]

# Drop any pre-existing label columns (including stale .x/.y duplicates)
# before re-merging the current label map.
label_bases <- setdiff(label_cols, "module_id")
drop_names <- unique(c(label_bases, paste0(label_bases, ".x"), paste0(label_bases, ".y")))
mods <- mods[, setdiff(colnames(mods), drop_names), drop = FALSE]

mods_labeled <- merge(mods, label_map, by = "module_id", all.x = TRUE, sort = FALSE)

# Keep core columns first if present.
core_first <- c("gene_id", "module_id", "kME", "hub_flag")
core_first <- core_first[core_first %in% colnames(mods_labeled)]
other_cols <- setdiff(colnames(mods_labeled), core_first)
mods_labeled <- mods_labeled[, c(core_first, other_cols), drop = FALSE]

write_tsv(mods_labeled, mod_file)
write_tsv(mods_labeled, mod_labeled_file)
cat("Wrote:\n")
cat(" - ", out_file, "\n", sep = "")
cat(" - ", mod_file, "\n", sep = "")
cat(" - ", mod_labeled_file, "\n", sep = "")
