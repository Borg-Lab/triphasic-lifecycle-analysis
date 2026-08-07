#!/usr/bin/env Rscript

source("R/io.R")
source("R/tai.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
contrasts <- read_contrast_table(Sys.getenv("CONTRASTS_CONFIG", unset = "config/contrasts.tsv"))

tab_dir <- file.path(paths$output_root, "tables")
ensure_dir(tab_dir)

read_tsv <- function(path) safe_read_table(path)
normalize_df_gene_id <- function(x) {
  if (nrow(x) == 0 || !"gene_id" %in% colnames(x)) return(x)
  x$gene_id <- normalize_gene_id(x$gene_id)
  x <- x[!is.na(x$gene_id) & nzchar(x$gene_id), , drop = FALSE]
  x
}
first_nonmissing <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) NA_character_ else as.character(x[1])
}
collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = ";")
}
logical_or <- function(x) any(as.logical(x), na.rm = TRUE)

means_path <- file.path(tab_dir, "group_tpm_means.tsv")
if (!file.exists(means_path)) stop("Missing ", means_path, call. = FALSE)
out <- normalize_df_gene_id(read_tsv(means_path))
if (nrow(out) == 0) stop("group_tpm_means.tsv is empty.", call. = FALSE)
out <- out[!duplicated(out$gene_id), , drop = FALSE]
base_gene_order <- out$gene_id

add_cols <- function(df, cols, prefix = NULL) {
  df <- normalize_df_gene_id(df)
  if (nrow(df) == 0) return(invisible(NULL))
  cols <- intersect(cols, colnames(df))
  if (!"gene_id" %in% cols) cols <- c("gene_id", cols)
  df <- df[, unique(cols), drop = FALSE]
  df <- df[!duplicated(df$gene_id), , drop = FALSE]
  if (!is.null(prefix)) {
    rename <- setdiff(colnames(df), "gene_id")
    colnames(df)[match(rename, colnames(df))] <- paste0(prefix, rename)
  }
  assign("out", merge(get("out", envir = parent.frame()), df, by = "gene_id", all.x = TRUE, sort = FALSE), envir = parent.frame())
  invisible(NULL)
}

# Coordinates and genomic context.
coords <- read_coordinates(paths)
add_cols(coords, c("gene_id", "chr", "start", "end", "strand"), prefix = "coord_")

# Stage-unique genes: detected only in one stage above the TPM threshold.
stage_unique <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "enrichment", "stage_unique_gene_sets.tsv")))
if (nrow(stage_unique) > 0) {
  su <- aggregate(marker_group ~ gene_id, stage_unique, collapse_unique)
  colnames(su) <- c("gene_id", "stage_unique_stage")
  su$is_stage_unique_gene <- TRUE
  out <- merge(out, su, by = "gene_id", all.x = TRUE, sort = FALSE)
}
out$is_stage_unique_gene[is.na(out$is_stage_unique_gene)] <- FALSE

# Marker genes and strict marker calls.
markers <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "markers", "marker_table.tsv")))
marker_cols <- c(
  "gene_id", "marker_group", "tau", "best_log2fc", "best_fdr", "de_available", "de_testable",
  "marker_specificity_flag", "marker_flag", "marker_class", "strict_tau", "strict_detected_stages",
  "strict_max_tpm", "strict_marker_group", "strict_specificity_flag", "strict_marker_flag",
  "strict_marker_class", "orig_marker_group", "orig_tau", "orig_marker_flag", "orig_marker_class"
)
add_cols(markers, marker_cols, prefix = "marker_")

# Generation/stage bias summary classification.
gen_class <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "generation", "generation_stage_bias_gene_classification.tsv")))
add_cols(gen_class, c("gene_id", "generation_biased", "stage_biased", "bias_class"), prefix = "genstage_")

# Triphasy biased/repressed sets and combined generation identity membership.
trip_biased <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "generation", "triphasy_biased_gene_set_membership.tsv")))
if (nrow(trip_biased) > 0) {
  x <- aggregate(set ~ gene_id, trip_biased, collapse_unique)
  colnames(x) <- c("gene_id", "triphasy_biased_sets")
  x$is_generation_biased_gene <- TRUE
  out <- merge(out, x, by = "gene_id", all.x = TRUE, sort = FALSE)
}
out$is_generation_biased_gene[is.na(out$is_generation_biased_gene)] <- FALSE

trip_down <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "generation", "triphasy_downregulated_gene_set_membership.tsv")))
if (nrow(trip_down) > 0) {
  x <- aggregate(set ~ gene_id, trip_down, collapse_unique)
  colnames(x) <- c("gene_id", "triphasy_repressed_sets")
  x$is_generation_repressed_gene <- TRUE
  out <- merge(out, x, by = "gene_id", all.x = TRUE, sort = FALSE)
}
out$is_generation_repressed_gene[is.na(out$is_generation_repressed_gene)] <- FALSE

trip_combined <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "sex_bias", "triphasy_sex_bias_gene_set_membership.tsv")))
if (nrow(trip_combined) > 0) {
  generation_sets <- trip_combined[!trip_combined$set %in% c("female_biased", "male_biased"), , drop = FALSE]
  if (nrow(generation_sets) > 0) {
    x <- aggregate(set ~ gene_id, generation_sets, collapse_unique)
    colnames(x) <- c("gene_id", "generation_identity_sets")
    out <- merge(out, x, by = "gene_id", all.x = TRUE, sort = FALSE)
  }
}

# Target specificity for generation-biased and repressed sets.
trip_spec <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "generation", "triphasy_biased_target_specificity.tsv")))
if (nrow(trip_spec) > 0) {
  keep <- c("gene_id", "target_group", "target_max_tpm", "opposing_max_tpm", "target_detected", "opposing_detected", "specificity_class")
  trip_spec <- trip_spec[, intersect(keep, colnames(trip_spec)), drop = FALSE]
  trip_spec <- trip_spec[!duplicated(trip_spec$gene_id), , drop = FALSE]
  names(trip_spec)[names(trip_spec) != "gene_id"] <- paste0("triphasy_biased_", names(trip_spec)[names(trip_spec) != "gene_id"])
  out <- merge(out, trip_spec, by = "gene_id", all.x = TRUE, sort = FALSE)
}

down_off <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "generation", "triphasy_downregulated_target_off_status.tsv")))
if (nrow(down_off) > 0) {
  keep <- c("gene_id", "target_group", "target_max_tpm", "opposing_max_tpm", "all_max_tpm", "target_generation_off", "target_off_opposing_detected", "all_stages_off")
  down_off <- down_off[, intersect(keep, colnames(down_off)), drop = FALSE]
  down_off <- down_off[!duplicated(down_off$gene_id), , drop = FALSE]
  names(down_off)[names(down_off) != "gene_id"] <- paste0("triphasy_repressed_", names(down_off)[names(down_off) != "gene_id"])
  out <- merge(out, down_off, by = "gene_id", all.x = TRUE, sort = FALSE)
}

# Sex-biased genes and specificity classes.
sbg <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "sex_bias", "sex_biased_gene_sets.tsv")))
if (nrow(sbg) > 0) {
  x <- sbg[, intersect(c("gene_id", "set", "log2fc", "fdr", "base_mean"), colnames(sbg)), drop = FALSE]
  x <- x[!duplicated(x$gene_id), , drop = FALSE]
  colnames(x) <- sub("^set$", "sex_bias_set", colnames(x))
  colnames(x) <- sub("^log2fc$", "sex_bias_log2fc", colnames(x))
  colnames(x) <- sub("^fdr$", "sex_bias_fdr", colnames(x))
  colnames(x) <- sub("^base_mean$", "sex_bias_base_mean", colnames(x))
  x$is_sex_biased_gene <- TRUE
  out <- merge(out, x, by = "gene_id", all.x = TRUE, sort = FALSE)
}
out$is_sex_biased_gene[is.na(out$is_sex_biased_gene)] <- FALSE

sex_spec <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "sex_bias", "sex_biased_target_specificity.tsv")))
if (nrow(sex_spec) > 0) {
  keep <- c(
    "gene_id", "region", "genomic_compartment", "specificity_class", "expression_specificity",
    "lifecycle_specificity_class", "target_core_lifecycle_specific", "target_detected",
    "opposing_detected", "non_gametophyte_detected", "outside_target_core_detected",
    "female_core_max_tpm", "male_core_max_tpm", "non_gametophyte_max_tpm"
  )
  sex_spec <- sex_spec[, intersect(keep, colnames(sex_spec)), drop = FALSE]
  sex_spec <- sex_spec[!duplicated(sex_spec$gene_id), , drop = FALSE]
  names(sex_spec)[names(sex_spec) != "gene_id"] <- paste0("sex_bias_", names(sex_spec)[names(sex_spec) != "gene_id"])
  out <- merge(out, sex_spec, by = "gene_id", all.x = TRUE, sort = FALSE)
}

# WGCNA module membership.
modules <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "modules", "module_membership_labeled.tsv")))
module_cols <- c(
  "gene_id", "module_id", "module_rank", "module_label", "module_label_generation", "module_label_context",
  "kME", "hub_flag", "dominant_stage", "dominant_generation", "second_stage", "second_generation",
  "dominant_stage_z", "second_stage_z", "specificity_gap"
)
add_cols(modules, module_cols, prefix = "module_")
out$is_module_member <- !is.na(out$module_module_label) & nzchar(as.character(out$module_module_label))
out$is_module_hub <- as.logical(out$module_hub_flag)
out$is_module_hub[is.na(out$is_module_hub)] <- FALSE

# DE-at-least-one contrast summary.
de <- normalize_df_gene_id(read_tsv(file.path(paths$output_root, "de", "de_results.tsv")))
if (nrow(de) > 0 && all(c("contrast_id", "log2fc", "fdr") %in% colnames(de))) {
  cthr <- contrasts[, c("contrast_id", "lfc_threshold", "fdr_threshold"), drop = FALSE]
  de <- merge(de, cthr, by = "contrast_id", all.x = TRUE, sort = FALSE)
  de$lfc_threshold[!is.finite(de$lfc_threshold)] <- 1
  de$fdr_threshold[!is.finite(de$fdr_threshold)] <- 0.05
  de$de_direction <- NA_character_
  de$de_direction[is.finite(de$fdr) & is.finite(de$log2fc) & de$fdr <= de$fdr_threshold & de$log2fc >= de$lfc_threshold] <- "up"
  de$de_direction[is.finite(de$fdr) & is.finite(de$log2fc) & de$fdr <= de$fdr_threshold & de$log2fc <= -de$lfc_threshold] <- "down"
  de_sig <- de[!is.na(de$de_direction), c("gene_id", "contrast_id", "de_direction"), drop = FALSE]
  if (nrow(de_sig) > 0) {
    de_sig$contrast_token <- paste0(de_sig$contrast_id, ":", de_sig$de_direction)
    de_sum <- data.frame(
      gene_id = names(table(de_sig$gene_id)),
      n_de_contrasts = as.integer(table(de_sig$gene_id)),
      stringsAsFactors = FALSE
    )
    de_con <- aggregate(contrast_token ~ gene_id, de_sig, collapse_unique)
    colnames(de_con) <- c("gene_id", "de_contrasts")
    de_sum <- merge(de_sum, de_con, by = "gene_id", all.x = TRUE, sort = FALSE)
    de_sum$is_de_in_any_contrast <- TRUE
    out <- merge(out, de_sum, by = "gene_id", all.x = TRUE, sort = FALSE)
  }
}
out$is_de_in_any_contrast[is.na(out$is_de_in_any_contrast)] <- FALSE
out$n_de_contrasts[is.na(out$n_de_contrasts)] <- 0L

# Dynamic expression class.
dyn <- normalize_df_gene_id(read_tsv(file.path(tab_dir, "constitutive_dynamic_gene_metrics.tsv")))
dyn_cols <- c("gene_id", "prevalence", "detected_n", "dynamic_range", "cv_log", "max_tpm", "min_tpm", "dominant_stage", "class")
add_cols(dyn, dyn_cols, prefix = "dynamic_")

# Phylostratum / gene age.
age <- read_gene_age_table(paths$gene_age_tsv)
age$gene_id <- normalize_gene_id(age$gene_id)
age <- age[!is.na(age$gene_id) & nzchar(age$gene_id), , drop = FALSE]
age$rank <- suppressWarnings(as.integer(age$rank))
age <- age[is.finite(age$rank), , drop = FALSE]
age <- age[!duplicated(age$gene_id), , drop = FALSE]
if (nrow(age) > 0) {
  tax_col <- if ("taxon" %in% colnames(age)) "taxon" else if ("age_name" %in% colnames(age)) "age_name" else NA_character_
  age_out <- data.frame(
    gene_id = age$gene_id,
    phylostratum_rank = age$rank,
    phylostratum_label = if (!is.na(tax_col)) paste0("PS", age$rank, ": ", age[[tax_col]]) else paste0("PS", age$rank),
    stringsAsFactors = FALSE
  )
  out <- merge(out, age_out, by = "gene_id", all.x = TRUE, sort = FALSE)
}

# High-level compact flags for quick filtering.
out$is_marker_gene <- as.logical(out$marker_marker_flag)
out$is_marker_gene[is.na(out$is_marker_gene)] <- FALSE
out$is_strict_marker_gene <- as.logical(out$marker_strict_marker_flag)
out$is_strict_marker_gene[is.na(out$is_strict_marker_gene)] <- FALSE
out$is_generation_identity_gene <- out$is_generation_biased_gene | out$is_generation_repressed_gene
out$is_identity_or_sbg_gene <- out$is_generation_identity_gene | out$is_sex_biased_gene

rename_map <- c(
  coord_chr = "chr",
  coord_start = "start",
  coord_end = "end",
  coord_strand = "strand",
  marker_marker_group = "marker_group",
  marker_marker_class = "marker_class",
  marker_marker_flag = "marker_flag_raw",
  marker_strict_marker_group = "strict_marker_group",
  marker_strict_marker_class = "strict_marker_class",
  marker_strict_marker_flag = "strict_marker_flag_raw",
  marker_orig_marker_group = "orig_marker_group",
  marker_orig_marker_flag = "orig_marker_flag",
  marker_orig_marker_class = "orig_marker_class",
  module_module_id = "module_id",
  module_module_rank = "module_rank",
  module_module_label = "module_label",
  module_module_label_generation = "module_label_generation",
  module_module_label_context = "module_label_context",
  module_hub_flag = "module_hub_flag_raw",
  module_dominant_stage = "module_dominant_stage",
  module_second_stage = "module_second_stage",
  dynamic_dynamic_range = "dynamic_range"
)
for (old_name in names(rename_map)) {
  if (old_name %in% colnames(out)) colnames(out)[colnames(out) == old_name] <- rename_map[[old_name]]
}

# Restore group TPM order and put interpretation columns before expression values.
out <- out[match(base_gene_order, out$gene_id), , drop = FALSE]
mean_cols <- grep("^mean_tpm_", colnames(out), value = TRUE)
annotation_cols <- intersect(c("Preferred_name", "Description", "TAP_family", "GOs", "KEGG_ko", "PFAMs", "COG_category"), colnames(out))
front_cols <- intersect(c(
  "gene_id",
  "chr", "start", "end", "strand",
  "is_stage_unique_gene", "stage_unique_stage",
  "is_marker_gene", "marker_group", "marker_class", "marker_tau", "marker_best_log2fc", "marker_best_fdr",
  "is_strict_marker_gene", "strict_marker_group", "strict_marker_class", "marker_strict_tau",
  "genstage_generation_biased", "genstage_stage_biased", "genstage_bias_class",
  "is_generation_biased_gene", "triphasy_biased_sets", "is_generation_repressed_gene", "triphasy_repressed_sets", "generation_identity_sets",
  "triphasy_biased_target_group", "triphasy_biased_specificity_class",
  "triphasy_repressed_target_group", "triphasy_repressed_target_generation_off", "triphasy_repressed_all_stages_off",
  "is_sex_biased_gene", "sex_bias_set", "sex_bias_log2fc", "sex_bias_fdr", "sex_bias_region", "sex_bias_genomic_compartment", "sex_bias_specificity_class", "sex_bias_lifecycle_specificity_class",
  "module_rank", "module_label", "module_label_generation", "module_id", "module_kME", "is_module_member", "is_module_hub", "module_dominant_stage", "module_second_stage",
  "is_de_in_any_contrast", "n_de_contrasts", "de_contrasts",
  "dynamic_class", "dynamic_prevalence", "dynamic_detected_n", "dynamic_range", "dynamic_cv_log", "dynamic_dominant_stage",
  "phylostratum_rank", "phylostratum_label"
), colnames(out))
other_cols <- setdiff(colnames(out), c(front_cols, mean_cols, annotation_cols))
out <- out[, c(front_cols, mean_cols, annotation_cols, other_cols), drop = FALSE]

out_path <- file.path(tab_dir, "group_tpm_means_annotated.tsv")
write_tsv(out, out_path)

simple_cols <- intersect(c(
  "gene_id",
  "stage_unique_stage",
  "marker_group",
  "strict_marker_group",
  "genstage_bias_class",
  "generation_identity_sets",
  "triphasy_biased_sets",
  "triphasy_repressed_sets",
  "sex_bias_set",
  "module_label",
  "module_rank",
  "is_module_hub",
  "dynamic_class",
  "phylostratum_label",
  "TAP_family",
  "Preferred_name"
), colnames(out))
simple <- out[, c(simple_cols, mean_cols), drop = FALSE]
simple_path <- file.path(tab_dir, "group_tpm_means_subgroups.tsv")
write_tsv(simple, simple_path)

summary <- data.frame(
  metric = c(
    "genes_total", "stage_unique_genes", "marker_genes", "strict_marker_genes", "generation_biased_genes",
    "generation_repressed_genes", "generation_identity_genes", "sex_biased_genes", "module_members",
    "module_hubs", "de_any_contrast", "genes_with_phylostratum"
  ),
  value = c(
    nrow(out),
    sum(out$is_stage_unique_gene, na.rm = TRUE),
    sum(out$is_marker_gene, na.rm = TRUE),
    sum(out$is_strict_marker_gene, na.rm = TRUE),
    sum(out$is_generation_biased_gene, na.rm = TRUE),
    sum(out$is_generation_repressed_gene, na.rm = TRUE),
    sum(out$is_generation_identity_gene, na.rm = TRUE),
    sum(out$is_sex_biased_gene, na.rm = TRUE),
    sum(out$is_module_member, na.rm = TRUE),
    sum(out$is_module_hub, na.rm = TRUE),
    sum(out$is_de_in_any_contrast, na.rm = TRUE),
    sum(!is.na(out$phylostratum_rank), na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
summary_path <- file.path(tab_dir, "group_tpm_means_annotated_summary.tsv")
write_tsv(summary, summary_path)

cat("Wrote:\n")
cat(" - ", out_path, "\n", sep = "")
cat(" - ", simple_path, "\n", sep = "")
cat(" - ", summary_path, "\n", sep = "")
