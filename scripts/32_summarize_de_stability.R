#!/usr/bin/env Rscript

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_dir <- file.path(paths$output_root, "tables")
ensure_dir(out_dir)

de_path <- file.path(paths$output_root, "de", "de_results.tsv")
contr_path <- "config/contrasts.tsv"
tpm_path <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")

if (!file.exists(de_path)) stop(sprintf("Missing: %s", de_path), call. = FALSE)
if (!file.exists(contr_path)) stop(sprintf("Missing: %s", contr_path), call. = FALSE)
if (!file.exists(tpm_path)) stop(sprintf("Missing: %s", tpm_path), call. = FALSE)

de <- read.delim(de_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
contr <- read.delim(contr_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
tpm <- read.delim(tpm_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

req_de <- c("gene_id", "contrast_id", "log2fc", "fdr")
if (!all(req_de %in% colnames(de))) {
  stop("de_results.tsv missing required columns.", call. = FALSE)
}
req_contr <- c("contrast_id", "lfc_threshold", "fdr_threshold")
if (!all(req_contr %in% colnames(contr))) {
  stop("config/contrasts.tsv missing required columns.", call. = FALSE)
}

de2 <- merge(
  de,
  contr[, req_contr],
  by = "contrast_id",
  all.x = TRUE,
  sort = FALSE
)

de2$de_hit <- with(
  de2,
  is.finite(fdr) & is.finite(log2fc) &
    (fdr <= fdr_threshold) &
    (abs(log2fc) >= lfc_threshold)
)

any_de <- tapply(de2$de_hit, de2$gene_id, function(x) any(x, na.rm = TRUE))
tested_any <- tapply(de2$fdr, de2$gene_id, function(x) any(is.finite(x)))
gene_stat <- data.frame(
  gene_id = names(any_de),
  de_any = as.logical(any_de),
  stringsAsFactors = FALSE
)
gene_stat$tested_any <- as.logical(tested_any[gene_stat$gene_id])

tpm_cols <- grep("^mean_tpm_", colnames(tpm), value = TRUE)
if (length(tpm_cols) == 0) stop("No mean_tpm_* columns in group_tpm_means.tsv", call. = FALSE)

mat <- as.matrix(tpm[, tpm_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- tpm$gene_id

# Keep stages with non-NA values.
mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
stage_n <- ncol(mat)

max_tpm <- apply(mat, 1, function(x) max(x, na.rm = TRUE))
min_tpm <- apply(mat, 1, function(x) min(x, na.rm = TRUE))
det1 <- apply(mat, 1, function(x) sum(x > 1, na.rm = TRUE))
det01 <- apply(mat, 1, function(x) sum(x > 0.1, na.rm = TRUE))
l2range <- log2(max_tpm + 1) - log2(min_tpm + 1)
q10 <- apply(mat, 1, function(x) as.numeric(stats::quantile(x, 0.10, na.rm = TRUE, type = 7)))
q90 <- apply(mat, 1, function(x) as.numeric(stats::quantile(x, 0.90, na.rm = TRUE, type = 7)))
iqr <- apply(mat, 1, function(x) IQR(x, na.rm = TRUE))

expr <- data.frame(
  gene_id = rownames(mat),
  max_tpm = max_tpm,
  min_tpm = min_tpm,
  detected_stages_tpm1 = det1,
  detected_stages_tpm0_1 = det01,
  log2_range = l2range,
  q10_tpm = q10,
  q90_tpm = q90,
  iqr_tpm = iqr,
  stringsAsFactors = FALSE
)

allg <- merge(expr, gene_stat, by = "gene_id", all.x = TRUE, sort = FALSE)
allg$de_any[is.na(allg$de_any)] <- FALSE
allg$tested_any[is.na(allg$tested_any)] <- FALSE
allg$class <- ifelse(allg$de_any, "DE_at_least_one_contrast", "stable_nonDE")
allg$expressed_any_tpm1 <- allg$max_tpm > 1
allg$constitutive_tpm1_all_stages <- allg$detected_stages_tpm1 == stage_n
allg$constitutive_tpm0_1_all_stages <- allg$detected_stages_tpm0_1 == stage_n

n_tested <- sum(allg$tested_any)
n_de <- sum(allg$de_any & allg$tested_any)
n_stable <- sum((!allg$de_any) & allg$tested_any)

summary_tbl <- data.frame(
  metric = c(
    "n_genes_in_group_tpm_means",
    "n_genes_tested_in_DE",
    "n_genes_DE_at_least_one_contrast",
    "n_genes_stable_nonDE",
    "pct_DE_of_tested",
    "pct_stable_of_tested",
    "n_expressed_any_tpm_gt1",
    "pct_expressed_any_tpm_gt1",
    "n_constitutive_tpm_gt1_all_stages",
    "n_constitutive_tpm_gt0_1_all_stages",
    "n_constitutive_tpm_gt1_and_stable_nonDE",
    "n_constitutive_tpm_gt1_and_DE_any"
  ),
  value = c(
    nrow(allg),
    n_tested,
    n_de,
    n_stable,
    100 * n_de / max(1, n_tested),
    100 * n_stable / max(1, n_tested),
    sum(allg$expressed_any_tpm1),
    100 * mean(allg$expressed_any_tpm1),
    sum(allg$constitutive_tpm1_all_stages),
    sum(allg$constitutive_tpm0_1_all_stages),
    sum(allg$constitutive_tpm1_all_stages & !allg$de_any),
    sum(allg$constitutive_tpm1_all_stages & allg$de_any)
  ),
  stringsAsFactors = FALSE
)

stats_by_class <- do.call(
  rbind,
  lapply(split(allg, allg$class), function(d) {
    data.frame(
      class = unique(d$class),
      n = nrow(d),
      log2_range_median = median(d$log2_range, na.rm = TRUE),
      log2_range_q25 = as.numeric(stats::quantile(d$log2_range, 0.25, na.rm = TRUE)),
      log2_range_q75 = as.numeric(stats::quantile(d$log2_range, 0.75, na.rm = TRUE)),
      iqr_tpm_median = median(d$iqr_tpm, na.rm = TRUE),
      iqr_tpm_q25 = as.numeric(stats::quantile(d$iqr_tpm, 0.25, na.rm = TRUE)),
      iqr_tpm_q75 = as.numeric(stats::quantile(d$iqr_tpm, 0.75, na.rm = TRUE)),
      detected_stages_tpm1_median = median(d$detected_stages_tpm1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

allg_tested <- allg[allg$tested_any, , drop = FALSE]
stats_by_class_tested <- do.call(
  rbind,
  lapply(split(allg_tested, allg_tested$class), function(d) {
    data.frame(
      class = unique(d$class),
      n = nrow(d),
      log2_range_median = median(d$log2_range, na.rm = TRUE),
      log2_range_q25 = as.numeric(stats::quantile(d$log2_range, 0.25, na.rm = TRUE)),
      log2_range_q75 = as.numeric(stats::quantile(d$log2_range, 0.75, na.rm = TRUE)),
      iqr_tpm_median = median(d$iqr_tpm, na.rm = TRUE),
      iqr_tpm_q25 = as.numeric(stats::quantile(d$iqr_tpm, 0.25, na.rm = TRUE)),
      iqr_tpm_q75 = as.numeric(stats::quantile(d$iqr_tpm, 0.75, na.rm = TRUE)),
      detected_stages_tpm1_median = median(d$detected_stages_tpm1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

const <- allg[allg$constitutive_tpm1_all_stages, , drop = FALSE]
const_summary <- data.frame(
  group = c("constitutive_all", "constitutive_stable_nonDE", "constitutive_DE_any"),
  n = c(nrow(const), sum(!const$de_any), sum(const$de_any)),
  log2_range_median = c(
    median(const$log2_range, na.rm = TRUE),
    median(const$log2_range[!const$de_any], na.rm = TRUE),
    median(const$log2_range[const$de_any], na.rm = TRUE)
  ),
  log2_range_q25 = c(
    as.numeric(stats::quantile(const$log2_range, 0.25, na.rm = TRUE)),
    as.numeric(stats::quantile(const$log2_range[!const$de_any], 0.25, na.rm = TRUE)),
    as.numeric(stats::quantile(const$log2_range[const$de_any], 0.25, na.rm = TRUE))
  ),
  log2_range_q75 = c(
    as.numeric(stats::quantile(const$log2_range, 0.75, na.rm = TRUE)),
    as.numeric(stats::quantile(const$log2_range[!const$de_any], 0.75, na.rm = TRUE)),
    as.numeric(stats::quantile(const$log2_range[const$de_any], 0.75, na.rm = TRUE))
  ),
  iqr_tpm_median = c(
    median(const$iqr_tpm, na.rm = TRUE),
    median(const$iqr_tpm[!const$de_any], na.rm = TRUE),
    median(const$iqr_tpm[const$de_any], na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write_tsv(summary_tbl, file.path(out_dir, "de_stability_summary.tsv"))
write_tsv(stats_by_class, file.path(out_dir, "de_dynamic_range_by_class.tsv"))
write_tsv(stats_by_class_tested, file.path(out_dir, "de_dynamic_range_by_class_tested_only.tsv"))
write_tsv(const_summary, file.path(out_dir, "constitutive_dynamic_range_summary.tsv"))
write_tsv(allg, file.path(out_dir, "gene_de_dynamic_profile.tsv"))

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "de_stability_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "de_dynamic_range_by_class.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "de_dynamic_range_by_class_tested_only.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "constitutive_dynamic_range_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "gene_de_dynamic_profile.tsv"), "\n", sep = "")
