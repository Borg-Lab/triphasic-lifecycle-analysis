#!/usr/bin/env Rscript

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"), include_all = TRUE)

out_fig_dir <- file.path(paths$output_root, "figures")
out_qc_dir <- file.path(paths$output_root, "qc")
ensure_dir(out_fig_dir)
ensure_dir(out_qc_dir)

rna_root <- dirname(dirname(paths$counts_tsv))
cov_path <- file.path(
  rna_root, "multiqc", "star_salmon", "multiqc_data",
  "mqc_qualimap_gene_coverage_profile_Normalised.txt"
)
if (!file.exists(cov_path)) stop(sprintf("Missing coverage file: %s", cov_path), call. = FALSE)

cov <- read.delim(cov_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!("Sample" %in% colnames(cov))) stop("Coverage file missing 'Sample' column.", call. = FALSE)

pos_cols <- setdiff(colnames(cov), "Sample")
if (length(pos_cols) < 10) stop("Coverage file has too few position columns.", call. = FALSE)
pos_num <- suppressWarnings(as.numeric(pos_cols))
if (any(!is.finite(pos_num))) stop("Coverage position columns are not numeric.", call. = FALSE)

cov <- cov[cov$Sample %in% design$sample_id, c("Sample", pos_cols), drop = FALSE]
if (nrow(cov) == 0) stop("No overlapping samples between coverage file and design.tsv", call. = FALSE)

meta <- design[match(cov$Sample, design$sample_id), , drop = FALSE]
meta$prep_type <- ifelse(meta$tissue == "thallus", "NEBNext_thallus", "SmartSeq3_nonthallus")

mat <- as.matrix(cov[, pos_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- cov$Sample

# 5'/3' and 5'/mid summary metrics
p5 <- mat[, 1]
p50 <- mat[, which.min(abs(pos_num - 50))]
p3 <- mat[, ncol(mat)]
ratio_5to3 <- p5 / p3
ratio_5to50 <- p5 / p50

sum_df <- data.frame(
  sample_id = cov$Sample,
  stage = meta$stage,
  tissue = meta$tissue,
  include = as.logical(meta$include),
  prep_type = meta$prep_type,
  cov_5prime = p5,
  cov_mid = p50,
  cov_3prime = p3,
  ratio_5_to_3 = ratio_5to3,
  ratio_5_to_mid = ratio_5to50,
  stringsAsFactors = FALSE
)
write_tsv(sum_df, file.path(out_qc_dir, "gene_body_coverage_5to3_by_sample.tsv"))

agg <- do.call(rbind, lapply(split(sum_df, sum_df$prep_type), function(d) {
  data.frame(
    prep_type = d$prep_type[1],
    n_samples = nrow(d),
    ratio_5_to_3_median = median(d$ratio_5_to_3, na.rm = TRUE),
    ratio_5_to_3_mean = mean(d$ratio_5_to_3, na.rm = TRUE),
    ratio_5_to_mid_median = median(d$ratio_5_to_mid, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_tsv(agg, file.path(out_qc_dir, "gene_body_coverage_5to3_summary.tsv"))

cols <- c(
  "SmartSeq3_nonthallus" = "#7b3294",
  "NEBNext_thallus" = "#af8e67"
)
prep <- meta$prep_type

ord <- order(prep, meta$stage, meta$replicate, na.last = TRUE)
mat_ord <- mat[ord, , drop = FALSE]
prep_ord <- prep[ord]

pdf(file.path(out_fig_dir, "qc_gene_body_coverage_5to3.pdf"), width = 9.5, height = 7.0, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1) + 0.1)

# Panel A: per-sample curves + group means
plot(
  NA, NA,
  xlim = c(min(pos_num), max(pos_num)),
  ylim = range(mat, na.rm = TRUE),
  xlab = "Gene body percentile (5' to 3')",
  ylab = "Normalized coverage",
  main = "Gene-body coverage profiles"
)
for (i in seq_len(nrow(mat_ord))) {
  lines(pos_num, mat_ord[i, ], col = grDevices::adjustcolor(cols[prep_ord[i]], alpha.f = 0.18), lwd = 1)
}
for (g in names(cols)) {
  idx <- which(prep == g)
  if (length(idx) == 0) next
  m <- colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  lines(pos_num, m, col = cols[g], lwd = 3)
}
legend("topright", legend = names(cols), col = cols, lwd = 3, bty = "n", cex = 0.9)
abline(v = c(0, 50, 99), lty = 3, col = "#bdbdbd")

# Panel B: 5'/3' ratio distribution
boxplot(
  sum_df$ratio_5_to_3 ~ sum_df$prep_type,
  col = cols[levels(factor(sum_df$prep_type))],
  border = "#4d4d4d",
  ylab = "5' / 3' normalized coverage ratio",
  xlab = "",
  main = "5'/3' coverage ratio"
)
stripchart(sum_df$ratio_5_to_3 ~ sum_df$prep_type, vertical = TRUE, method = "jitter", pch = 16, add = TRUE, col = "#333333")
abline(h = 1, lty = 2, col = "#888888")

dev.off()

cat("Wrote:\n")
cat(" - ", file.path(out_fig_dir, "qc_gene_body_coverage_5to3.pdf"), "\n", sep = "")
cat(" - ", file.path(out_qc_dir, "gene_body_coverage_5to3_by_sample.tsv"), "\n", sep = "")
cat(" - ", file.path(out_qc_dir, "gene_body_coverage_5to3_summary.tsv"), "\n", sep = "")
