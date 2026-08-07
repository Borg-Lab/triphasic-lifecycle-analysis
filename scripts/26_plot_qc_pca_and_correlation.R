#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

source("R/io.R")
source("R/qc.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))

design <- design[as.logical(design$include), , drop = FALSE]

expr <- prepare_expression_inputs(paths, design, analysis)
counts <- expr$counts

out_qc_dir <- file.path(paths$output_root, "qc")
out_fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_qc_dir)
ensure_dir(out_fig_dir)

# --------------------------
# 1) MultiQC summary parsing
# --------------------------
rna_root <- dirname(dirname(paths$counts_tsv))
multiqc_general <- file.path(rna_root, "multiqc", "star_salmon", "multiqc_data", "multiqc_general_stats.txt")
if (!file.exists(multiqc_general)) {
  stop(sprintf("Missing MultiQC general stats file: %s", multiqc_general), call. = FALSE)
}

mq <- read.delim(multiqc_general, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Sample" %in% colnames(mq)) {
  stop("Column 'Sample' not found in multiqc_general_stats.txt", call. = FALSE)
}

# Keep only exact sample IDs from design to avoid read-level entries like _1/_2 rows.
mq <- mq[mq$Sample %in% design$sample_id, , drop = FALSE]

pick_col <- function(df, name) {
  if (name %in% colnames(df)) return(df[[name]])
  rep(NA_real_, nrow(df))
}

mq_small <- data.frame(
  sample_id = mq$Sample,
  fastqc_raw_total_sequences = as.numeric(pick_col(mq, "FastQC (raw)_mqc-generalstats-fastqc_raw-total_sequences")),
  fastqc_trimmed_total_sequences = as.numeric(pick_col(mq, "FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-total_sequences")),
  samtools_raw_total_sequences = as.numeric(pick_col(mq, "Samtools_mqc-generalstats-samtools-raw_total_sequences")),
  star_uniquely_mapped = as.numeric(pick_col(mq, "STAR_mqc-generalstats-star-uniquely_mapped")),
  star_uniquely_mapped_percent = as.numeric(pick_col(mq, "STAR_mqc-generalstats-star-uniquely_mapped_percent")),
  salmon_num_mapped = as.numeric(pick_col(mq, "Salmon_mqc-generalstats-salmon-num_mapped")),
  salmon_percent_mapped = as.numeric(pick_col(mq, "Salmon_mqc-generalstats-salmon-percent_mapped")),
  samtools_reads_mapped = as.numeric(pick_col(mq, "Samtools_mqc-generalstats-samtools-reads_mapped")),
  samtools_reads_mapped_percent = as.numeric(pick_col(mq, "Samtools_mqc-generalstats-samtools-reads_mapped_percent")),
  picard_percent_duplication = as.numeric(pick_col(mq, "Picard_mqc-generalstats-picard-PERCENT_DUPLICATION")),
  dupradar_intercept = as.numeric(pick_col(mq, "DupRadar_mqc-generalstats-dupradar-dupRadar_intercept")),
  sortmerna_rrna_percent = as.numeric(pick_col(mq, "SortMeRNA_mqc-generalstats-sortmerna-rRNA_pct")),
  stringsAsFactors = FALSE
)

qc_base <- build_qc_metrics(expr$counts, expr$tpm, design, analysis)
qc_merged <- merge(design, qc_base[, c("sample_id", "library_size", "detected_genes")], by = "sample_id", all.x = TRUE)
qc_merged <- merge(qc_merged, mq_small, by = "sample_id", all.x = TRUE)

qc_merged <- qc_merged[order(qc_merged$stage, qc_merged$replicate), , drop = FALSE]
write_tsv(qc_merged, file.path(out_qc_dir, "nfcore_qc_summary.tsv"))

num_cols <- c(
  "library_size", "detected_genes",
  "fastqc_raw_total_sequences", "fastqc_trimmed_total_sequences",
  "samtools_raw_total_sequences", "samtools_reads_mapped", "samtools_reads_mapped_percent",
  "star_uniquely_mapped", "star_uniquely_mapped_percent",
  "salmon_num_mapped", "salmon_percent_mapped",
  "picard_percent_duplication", "dupradar_intercept", "sortmerna_rrna_percent"
)
num_cols <- intersect(num_cols, colnames(qc_merged))

stage_split <- split(qc_merged, qc_merged$stage)
stage_summary <- do.call(
  rbind,
  lapply(names(stage_split), function(stg) {
    d <- stage_split[[stg]]
    out <- data.frame(stage = stg, n_samples = nrow(d), stringsAsFactors = FALSE)
    for (nc in num_cols) {
      out[[paste0(nc, "_mean")]] <- mean(as.numeric(d[[nc]]), na.rm = TRUE)
      out[[paste0(nc, "_sd")]] <- stats::sd(as.numeric(d[[nc]]), na.rm = TRUE)
    }
    out
  })
)
write_tsv(stage_summary, file.path(out_qc_dir, "nfcore_qc_by_stage.tsv"))

# --------------------------
# 2) DESeq2 rlog PCA plot
# --------------------------
if (nrow(counts) < 2) {
  stop("Not enough genes in counts matrix for PCA.", call. = FALSE)
}

samples <- intersect(design$sample_id, colnames(counts))
if (length(samples) < 2) {
  stop("Not enough samples for PCA.", call. = FALSE)
}

counts <- counts[, c("gene_id", samples), drop = FALSE]
count_mat <- as.matrix(counts[, -1, drop = FALSE])
storage.mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
rownames(count_mat) <- counts$gene_id

coldata <- design[match(samples, design$sample_id), , drop = FALSE]
rownames(coldata) <- coldata$sample_id

stage_order <- c(
  "Bm3235f_vegtissue",
  "Bm3235f_procarps",
  "Bm3235m_vegtissue",
  "Bm3235m_spermatangia",
  "Bm3235m_spermatia",
  "Bm3235f_pericarp",
  "Bm3235f_cystocarps_imm",
  "Bm3235f_cystocarps_mat",
  "Bm3235c_carposporophyte",
  "Bm3235t_carpospores",
  "Bm3235t_vegtissue",
  "Bm3235t_stichidia",
  "Bm3235t_tetraspores"
)

stage_levels <- c(intersect(stage_order, unique(coldata$stage)), setdiff(unique(coldata$stage), stage_order))
coldata$stage <- factor(coldata$stage, levels = stage_levels)
coldata$generation <- factor(coldata$generation)

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ stage)
rld <- rlog(dds, blind = FALSE)
pca_df <- plotPCA(rld, intgroup = c("stage", "generation"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"), 1)

# grouped palette per request:
# - cystocarps + carpospores: purple
# - female tissues (veg + procarps): maroon
# - tetrasporophyte tissues: straw green
# - male tissues: teal shades
stage_cols <- c(
  "Bm3235f_cystocarps_imm" = "#B2ABD2", # purple light
  "Bm3235f_cystocarps_mat" = "#8073AC", # purple mid
  "Bm3235f_pericarp"       = "#C77C7C", # pericarp pink
  "Bm3235c_carposporophyte" = "#3F1F78", # carposporophyte purple deep
  "Bm3235t_carpospores"    = "#5E3C99", # purple dark
  "Bm3235f_vegtissue"      = "#B55454", # maroon light
  "Bm3235f_procarps"       = "#7F1D1D", # maroon dark
  "Bm3235t_vegtissue"      = "#C7DFA3", # straw green light
  "Bm3235t_stichidia"      = "#A8C77A", # straw green mid
  "Bm3235t_tetraspores"    = "#7FA35A", # straw green dark
  "Bm3235m_vegtissue"      = "#7FB8B2", # teal light
  "Bm3235m_spermatangia"   = "#3F9188", # teal mid
  "Bm3235m_spermatia"      = "#0F766E"  # teal dark
)

missing_stage_cols <- setdiff(levels(coldata$stage), names(stage_cols))
if (length(missing_stage_cols) > 0) {
  extra_cols <- grDevices::hcl.colors(length(missing_stage_cols), palette = "Dark 3")
  names(extra_cols) <- missing_stage_cols
  stage_cols <- c(stage_cols, extra_cols)
}
stage_cols <- stage_cols[levels(coldata$stage)]

p <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = stage, shape = generation), size = 7, stroke = 1.4) +
  scale_color_manual(values = stage_cols, drop = FALSE) +
  coord_fixed(ratio = 1) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 8)),
    shape = guide_legend(override.aes = list(size = 8))
  )

pca_pdf <- file.path(out_fig_dir, "qc_pca_stage_generation.pdf")
ggplot2::ggsave(pca_pdf, plot = p, width = 8, height = 8)

pca_out <- pca_df
pca_out$sample_id <- rownames(pca_df)
write_tsv(pca_out, file.path(out_qc_dir, "pca_scores.tsv"))

# --------------------------
# 3) Sample correlation heatmap (rlog)
# --------------------------
rlg_mat <- SummarizedExperiment::assay(rld)
sample_cor <- stats::cor(rlg_mat, method = "pearson", use = "pairwise.complete.obs")

ord <- order(coldata$stage, coldata$replicate, na.last = TRUE)
sample_cor <- sample_cor[ord, ord, drop = FALSE]
col_annot <- data.frame(
  stage = coldata$stage[ord],
  generation = coldata$generation[ord],
  tissue = coldata$tissue[ord],
  row.names = rownames(coldata)[ord],
  stringsAsFactors = FALSE
)

cor_mat_df <- data.frame(sample_id = rownames(sample_cor), sample_cor, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv(cor_mat_df, file.path(out_qc_dir, "sample_correlation_matrix.tsv"))

cor_pdf <- file.path(out_fig_dir, "qc_sample_correlation_heatmap.pdf")
if (requireNamespace("pheatmap", quietly = TRUE)) {
  n_samp <- ncol(sample_cor)
  cell_px <- 10
  dev_size <- max(8, min(18, (n_samp * cell_px) / 72 + 4))
  grDevices::pdf(cor_pdf, width = dev_size, height = dev_size)
  ann_stage_cols <- stage_cols[unique(as.character(col_annot$stage))]
  missing_ann_stage <- setdiff(unique(as.character(col_annot$stage)), names(ann_stage_cols))
  if (length(missing_ann_stage) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_ann_stage), palette = "Dark 3")
    names(extra_cols) <- missing_ann_stage
    ann_stage_cols <- c(ann_stage_cols, extra_cols)
  }
  ann_gen_cols <- c("1n" = "#4D4D4D", "2n" = "#BDBDBD")
  pheatmap::pheatmap(
    sample_cor,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = col_annot,
    annotation_row = col_annot,
    annotation_colors = list(
      stage = ann_stage_cols,
      generation = ann_gen_cols
    ),
    cellwidth = cell_px,
    cellheight = cell_px,
    main = "Sample-Sample Correlation (rlog Pearson)",
    color = grDevices::colorRampPalette(c("#C7A76C", "#F7EFE5", "#5E3C99"))(100),
    border_color = NA
  )
  grDevices::dev.off()
} else {
  grDevices::pdf(cor_pdf, width = 10, height = 10)
  stats::heatmap(
    sample_cor,
    symm = TRUE,
    scale = "none",
    col = grDevices::colorRampPalette(c("#C7A76C", "#F7EFE5", "#5E3C99"))(100),
    margins = c(10, 10),
    main = "Sample-Sample Correlation (rlog Pearson)"
  )
  grDevices::dev.off()
}

cat("Wrote:\n")
cat(" - ", file.path(out_qc_dir, "nfcore_qc_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_qc_dir, "nfcore_qc_by_stage.tsv"), "\n", sep = "")
cat(" - ", file.path(out_qc_dir, "pca_scores.tsv"), "\n", sep = "")
cat(" - ", file.path(out_qc_dir, "sample_correlation_matrix.tsv"), "\n", sep = "")
cat(" - ", pca_pdf, "\n", sep = "")
cat(" - ", cor_pdf, "\n", sep = "")
