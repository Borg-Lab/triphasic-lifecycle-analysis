#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))

in_file <- file.path(paths$output_root, "tables", "group_tpm_means.tsv")
out_plot <- file.path(paths$output_root, "figures", "expr_detectability_barplot.pdf")
out_tsv <- file.path(paths$output_root, "tables", "expr_detectability_summary.tsv")

if (!file.exists(in_file)) stop(sprintf("Missing input: %s", in_file), call. = FALSE)

ensure_dir(dirname(out_plot))
ensure_dir(dirname(out_tsv))

dat <- read.delim(in_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(dat) == 0) stop("group_tpm_means.tsv is empty.", call. = FALSE)

stage_cols <- grep("^mean_tpm_", colnames(dat), value = TRUE)
if (length(stage_cols) == 0) stop("No mean_tpm_ columns found.", call. = FALSE)

keep <- vapply(stage_cols, function(cc) any(is.finite(as.numeric(dat[[cc]]))), logical(1))
stage_cols <- stage_cols[keep]
if (length(stage_cols) == 0) stop("No usable stage columns after NA filtering.", call. = FALSE)

stage_order <- c(
  "Bm3235f_vegtissue", "Bm3235f_procarps",
  "Bm3235m_vegtissue", "Bm3235m_spermatangia", "Bm3235m_spermatia",
  "Bm3235f_pericarp", "Bm3235f_cystocarps_imm", "Bm3235f_cystocarps_mat",
  "Bm3235c_carposporophyte", "Bm3235t_carpospores",
  "Bm3235t_vegtissue", "Bm3235t_stichidia", "Bm3235t_tetraspores"
)
pref <- paste0("mean_tpm_", stage_order)
stage_cols <- c(intersect(pref, stage_cols), setdiff(stage_cols, pref))
stages <- gsub("^mean_tpm_", "", stage_cols)

stage_palette <- c(
  Bm3235f_vegtissue = "#C23B4E",
  Bm3235f_procarps = "#8E2F3F",
  Bm3235m_vegtissue = "#2B7A78",
  Bm3235m_spermatangia = "#1F5F61",
  Bm3235m_spermatia = "#143F4A",
  Bm3235f_pericarp = "#C77C7C",
  Bm3235f_cystocarps_imm = "#7E5AA7",
  Bm3235f_cystocarps_mat = "#5E3C99",
  Bm3235c_carposporophyte = "#3F1F78",
  Bm3235t_carpospores = "#D95F02",
  Bm3235t_vegtissue = "#B8860B",
  Bm3235t_stichidia = "#E6A23C",
  Bm3235t_tetraspores = "#F0C36A"
)

mat <- as.matrix(dat[, stage_cols, drop = FALSE])
storage.mode(mat) <- "numeric"

tpm1 <- colSums(mat > 1, na.rm = TRUE)
tpm3 <- colSums(mat > 3, na.rm = TRUE)
tpm5 <- colSums(mat > 5, na.rm = TRUE)

sum_df <- data.frame(
  stage = stages,
  n_genes_tpm_gt_1 = as.integer(tpm1[stage_cols]),
  n_genes_tpm_gt_3 = as.integer(tpm3[stage_cols]),
  n_genes_tpm_gt_5 = as.integer(tpm5[stage_cols]),
  stringsAsFactors = FALSE
)
write_tsv(sum_df, out_tsv)

plot_df <- rbind(
  data.frame(stage = stages, threshold = "TPM > 1", n_genes = as.integer(tpm1[stage_cols]), stringsAsFactors = FALSE),
  data.frame(stage = stages, threshold = "TPM > 3", n_genes = as.integer(tpm3[stage_cols]), stringsAsFactors = FALSE),
  data.frame(stage = stages, threshold = "TPM > 5", n_genes = as.integer(tpm5[stage_cols]), stringsAsFactors = FALSE)
)
plot_df$stage <- factor(plot_df$stage, levels = stages)
plot_df$threshold <- factor(plot_df$threshold, levels = c("TPM > 1", "TPM > 3", "TPM > 5"))
plot_df$n_genes_k <- plot_df$n_genes / 1000

p <- ggplot(plot_df, aes(x = stage, y = n_genes_k, fill = stage, alpha = threshold)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, color = NA) +
  scale_fill_manual(values = stage_palette, drop = FALSE) +
  scale_alpha_manual(values = c("TPM > 1" = 0.95, "TPM > 3" = 0.62, "TPM > 5" = 0.34), drop = FALSE) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(fill = "none") +
  labs(x = NULL, y = "Detected genes (x1000)")

ggsave(out_plot, p, width = 8.6, height = 5.0, useDingbats = FALSE)

cat("Wrote:\n")
cat(" - ", out_plot, "\n", sep = "")
cat(" - ", out_tsv, "\n", sep = "")
