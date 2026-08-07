#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_tai <- file.path(paths$output_root, "tai")
ensure_dir(out_tai)
fig_w <- 8.8
fig_h <- 5.4

inp_main <- file.path(out_tai, "tai_plot_input.tsv")
inp_raw <- file.path(out_tai, "tai_signature_values.tsv")
inp_trans <- file.path(out_tai, "tai_signature_transformed_values.tsv")
inp_rank <- file.path(out_tai, "tai_signature_rank_transformed_values.tsv")
inp_quant <- file.path(out_tai, "tai_signature_gene_quantiles_values.tsv")
inp_flt <- file.path(out_tai, "tai_flatline_test.tsv")

need <- c(inp_main, inp_raw, inp_trans, inp_rank, inp_quant)
miss <- need[!file.exists(need)]
if (length(miss) > 0) {
  stop(sprintf("Missing required TAI inputs: %s", paste(miss, collapse = ", ")), call. = FALSE)
}

main_df <- read.delim(inp_main, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
raw_df <- read.delim(inp_raw, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
trans_df <- read.delim(inp_trans, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rank_df <- read.delim(inp_rank, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
quant_df <- read.delim(inp_quant, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
flt_df <- if (file.exists(inp_flt)) read.delim(inp_flt, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) else data.frame()

if (!all(c("stage_order", "stage_id", "TAI") %in% colnames(main_df))) {
  stop("tai_plot_input.tsv missing required columns: stage_order, stage_id, TAI", call. = FALSE)
}

main_df <- main_df[order(main_df$stage_order), , drop = FALSE]
stage_levels <- as.character(main_df$stage_id)
main_df$stage_id <- factor(main_df$stage_id, levels = stage_levels)

stage_cols <- c(
  "Bm3235f_vegtissue"       = "#AC5858",
  "Bm3235f_procarps"        = "#78211C",
  "Bm3235m_vegtissue"       = "#8EB8B2",
  "Bm3235m_spermatangia"    = "#5B9388",
  "Bm3235m_spermatia"       = "#387568",
  "Bm3235f_pericarp"        = "#B87B7E",
  "Bm3235f_cystocarps_imm"  = "#ACA4C7",
  "Bm3235f_cystocarps_mat"  = "#8172AA",
  "Bm3235c_carposporophyte" = "#5E3895",
  "Bm3235t_carpospores"     = "#442076",
  "Bm3235t_vegtissue"       = "#C7DA9F",
  "Bm3235t_stichidia"       = "#A2BE72",
  "Bm3235t_tetraspores"     = "#7D9E59"
)
missing_stage_cols <- setdiff(stage_levels, names(stage_cols))
if (length(missing_stage_cols) > 0) {
  fallback <- rep("#9E9E9E", length(missing_stage_cols))
  names(fallback) <- missing_stage_cols
  stage_cols <- c(stage_cols, fallback)
}
stage_cols <- stage_cols[stage_levels]

flt_lab <- ""
if (nrow(flt_df) > 0 && "p_value" %in% colnames(flt_df)) {
  pv <- suppressWarnings(as.numeric(flt_df$p_value[1]))
  if (is.finite(pv)) flt_lab <- sprintf("Flatline p = %.2e", pv)
}

# 1) Main ordered signature omitted; manuscript uses transformed and quantile panels.

# 2) Transformed signatures
# 2) Transformed signatures
trans_df$stage <- sub("^mean_tpm_", "", as.character(trans_df$stage))
trans_df <- trans_df[trans_df$stage %in% stage_levels, , drop = FALSE]
trans_df$stage <- factor(trans_df$stage, levels = stage_levels)
t_levels <- unique(as.character(trans_df$transformation))
t_levels <- t_levels[order(t_levels)]
if ("none" %in% t_levels) {
  t_levels <- c("none", setdiff(t_levels, "none"))
}
trans_df$transformation <- factor(trans_df$transformation, levels = t_levels)
t_others <- setdiff(t_levels, "none")
t_alpha <- setNames(rep(1, length(t_levels)), t_levels)
t_lt <- setNames(rep("solid", length(t_levels)), t_levels)
t_cols <- setNames(rep("#000000", length(t_levels)), t_levels)
if (length(t_others) > 0) {
  t_alpha[t_others] <- seq(0.35, 0.90, length.out = length(t_others))
  t_cols[t_others] <- grDevices::colorRampPalette(c("#BDBDBD", "#000000"))(length(t_others))
  lt_base <- c("11", "21", "31", "41", "solid")
  t_lt[t_others] <- lt_base[pmax(1, round(seq(1, length(lt_base), length.out = length(t_others))))]
}
if ("none" %in% names(t_cols)) t_cols["none"] <- "#AC5858"
if ("none" %in% names(t_alpha)) t_alpha["none"] <- 1
if ("none" %in% names(t_lt)) t_lt["none"] <- "solid"
trans_df$line_w <- ifelse(as.character(trans_df$transformation) == "none", 1.2, 0.9)
p_trans <- ggplot(trans_df, aes(x = stage, y = value, group = transformation)) +
  geom_line(aes(alpha = transformation, linetype = transformation, color = transformation, linewidth = line_w)) +
  geom_point(aes(fill = stage), shape = 21, size = 3.0, stroke = 0, color = "#00000000") +
  scale_linewidth_identity() +
  scale_color_manual(values = t_cols, drop = FALSE) +
  scale_fill_manual(values = stage_cols, drop = FALSE) +
  scale_alpha_manual(values = t_alpha, drop = FALSE) +
  scale_linetype_manual(values = t_lt, drop = FALSE) +
  labs(
    title = "TAI Signature by Count Transformation",
    x = NULL,
    y = "TAI",
    color = "Transformation",
    fill = "Stage",
    alpha = "Transformation",
    linetype = "Transformation"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(
  file.path(out_tai, "tai_signature_transformed_replot_with_legend.pdf"),
  p_trans + theme(legend.position = "right"),
  width = fig_w,
  height = fig_h
)

# 3) Rank transformed vs raw omitted from publication workflow.

# 4) Gene quantile curves
# 4) Gene quantile curves
quant_df$stage <- sub("^mean_tpm_", "", as.character(quant_df$stage))
quant_df <- quant_df[quant_df$stage %in% stage_levels, , drop = FALSE]
quant_df$stage <- factor(quant_df$stage, levels = stage_levels)
quant_df$q_lab <- paste0("q=", quant_df$quantile)
quant_df$line_w <- ifelse(quant_df$q_lab == "q=1", 1.2, 0.9)

q_levels <- unique(quant_df$q_lab[order(quant_df$quantile)])
quant_df$q_lab <- factor(quant_df$q_lab, levels = q_levels)
# Reversed lightening order: lower quantiles lighter, higher quantiles darker.
q_alpha <- seq(0.25, 1.00, length.out = length(q_levels))
names(q_alpha) <- q_levels
# Harmonious accent for q=1 while keeping other quantiles grayscale.
q_line_cols <- rep("#000000", length(q_levels))
names(q_line_cols) <- q_levels
if ("q=1" %in% names(q_line_cols)) q_line_cols["q=1"] <- "#AC5858"
# Dash density progression: many dashes -> fewer dashes -> solid.
lt_base <- c("11", "21", "31", "41", "solid")
lt_vals <- lt_base[pmax(1, round(seq(1, length(lt_base), length.out = length(q_levels))))]
names(lt_vals) <- q_levels

p_quant <- ggplot(quant_df, aes(x = stage, y = value, group = q_lab)) +
  geom_line(aes(alpha = q_lab, linetype = q_lab, linewidth = line_w, color = q_lab)) +
  geom_point(aes(fill = stage), shape = 21, size = 3.0, stroke = 0, color = "#00000000") +
  labs(
    title = "TAI Signature Across Gene-Quantile Cuts",
    x = NULL,
    y = "TAI",
    color = "Quantile",
    fill = "Stage",
    alpha = "Quantile",
    linetype = "Quantile"
  ) +
  scale_linewidth_identity() +
  scale_alpha_manual(values = q_alpha, drop = FALSE) +
  scale_linetype_manual(values = lt_vals, drop = FALSE) +
  scale_color_manual(values = q_line_cols, drop = FALSE) +
  scale_fill_manual(values = stage_cols, drop = FALSE) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(
  file.path(out_tai, "tai_signature_gene_quantiles_replot_with_legend.pdf"),
  p_quant + theme(legend.position = "right"),
  width = fig_w,
  height = fig_h
)

cat("Wrote:\n")
cat(" - ", file.path(out_tai, "tai_signature_transformed_replot_with_legend.pdf"), "\n", sep = "")
cat(" - ", file.path(out_tai, "tai_signature_gene_quantiles_replot_with_legend.pdf"), "\n", sep = "")
