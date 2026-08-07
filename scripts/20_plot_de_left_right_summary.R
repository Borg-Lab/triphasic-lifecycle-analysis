#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_root <- paths$output_root
in_file <- file.path(out_root, "de", "de_contrast_summary.tsv")
out_dir <- file.path(out_root, "figures")
out_file <- file.path(out_dir, "de_left_right_bar_summary.pdf")

if (!file.exists(in_file)) stop(sprintf("Missing input: %s", in_file), call. = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.delim(in_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(df) == 0) stop("de_contrast_summary.tsv is empty.", call. = FALSE)

df$contrast_id <- as.character(df$contrast_id)
df$up_sig <- as.numeric(df$up_sig)
df$down_sig <- as.numeric(df$down_sig)
df$total_sig <- as.numeric(df$total_sig)

# Order contrasts by total significant genes for readability.
ord <- order(df$total_sig, decreasing = TRUE)
df <- df[ord, , drop = FALSE]
df$contrast_id <- factor(df$contrast_id, levels = rev(df$contrast_id))

plot_df <- rbind(
  data.frame(contrast_id = df$contrast_id, direction = "up", count = df$up_sig, signed_count = df$up_sig, stringsAsFactors = FALSE),
  data.frame(contrast_id = df$contrast_id, direction = "down", count = df$down_sig, signed_count = -df$down_sig, stringsAsFactors = FALSE)
)

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = signed_count, y = contrast_id, fill = direction)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
  ggplot2::scale_fill_manual(values = c(up = "#D73027", down = "#4575B4")) +
  ggplot2::scale_x_continuous(labels = function(x) abs(x)) +
  ggplot2::labs(
    title = "Differential Expression Summary by Contrast",
    subtitle = "Left = downregulated, Right = upregulated (|log2FC| >= 1, FDR <= 0.05)",
    x = "Number of significant DE genes",
    y = "Contrast",
    fill = "Direction"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 8),
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(out_file, p, width = 11, height = 8.5)

cat("Wrote:\n")
cat(" - ", out_file, "\n", sep = "")
