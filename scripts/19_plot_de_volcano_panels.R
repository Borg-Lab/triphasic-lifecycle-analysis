#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

if (!requireNamespace("ggrastr", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("ggrastr", ask = FALSE, update = FALSE)
}


source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
contrasts <- read_contrast_table(Sys.getenv("CONTRASTS_CONFIG", unset = "config/contrasts.tsv"))

de_path <- file.path(paths$output_root, "de", "de_results.tsv")
if (!file.exists(de_path)) stop(sprintf("Missing DE file: %s", de_path), call. = FALSE)
de <- read.delim(de_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(de) == 0) stop("de_results.tsv is empty.", call. = FALSE)

out_de_dir <- file.path(paths$output_root, "de")
out_fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_de_dir)
ensure_dir(out_fig_dir)

de <- merge(
  de,
  contrasts[, c("contrast_id", "lfc_threshold", "fdr_threshold")],
  by = "contrast_id",
  all.x = TRUE,
  sort = FALSE
)

de$lfc_threshold[!is.finite(de$lfc_threshold)] <- 1
de$fdr_threshold[!is.finite(de$fdr_threshold)] <- 0.05
de$neg_log10_fdr <- -log10(pmax(de$fdr, 1e-300))
de$direction <- "Not significant"
de$direction[is.finite(de$fdr) & is.finite(de$log2fc) & de$fdr <= de$fdr_threshold & de$log2fc <= -de$lfc_threshold] <- "Down"
de$direction[is.finite(de$fdr) & is.finite(de$log2fc) & de$fdr <= de$fdr_threshold & de$log2fc >=  de$lfc_threshold] <- "Up"
de$contrast_id <- factor(de$contrast_id, levels = contrasts$contrast_id)

summary_rows <- lapply(split(de, de$contrast_id, drop = TRUE), function(d) {
  data.frame(
    contrast_id = as.character(d$contrast_id[1]),
    tested_genes = sum(is.finite(d$fdr)),
    up_sig = sum(d$direction == "Up", na.rm = TRUE),
    down_sig = sum(d$direction == "Down", na.rm = TRUE),
    total_sig = sum(d$direction %in% c("Up", "Down"), na.rm = TRUE),
    lfc_threshold = unique(d$lfc_threshold)[1],
    fdr_threshold = unique(d$fdr_threshold)[1],
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[match(contrasts$contrast_id, summary_df$contrast_id), , drop = FALSE]

summary_file <- file.path(out_de_dir, "de_contrast_summary.tsv")
write_tsv(summary_df, summary_file)

# Cap only the drawing range so extreme adjusted p-values do not flatten all panels.
de$neg_log10_fdr_plot <- pmin(de$neg_log10_fdr, 80)
de$log2fc_plot <- pmax(pmin(de$log2fc, 15), -15)

label_df <- summary_df
label_df$contrast_id <- factor(label_df$contrast_id, levels = contrasts$contrast_id)
label_df$x <- -14
label_df$y <- 76
label_df$label <- sprintf("Down: %s\nUp: %s", label_df$down_sig, label_df$up_sig)

volcano_file <- file.path(out_fig_dir, "de_volcanoplots.pdf")
p <- ggplot(de, aes(x = log2fc_plot, y = neg_log10_fdr_plot, color = direction)) +
  ggrastr::geom_point_rast(size = 0.75, alpha = 0.65, stroke = 0, raster.dpi = 600) +
  geom_vline(aes(xintercept = lfc_threshold), linewidth = 0.15, linetype = "dashed", color = "grey35") +
  geom_vline(aes(xintercept = -lfc_threshold), linewidth = 0.15, linetype = "dashed", color = "grey35") +
  geom_hline(aes(yintercept = -log10(fdr_threshold)), linewidth = 0.15, linetype = "dashed", color = "grey35") +
  geom_text(data = label_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.2, color = "grey15") +
  facet_wrap(~ contrast_id, ncol = 5) +
  scale_color_manual(
    values = c("Down" = "#2B7A78", "Not significant" = "#D6D0C4", "Up" = "#663333"),
    breaks = c("Up", "Down", "Not significant")
  ) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0, 80), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(aspect.ratio = 1) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.2),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8)
  ) +
  labs(x = "log2 fold-change (clipped at +/-15)", y = "-log10 FDR (clipped at 80)", color = NULL)

ggsave(volcano_file, p, width = 9.5, height = 8.2, useDingbats = FALSE)

cat("Wrote:\n")
cat(" - ", summary_file, "\n", sep = "")
cat(" - ", volcano_file, "\n", sep = "")
