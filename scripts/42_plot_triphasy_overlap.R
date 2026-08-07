#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source("R/io.R")
source("R/qc.R")

pairwise_overlap_stats <- function(set_map, universe) {
  nm <- names(set_map)
  out <- data.frame(
    set_a = character(), set_b = character(),
    n_a = integer(), n_b = integer(), n_intersection = integer(),
    jaccard = numeric(), fisher_or = numeric(), fisher_p = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(nm) < 2) return(out)
  U <- length(universe)
  for (i in seq_len(length(nm) - 1)) {
    for (j in seq((i + 1), length(nm))) {
      a <- unique(set_map[[nm[i]]])
      b <- unique(set_map[[nm[j]]])
      ia <- intersect(a, b)
      only_a <- setdiff(a, b)
      only_b <- setdiff(b, a)
      neither <- max(0L, U - length(ia) - length(only_a) - length(only_b))
      m <- matrix(c(length(ia), length(only_a), length(only_b), neither), nrow = 2, byrow = TRUE)
      ft <- tryCatch(fisher.test(m), error = function(e) NULL)
      out <- rbind(out, data.frame(
        set_a = nm[i],
        set_b = nm[j],
        n_a = length(a),
        n_b = length(b),
        n_intersection = length(ia),
        jaccard = if ((length(a) + length(b) - length(ia)) > 0) {
          length(ia) / (length(a) + length(b) - length(ia))
        } else {
          0
        },
        fisher_or = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
        fisher_p = if (!is.null(ft)) as.numeric(ft$p.value) else NA_real_,
        stringsAsFactors = FALSE
      ))
    }
  }
  out
}

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))

out_dir <- file.path(paths$output_root, "generation")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

biased_file <- file.path(out_dir, "triphasy_biased_gene_set_membership.tsv")
repressed_file <- file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv")
if (!file.exists(biased_file) || !file.exists(repressed_file)) {
  stop("Missing triphasy biased/repressed membership files. Run 73_generation_transcriptome_overlap.R first.", call. = FALSE)
}

set_levels <- c(
  "GAM_core_biased", "GAM_core_downregulated",
  "CSP_biased", "CSP_downregulated",
  "TETRA_biased", "TETRA_downregulated"
)
set_labels <- c(
  GAM_core_biased = "GAM biased",
  GAM_core_downregulated = "GAM repressed",
  CSP_biased = "Carposporophyte biased",
  CSP_downregulated = "Carposporophyte repressed",
  TETRA_biased = "Tetrasporophyte biased",
  TETRA_downregulated = "Tetrasporophyte repressed"
)

expr <- prepare_expression_inputs(paths, design, analysis)
universe <- normalize_gene_id(unique(expr$tpm$gene_id))

membership <- bind_rows(
  read.delim(biased_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(direction = "biased"),
  read.delim(repressed_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(direction = "repressed")
) %>%
  mutate(gene_id = normalize_gene_id(gene_id)) %>%
  filter(set %in% set_levels, gene_id %in% universe) %>%
  distinct(gene_id, set, direction)

set_map <- split(membership$gene_id, factor(membership$set, levels = set_levels))
set_map <- set_map[set_levels]
set_map <- set_map[!vapply(set_map, is.null, logical(1))]

stats <- pairwise_overlap_stats(set_map, universe) %>%
  mutate(
    fisher_fdr = p.adjust(fisher_p, method = "BH"),
    log2_fisher_or = log2(fisher_or),
    set_a_label = unname(set_labels[set_a]),
    set_b_label = unname(set_labels[set_b])
  ) %>%
  arrange(fisher_fdr, desc(jaccard))
write_tsv(stats, file.path(out_dir, "triphasy_up_down_pairwise_overlap_stats.tsv"))

heat_df <- stats %>%
  transmute(set_a, set_b, jaccard, fisher_or, fisher_fdr, log2_fisher_or) %>%
  bind_rows(
    stats %>%
      select(set_a_orig = set_a, set_b_orig = set_b, jaccard, fisher_or, fisher_fdr, log2_fisher_or) %>%
      transmute(set_a = set_b_orig, set_b = set_a_orig, jaccard, fisher_or, fisher_fdr, log2_fisher_or)
  ) %>%
  bind_rows(data.frame(
    set_a = set_levels, set_b = set_levels, jaccard = NA_real_,
    fisher_or = NA_real_, fisher_fdr = NA_real_, log2_fisher_or = NA_real_,
    stringsAsFactors = FALSE
  )) %>%
  mutate(
    set_a_label = factor(unname(set_labels[set_a]), levels = unname(set_labels[set_levels])),
    set_b_label = factor(unname(set_labels[set_b]), levels = rev(unname(set_labels[set_levels]))),
    jaccard_bin = cut(
      jaccard,
      breaks = c(-Inf, 0, 0.005, 0.01, 0.03, 0.07, 0.12, Inf),
      labels = c("0", "0-0.005", "0.005-0.01", "0.01-0.03", "0.03-0.07", "0.07-0.12", ">0.12"),
      include.lowest = TRUE,
      right = TRUE
    ),
    jaccard_label = if_else(is.na(jaccard), "", sprintf("%.3f", jaccard)),
    fdr_sig = is.finite(fisher_fdr) & fisher_fdr <= 0.05,
    log2_or_label = if_else(is.na(log2_fisher_or), "", sprintf("%.2f", log2_fisher_or))
  )
write_tsv(heat_df, file.path(out_dir, "triphasy_up_down_overlap_heatmap_values.tsv"))

p_jaccard <- ggplot(heat_df, aes(x = set_a_label, y = set_b_label, fill = jaccard_bin)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = jaccard_label), size = 2.4, color = "black") +
  scale_fill_manual(
    values = c(
      "0" = "#F5F5F5",
      "0-0.005" = "#EFE4D4",
      "0.005-0.01" = "#E7D5BC",
      "0.01-0.03" = "#D3B089",
      "0.03-0.07" = "#BE946C",
      "0.07-0.12" = "#8F3A3A",
      ">0.12" = "#663333"
    ),
    na.value = "#D9D9D9",
    drop = FALSE
  ) +
  labs(
    title = "Overlap Between Triphasy Biased and Repressed Gene Sets",
    subtitle = "Jaccard index; diagonal omitted",
    x = NULL,
    y = NULL,
    fill = "Jaccard"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_overlap_jaccard_heatmap.pdf"), p_jaccard, width = 7.4, height = 6.6)

or_limit <- max(abs(heat_df$log2_fisher_or[is.finite(heat_df$log2_fisher_or)]), na.rm = TRUE)
or_limit <- min(or_limit, 8)
p_or <- ggplot(heat_df, aes(x = set_a_label, y = set_b_label, fill = pmax(pmin(log2_fisher_or, or_limit), -or_limit))) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = log2_or_label, fontface = if_else(fdr_sig, "bold", "plain")), size = 2.4, color = "black") +
  scale_fill_gradientn(
    colors = c("#6B6B6B", "white", "#CC9966", "#663333"),
    values = rescale(c(-or_limit, 0, or_limit * 0.45, or_limit)),
    limits = c(-or_limit, or_limit),
    oob = squish,
    na.value = "#D9D9D9",
    name = "log2 odds\nratio"
  ) +
  labs(
    title = "Enrichment/Depletion Between Triphasy Biased and Repressed Gene Sets",
    subtitle = "Positive = enriched overlap; negative = depleted overlap; bold labels indicate Fisher FDR <= 0.05",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(fig_dir, "triphasy_up_down_overlap_log2_oddsratio_heatmap.pdf"), p_or, width = 7.6, height = 6.6)

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "triphasy_up_down_pairwise_overlap_stats.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_overlap_heatmap_values.tsv"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_overlap_jaccard_heatmap.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_overlap_log2_oddsratio_heatmap.pdf"), "\n", sep = "")
