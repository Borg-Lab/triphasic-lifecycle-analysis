#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_dir <- file.path(paths$output_root, "generation")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

set.seed(3235)
n_perm <- 2000L
pseudo <- 0.5

membership_file <- file.path(out_dir, "triphasy_biased_gene_set_membership.tsv")
down_membership_file <- file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv")
age_file <- paths$gene_age_tsv
if (!file.exists(membership_file)) stop("Missing ", membership_file, ". Run 73_generation_transcriptome_overlap.R first.", call. = FALSE)
if (!file.exists(down_membership_file)) stop("Missing ", down_membership_file, ". Run 73_generation_transcriptome_overlap.R first.", call. = FALSE)
if (!file.exists(age_file)) stop("Missing gene_age_tsv: ", age_file, call. = FALSE)

set_levels <- c(
  "GAM_core_biased",
  "GAM_core_downregulated",
  "CSP_biased",
  "CSP_downregulated",
  "TETRA_biased",
  "TETRA_downregulated"
)
set_labels <- c(
  GAM_core_biased = "GAM biased",
  GAM_core_downregulated = "GAM repressed",
  CSP_biased = "Carposporophyte biased",
  CSP_downregulated = "Carposporophyte repressed",
  TETRA_biased = "Tetrasporophyte biased",
  TETRA_downregulated = "Tetrasporophyte repressed"
)

ages <- read.delim(age_file, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ages) < 3) stop("gene_age_tsv must have >=3 columns", call. = FALSE)
colnames(ages)[1:3] <- c("transcript_id", "age_name", "rank")
ages$gene_id <- normalize_gene_id(as.character(ages$transcript_id))
ages$rank <- suppressWarnings(as.integer(ages$rank))
ages <- ages[!is.na(ages$gene_id) & nzchar(ages$gene_id) & is.finite(ages$rank), c("gene_id", "rank", "age_name"), drop = FALSE]
ages <- unique(ages)
pc_genes <- get_protein_coding_gene_ids(paths)
ages <- ages[ages$gene_id %in% pc_genes, , drop = FALSE]

rank_vals <- sort(unique(ages$rank))
rank_name_map <- tapply(ages$age_name, ages$rank, function(x) names(sort(table(x), decreasing = TRUE))[1])
rank_label_map <- setNames(
  paste0("PS", rank_vals, ": ", unname(rank_name_map[as.character(rank_vals)])),
  rank_vals
)
ages$rank_label <- unname(rank_label_map[as.character(ages$rank)])
rank_levels <- unname(rank_label_map[as.character(rank_vals)])

# Main triphasy phylorank analysis: background is protein-coding genes with a phylostratum assignment.
ages_u <- unique(ages[, c("gene_id", "rank", "rank_label"), drop = FALSE])
ages_u <- ages_u[!is.na(ages_u$gene_id) & nzchar(ages_u$gene_id) & is.finite(ages_u$rank) & !is.na(ages_u$rank_label), , drop = FALSE]
if (nrow(ages_u) < 10) stop("Too few genes in phylorank universe.", call. = FALSE)

membership <- rbind(
  read.delim(membership_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)[, c("gene_id", "set"), drop = FALSE],
  read.delim(down_membership_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)[, c("gene_id", "set"), drop = FALSE]
)
membership$gene_id <- normalize_gene_id(membership$gene_id)
membership <- membership[membership$set %in% set_levels, c("gene_id", "set"), drop = FALSE]
membership <- unique(membership)
membership <- membership[membership$gene_id %in% ages_u$gene_id, , drop = FALSE]
if (nrow(membership) < 10) stop("Too few triphasy up/down genes with phylorank assignment.", call. = FALSE)

# Triphasy up/down sets are tested against the same gene-age background.
dat <- merge(membership, ages_u, by = "gene_id", all = FALSE)
dat <- unique(dat[, c("gene_id", "set", "rank", "rank_label"), drop = FALSE])
set_levels <- set_levels[set_levels %in% unique(dat$set)]

dat$set <- factor(dat$set, levels = set_levels)
dat$rank_label <- factor(dat$rank_label, levels = rank_levels)
ages_u$rank_label <- factor(ages_u$rank_label, levels = rank_levels)

obs_tab <- table(dat$set, dat$rank_label)
obs_mat <- unclass(obs_tab)

N <- nrow(ages_u)
set_tot <- rowSums(obs_mat)
rank_tot <- as.numeric(table(ages_u$rank_label))
names(rank_tot) <- levels(ages_u$rank_label)
exp_mat <- outer(set_tot, rank_tot, FUN = function(a, b) a * b / N)
log2oe <- log2((obs_mat + pseudo) / (exp_mat + pseudo))

chi_asym <- suppressWarnings(chisq.test(x = obs_mat, p = rank_tot / N, rescale.p = TRUE))
chi_sim <- suppressWarnings(chisq.test(x = obs_mat, p = rank_tot / N, rescale.p = TRUE, simulate.p.value = TRUE, B = 20000))
global_df <- data.frame(
  n_genes_tested_background = N,
  n_triphasy_up_down_genes_tested = sum(set_tot),
  n_sets = nrow(obs_mat),
  n_phyloranks = ncol(obs_mat),
  chisq_statistic = unname(chi_asym$statistic),
  chisq_df = unname(chi_asym$parameter),
  chisq_pvalue = chi_asym$p.value,
  chisq_simulated_pvalue = chi_sim$p.value,
  stringsAsFactors = FALSE
)

res <- vector("list", length = nrow(obs_mat) * ncol(obs_mat))
k <- 1L
for (i in seq_len(nrow(obs_mat))) {
  for (j in seq_len(ncol(obs_mat))) {
    a <- obs_mat[i, j]
    b <- set_tot[i] - a
    c <- rank_tot[j] - a
    d <- N - a - b - c
    ft_enrich <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    ft_deplete <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "less")
    res[[k]] <- data.frame(
      set = rownames(obs_mat)[i],
      set_label = unname(set_labels[rownames(obs_mat)[i]]),
      rank_label = colnames(obs_mat)[j],
      observed = as.integer(a),
      expected = as.numeric(exp_mat[i, j]),
      background_total_rank = as.integer(rank_tot[j]),
      set_total = as.integer(set_tot[i]),
      background_total = as.integer(N),
      set_proportion = as.numeric(a / set_tot[i]),
      background_proportion = as.numeric(rank_tot[j] / N),
      log2_obs_exp = as.numeric(log2oe[i, j]),
      odds_ratio = unname(ft_enrich$estimate),
      fisher_enrich_p = ft_enrich$p.value,
      fisher_deplete_p = ft_deplete$p.value,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
res_df <- do.call(rbind, res)
res_df$fisher_enrich_fdr_global <- p.adjust(res_df$fisher_enrich_p, method = "BH")
res_df$fisher_deplete_fdr_global <- p.adjust(res_df$fisher_deplete_p, method = "BH")
res_df$fisher_enrich_fdr_within_set <- ave(res_df$fisher_enrich_p, res_df$set, FUN = function(p) p.adjust(p, method = "BH"))
res_df$fisher_deplete_fdr_within_set <- ave(res_df$fisher_deplete_p, res_df$set, FUN = function(p) p.adjust(p, method = "BH"))

ge_counts <- matrix(0L, nrow = nrow(obs_mat), ncol = ncol(obs_mat), dimnames = dimnames(obs_mat))
le_counts <- matrix(0L, nrow = nrow(obs_mat), ncol = ncol(obs_mat), dimnames = dimnames(obs_mat))
gene_to_idx <- setNames(seq_len(nrow(ages_u)), ages_u$gene_id)
set_gene_idx <- lapply(set_levels, function(s) as.integer(gene_to_idx[dat$gene_id[dat$set == s]]))
names(set_gene_idx) <- set_levels
ranks_all <- as.character(ages_u$rank_label)

for (b in seq_len(n_perm)) {
  perm_ranks <- sample(ranks_all, replace = FALSE)
  perm_mat <- matrix(0L, nrow = length(set_levels), ncol = length(rank_levels), dimnames = list(set_levels, rank_levels))
  for (i in seq_along(set_levels)) {
    idx <- set_gene_idx[[i]]
    if (length(idx) == 0) next
    perm_mat[i, ] <- as.integer(table(factor(perm_ranks[idx], levels = rank_levels)))
  }
  perm_log2oe <- log2((perm_mat + pseudo) / (exp_mat + pseudo))
  ge_counts <- ge_counts + (perm_log2oe >= log2oe)
  le_counts <- le_counts + (perm_log2oe <= log2oe)
}

emp_enrich_p <- (ge_counts + 1) / (n_perm + 1)
emp_deplete_p <- (le_counts + 1) / (n_perm + 1)
emp_df <- data.frame(
  set = rep(rownames(emp_enrich_p), times = ncol(emp_enrich_p)),
  rank_label = rep(colnames(emp_enrich_p), each = nrow(emp_enrich_p)),
  perm_enrich_p = as.vector(emp_enrich_p),
  perm_deplete_p = as.vector(emp_deplete_p),
  stringsAsFactors = FALSE
)
res_df <- merge(res_df, emp_df, by = c("set", "rank_label"), all.x = TRUE, sort = FALSE)
res_df$perm_enrich_fdr_global <- p.adjust(res_df$perm_enrich_p, method = "BH")
res_df$perm_deplete_fdr_global <- p.adjust(res_df$perm_deplete_p, method = "BH")
res_df$perm_enrich_fdr_within_set <- ave(res_df$perm_enrich_p, res_df$set, FUN = function(p) p.adjust(p, method = "BH"))
res_df$perm_deplete_fdr_within_set <- ave(res_df$perm_deplete_p, res_df$set, FUN = function(p) p.adjust(p, method = "BH"))

res_df$significant_enriched_global <- res_df$fisher_enrich_fdr_global <= 0.05
res_df$significant_depleted_global <- res_df$fisher_deplete_fdr_global <= 0.05
res_df$significant_enriched_within_set <- res_df$fisher_enrich_fdr_within_set <= 0.05
res_df$significant_depleted_within_set <- res_df$fisher_deplete_fdr_within_set <= 0.05
res_df$perm_significant_enriched_global <- res_df$perm_enrich_fdr_global <= 0.05
res_df$perm_significant_depleted_global <- res_df$perm_deplete_fdr_global <= 0.05
res_df$perm_significant_enriched_within_set <- res_df$perm_enrich_fdr_within_set <= 0.05
res_df$perm_significant_depleted_within_set <- res_df$perm_deplete_fdr_within_set <= 0.05
res_df <- res_df[order(factor(res_df$set, levels = set_levels), factor(res_df$rank_label, levels = rank_levels)), , drop = FALSE]

summary_df <- aggregate(
  cbind(
    n_enriched_fisher_global = as.integer(res_df$significant_enriched_global),
    n_enriched_fisher_within_set = as.integer(res_df$significant_enriched_within_set),
    n_depleted_fisher_global = as.integer(res_df$significant_depleted_global),
    n_depleted_fisher_within_set = as.integer(res_df$significant_depleted_within_set),
    n_enriched_perm_global = as.integer(res_df$perm_significant_enriched_global),
    n_enriched_perm_within_set = as.integer(res_df$perm_significant_enriched_within_set),
    n_depleted_perm_global = as.integer(res_df$perm_significant_depleted_global),
    n_depleted_perm_within_set = as.integer(res_df$perm_significant_depleted_within_set)
  ) ~ set + set_label,
  data = res_df,
  FUN = sum
)

top_hits <- do.call(rbind, lapply(split(res_df, res_df$set), function(d) {
  set_id <- unique(d$set)[1]
  set_label <- unique(d$set_label)[1]
  enr <- d[d$fisher_enrich_fdr_within_set <= 0.05, , drop = FALSE]
  enr <- enr[order(enr$fisher_enrich_fdr_within_set, -enr$log2_obs_exp), , drop = FALSE]
  dep <- d[d$fisher_deplete_fdr_within_set <= 0.05, , drop = FALSE]
  dep <- dep[order(dep$fisher_deplete_fdr_within_set, dep$log2_obs_exp), , drop = FALSE]
  data.frame(
    set = set_id,
    set_label = set_label,
    top_enriched_ranks = if (nrow(enr) == 0) NA_character_ else paste(head(enr$rank_label, 3), collapse = " | "),
    top_depleted_ranks = if (nrow(dep) == 0) NA_character_ else paste(head(dep$rank_label, 3), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
summary_df <- merge(summary_df, top_hits, by = c("set", "set_label"), all.x = TRUE, sort = FALSE)
summary_df <- summary_df[order(factor(summary_df$set, levels = set_levels)), , drop = FALSE]

category_audit <- data.frame(
  set = set_levels,
  set_label = unname(set_labels[set_levels]),
  n_genes_with_phylorank = as.integer(set_tot[set_levels]),
  background = "Genome_average",
  n_background_genes_with_phylorank = N,
  membership_source = file.path(out_dir, "triphasy_biased_gene_set_membership.tsv"),
  downregulated_membership_source = file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv"),
  matches_plot = file.path(fig_dir, "triphasy_up_down_phylorank_stacked_bar_relative.pdf"),
  stringsAsFactors = FALSE
)

write_tsv(category_audit, file.path(out_dir, "triphasy_up_down_phylorank_category_audit.tsv"))
write_tsv(global_df, file.path(out_dir, "triphasy_up_down_phylorank_global_association.tsv"))
write_tsv(res_df, file.path(out_dir, "triphasy_up_down_phylorank_enrichment_stats.tsv"))
write_tsv(summary_df, file.path(out_dir, "triphasy_up_down_phylorank_enrichment_summary.tsv"))

box_values <- rbind(
  data.frame(
    gene_id = dat$gene_id,
    set = as.character(dat$set),
    set_label = unname(set_labels[as.character(dat$set)]),
    rank = dat$rank,
    rank_label = as.character(dat$rank_label),
    stringsAsFactors = FALSE
  ),
  data.frame(
    gene_id = ages_u$gene_id,
    set = "Genome_average",
    set_label = "Genome average",
    rank = ages_u$rank,
    rank_label = as.character(ages_u$rank_label),
    stringsAsFactors = FALSE
  )
)
box_set_levels <- c(set_levels, "Genome_average")
box_label_levels <- c(unname(set_labels[set_levels]), "Genome average")
box_values$set <- factor(box_values$set, levels = box_set_levels)
box_values$set_label <- factor(box_values$set_label, levels = rev(box_label_levels))
box_values$rank_label <- factor(box_values$rank_label, levels = rank_levels)
write_tsv(box_values, file.path(out_dir, "triphasy_up_down_phylorank_boxplot_values.tsv"))

plot_df <- res_df
plot_df$set_label <- factor(plot_df$set_label, levels = rev(unname(set_labels[set_levels])))
plot_df$rank_label <- factor(plot_df$rank_label, levels = rank_levels)
plot_df$sig_mark <- ifelse(plot_df$fisher_enrich_fdr_within_set <= 0.05, "*", ifelse(plot_df$fisher_deplete_fdr_within_set <= 0.05, "-", ""))

p <- ggplot(plot_df, aes(x = rank_label, y = set_label, fill = log2_obs_exp)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig_mark), size = 3.6, color = "black") +
  scale_fill_gradient2(
    low = "#bdbdbd",
    mid = "white",
    high = "#663333",
    midpoint = 0,
    name = "log2(obs/exp)"
  ) +
  labs(
    title = "Triphasy Biased and Repressed Gene Set Phylorank Enrichment",
    subtitle = "Same categories as triphasy_up_down_phylorank_stacked_bar_relative; * enriched, - depleted; within-set Fisher BH FDR<=0.05",
    x = "Phylorank",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(file.path(fig_dir, "triphasy_up_down_phylorank_enrichment_heatmap.pdf"), p, width = 9.2, height = 5.2)

p_box <- ggplot(box_values, aes(x = set_label, y = rank)) +
  geom_boxplot(fill = "#D9D2C8", color = "#333333", width = 0.62, outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "#663333", color = "white", stroke = 0.25) +
  coord_flip() +
  scale_y_continuous(
    breaks = sort(unique(ages_u$rank)),
    labels = rank_levels[match(sort(unique(ages_u$rank)), rank_vals)]
  ) +
  labs(
    title = "Triphasy Biased and Repressed Gene Sets: Phylorank Age Distribution",
    x = NULL,
    y = "Phylorank"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_phylorank_boxplot.pdf"), p_box, width = 9.2, height = 5.4)

cat("Wrote:
")
cat(" - ", file.path(out_dir, "triphasy_up_down_phylorank_category_audit.tsv"), "
", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_phylorank_global_association.tsv"), "
", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_phylorank_enrichment_stats.tsv"), "
", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_phylorank_enrichment_summary.tsv"), "
", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_phylorank_boxplot_values.tsv"), "
", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_phylorank_enrichment_heatmap.pdf"), "
", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_phylorank_boxplot.pdf"), "
", sep = "")
