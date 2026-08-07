#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_root <- paths$output_root
out_tab <- file.path(out_root, "tables")
out_fig <- file.path(out_root, "figures")
ensure_dir(out_tab)
ensure_dir(out_fig)

set.seed(3235)
n_perm <- 2000L
pseudo <- 0.5

age_file <- paths$gene_age_tsv
module_file <- file.path(out_root, "modules", "module_membership_labeled.tsv")
if (!file.exists(module_file)) module_file <- file.path(out_root, "modules", "module_membership.tsv")
if (!file.exists(age_file)) stop("Missing gene_age_tsv: ", age_file, call. = FALSE)
if (!file.exists(module_file)) stop("Missing module membership table: ", module_file, call. = FALSE)

ages <- read.delim(age_file, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ages) < 3) stop("gene_age_tsv must have >=3 columns", call. = FALSE)
colnames(ages)[1:3] <- c("transcript_id", "age_name", "rank")
ages$gene_id <- sub("\\.t[0-9]+$", "", as.character(ages$transcript_id))
ages$rank <- suppressWarnings(as.integer(ages$rank))
ages <- ages[!is.na(ages$gene_id) & nzchar(ages$gene_id) & is.finite(ages$rank), c("gene_id", "rank", "age_name"), drop = FALSE]
ages <- unique(ages)

rank_vals <- sort(unique(ages$rank))
rank_name_map <- tapply(ages$age_name, ages$rank, function(x) names(sort(table(x), decreasing = TRUE))[1])
rank_label_map <- setNames(
  paste0("PS", rank_vals, ": ", unname(rank_name_map[as.character(rank_vals)])),
  rank_vals
)
ages$rank_label <- unname(rank_label_map[as.character(ages$rank)])

mods <- read.delim(module_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rank_cols <- grep("^module_rank($|\\.)", colnames(mods), value = TRUE)
label_cols <- grep("^module_label($|\\.)", colnames(mods), value = TRUE)
module_col <- if (length(rank_cols) > 0) rank_cols[1] else "module_id"
label_col <- if (length(label_cols) > 0) label_cols[1] else NA_character_
if (!is.na(label_col)) {
  mod <- unique(mods[, c("gene_id", module_col, label_col), drop = FALSE])
  colnames(mod) <- c("gene_id", "module_rank", "module_label")
} else {
  mod <- unique(mods[, c("gene_id", module_col), drop = FALSE])
  colnames(mod) <- c("gene_id", "module_rank")
  mod$module_label <- mod$module_rank
}
module_label_lookup <- unique(mod[, c("module_rank", "module_label"), drop = FALSE])
module_label_lookup <- module_label_lookup[!duplicated(module_label_lookup$module_rank), , drop = FALSE]

# Universe matches the sex-bias phylorank analysis: protein-coding genes with
# genomic coordinates, excluding organellar chromosomes. This avoids comparing
# module genes against unplaced/non-coordinate or organellar genes not available
# to chromosomal analyses.
pc_genes <- get_protein_coding_gene_ids(paths)
coords_all <- read_coordinates(paths)
organellar_chrs <- c("ChrMt", "ChrPt")
coord_universe <- unique(coords_all$gene_id[coords_all$gene_id %in% pc_genes & !coords_all$chr %in% organellar_chrs])

ages_u <- unique(ages[, c("gene_id", "rank_label"), drop = FALSE])
ages_u <- ages_u[
  !is.na(ages_u$gene_id) & nzchar(ages_u$gene_id) & !is.na(ages_u$rank_label) & ages_u$gene_id %in% coord_universe,
  ,
  drop = FALSE
]
if (nrow(ages_u) < 10) stop("Too few genes in phylorank universe.", call. = FALSE)

universe_audit <- data.frame(
  metric = c(
    "phylorank_background",
    "protein_coding_genes",
    "protein_coding_genes_with_coordinates_nonorganellar",
    "organellar_chromosomes_excluded",
    "genes_with_phylorank_in_filtered_background"
  ),
  value = c(
    "protein_coding_with_coordinates_excluding_ChrMt_ChrPt",
    length(pc_genes),
    length(coord_universe),
    paste(organellar_chrs, collapse = ","),
    nrow(ages_u)
  ),
  stringsAsFactors = FALSE
)
write_tsv(universe_audit, file.path(out_tab, "module_phylorank_universe_audit.tsv"))

# Module genes intersected with universe
dat <- merge(mod, ages_u, by = "gene_id", all = FALSE)
dat <- unique(dat[, c("gene_id", "module_rank", "module_label", "rank_label"), drop = FALSE])
if (nrow(dat) < 10) stop("Too few module genes with phylorank assignment.", call. = FALSE)

module_levels <- sort(unique(dat$module_rank))
if (all(grepl("^M[0-9]+$", module_levels))) {
  module_levels <- module_levels[order(as.integer(sub("^M", "", module_levels)))]
}
rank_levels <- unname(rank_label_map[as.character(rank_vals)])

dat$module_rank <- factor(dat$module_rank, levels = module_levels)
dat$rank_label <- factor(dat$rank_label, levels = rank_levels)
ages_u$rank_label <- factor(ages_u$rank_label, levels = rank_levels)

obs_tab <- table(dat$module_rank, dat$rank_label) # module x rank observed
obs_mat <- unclass(obs_tab)

N <- nrow(ages_u)
mod_tot <- rowSums(obs_mat)
rank_tot <- as.numeric(table(ages_u$rank_label))
names(rank_tot) <- levels(ages_u$rank_label)
exp_mat <- outer(mod_tot, rank_tot, FUN = function(a, b) a * b / N)
log2oe <- log2((obs_mat + pseudo) / (exp_mat + pseudo))

# Global association test
# Global test against genome background: module x rank observed vs expected
chi_asym <- suppressWarnings(chisq.test(x = obs_mat, p = rank_tot / N, rescale.p = TRUE))
chi_sim <- suppressWarnings(chisq.test(x = obs_mat, p = rank_tot / N, rescale.p = TRUE, simulate.p.value = TRUE, B = 20000))
global_df <- data.frame(
  n_genes_tested = N,
  n_modules = nrow(obs_mat),
  n_phyloranks = ncol(obs_mat),
  chisq_statistic = unname(chi_asym$statistic),
  chisq_df = unname(chi_asym$parameter),
  chisq_pvalue = chi_asym$p.value,
  chisq_simulated_pvalue = chi_sim$p.value,
  stringsAsFactors = FALSE
)

# Per-module x phylorank one-sided Fisher enrichment/depletion
res <- vector("list", length = nrow(obs_mat) * ncol(obs_mat))
k <- 1L
for (i in seq_len(nrow(obs_mat))) {
  for (j in seq_len(ncol(obs_mat))) {
    a <- obs_mat[i, j]
    b <- mod_tot[i] - a
    c <- rank_tot[j] - a
    d <- N - a - b - c
    contingency <- matrix(c(a, b, c, d), nrow = 2)
    ft_enrich <- fisher.test(contingency, alternative = "greater")
    ft_deplete <- fisher.test(contingency, alternative = "less")
    res[[k]] <- data.frame(
      module_rank = rownames(obs_mat)[i],
      rank_label = colnames(obs_mat)[j],
      observed = a,
      expected = exp_mat[i, j],
      log2_obs_exp = log2oe[i, j],
      odds_ratio = unname(ft_enrich$estimate),
      fisher_enrich_p = ft_enrich$p.value,
      fisher_deplete_p = ft_deplete$p.value,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
res_df <- do.call(rbind, res)
res_df$fisher_enrich_fdr <- p.adjust(res_df$fisher_enrich_p, method = "BH")
res_df$fisher_enrich_fdr_within_module <- ave(res_df$fisher_enrich_p, res_df$module_rank, FUN = function(p) {
  p.adjust(p, method = "BH")
})
res_df$fisher_deplete_fdr <- p.adjust(res_df$fisher_deplete_p, method = "BH")
res_df$fisher_deplete_fdr_within_module <- ave(res_df$fisher_deplete_p, res_df$module_rank, FUN = function(p) {
  p.adjust(p, method = "BH")
})

# Permutation p-values (one-sided: enrichment perm log2OE >= observed, depletion perm log2OE <= observed)
ge_counts <- matrix(0L, nrow = nrow(obs_mat), ncol = ncol(obs_mat), dimnames = dimnames(obs_mat))
le_counts <- matrix(0L, nrow = nrow(obs_mat), ncol = ncol(obs_mat), dimnames = dimnames(obs_mat))
# Build module membership indices in full genome universe for permutation.
gene_to_idx <- setNames(seq_len(nrow(ages_u)), ages_u$gene_id)
module_gene_idx <- lapply(module_levels, function(m) {
  gids <- dat$gene_id[dat$module_rank == m]
  as.integer(gene_to_idx[gids])
})
names(module_gene_idx) <- module_levels

ranks_all <- as.character(ages_u$rank_label)

for (b in seq_len(n_perm)) {
  perm_ranks <- sample(ranks_all, replace = FALSE)
  perm_mat <- matrix(0L, nrow = length(module_levels), ncol = length(rank_levels),
                     dimnames = list(module_levels, rank_levels))
  for (i in seq_along(module_levels)) {
    idx <- module_gene_idx[[i]]
    if (length(idx) == 0) next
    tab_i <- table(factor(perm_ranks[idx], levels = rank_levels))
    perm_mat[i, ] <- as.integer(tab_i)
  }
  perm_log2oe <- log2((perm_mat + pseudo) / (exp_mat + pseudo))
  ge_counts <- ge_counts + (perm_log2oe >= log2oe)
  le_counts <- le_counts + (perm_log2oe <= log2oe)
}

emp_enrich_p <- (ge_counts + 1) / (n_perm + 1)
emp_deplete_p <- (le_counts + 1) / (n_perm + 1)
emp_df <- data.frame(
  module_rank = rep(rownames(emp_enrich_p), times = ncol(emp_enrich_p)),
  rank_label = rep(colnames(emp_enrich_p), each = nrow(emp_enrich_p)),
  perm_enrich_p = as.vector(emp_enrich_p),
  perm_deplete_p = as.vector(emp_deplete_p),
  stringsAsFactors = FALSE
)
res_df <- merge(res_df, emp_df, by = c("module_rank", "rank_label"), all.x = TRUE, sort = FALSE)
res_df$module_label <- module_label_lookup$module_label[match(res_df$module_rank, module_label_lookup$module_rank)]
res_df$module_label[is.na(res_df$module_label) | !nzchar(res_df$module_label)] <- res_df$module_rank[is.na(res_df$module_label) | !nzchar(res_df$module_label)]
res_df$perm_enrich_fdr <- p.adjust(res_df$perm_enrich_p, method = "BH")
res_df$perm_enrich_fdr_within_module <- ave(res_df$perm_enrich_p, res_df$module_rank, FUN = function(p) {
  p.adjust(p, method = "BH")
})
res_df$perm_deplete_fdr <- p.adjust(res_df$perm_deplete_p, method = "BH")
res_df$perm_deplete_fdr_within_module <- ave(res_df$perm_deplete_p, res_df$module_rank, FUN = function(p) {
  p.adjust(p, method = "BH")
})

# Backward-compatible aliases: these retain the previous enrichment-only column names.
res_df$fisher_p <- res_df$fisher_enrich_p
res_df$fisher_fdr <- res_df$fisher_enrich_fdr
res_df$fisher_fdr_within_module <- res_df$fisher_enrich_fdr_within_module

res_df$significant_enriched_global <- res_df$fisher_enrich_fdr <= 0.05
res_df$significant_enriched_within_module <- res_df$fisher_enrich_fdr_within_module <= 0.05
res_df$significant_depleted_global <- res_df$fisher_deplete_fdr <= 0.05
res_df$significant_depleted_within_module <- res_df$fisher_deplete_fdr_within_module <= 0.05
res_df$perm_significant_enriched_global <- res_df$perm_enrich_fdr <= 0.05
res_df$perm_significant_enriched_within_module <- res_df$perm_enrich_fdr_within_module <= 0.05
res_df$perm_significant_depleted_global <- res_df$perm_deplete_fdr <= 0.05
res_df$perm_significant_depleted_within_module <- res_df$perm_deplete_fdr_within_module <= 0.05

res_df <- res_df[order(res_df$module_rank, res_df$rank_label), , drop = FALSE]

# Summary per module
summary_df <- aggregate(
  cbind(
    n_enriched_fisher_global = as.integer(res_df$significant_enriched_global),
    n_enriched_fisher_within_module = as.integer(res_df$significant_enriched_within_module),
    n_depleted_fisher_global = as.integer(res_df$significant_depleted_global),
    n_depleted_fisher_within_module = as.integer(res_df$significant_depleted_within_module),
    n_enriched_perm_global = as.integer(res_df$perm_significant_enriched_global),
    n_enriched_perm_within_module = as.integer(res_df$perm_significant_enriched_within_module),
    n_depleted_perm_global = as.integer(res_df$perm_significant_depleted_global),
    n_depleted_perm_within_module = as.integer(res_df$perm_significant_depleted_within_module)
  ) ~ module_rank,
  data = res_df,
  FUN = sum
)
top_hits <- do.call(rbind, lapply(split(res_df, res_df$module_rank), function(d) {
  module_id <- unique(d$module_rank)[1]
  enr <- d[order(d$fisher_enrich_fdr_within_module, -d$log2_obs_exp), , drop = FALSE]
  enr <- enr[enr$fisher_enrich_fdr_within_module <= 0.05, , drop = FALSE]
  dep <- d[order(d$fisher_deplete_fdr_within_module, d$log2_obs_exp), , drop = FALSE]
  dep <- dep[dep$fisher_deplete_fdr_within_module <= 0.05, , drop = FALSE]
  data.frame(
    module_rank = module_id,
    top_enriched_ranks = if (nrow(enr) == 0) NA_character_ else paste(head(enr$rank_label, 3), collapse = " | "),
    top_depleted_ranks = if (nrow(dep) == 0) NA_character_ else paste(head(dep$rank_label, 3), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
summary_df$module_label <- module_label_lookup$module_label[match(summary_df$module_rank, module_label_lookup$module_rank)]
summary_df$module_label[is.na(summary_df$module_label) | !nzchar(summary_df$module_label)] <- summary_df$module_rank[is.na(summary_df$module_label) | !nzchar(summary_df$module_label)]
summary_df <- merge(summary_df, top_hits, by = "module_rank", all.x = TRUE, sort = FALSE)
summary_df <- summary_df[order(summary_df$module_rank), , drop = FALSE]

write_tsv(res_df, file.path(out_tab, "module_phylorank_enrichment_stats.tsv"))
write_tsv(summary_df, file.path(out_tab, "module_phylorank_enrichment_summary.tsv"))
write_tsv(global_df, file.path(out_tab, "module_phylorank_global_association.tsv"))

# Heatmap with significance marks
plot_df <- res_df
module_label_levels <- module_label_lookup$module_label[match(module_levels, module_label_lookup$module_rank)]
module_label_levels[is.na(module_label_levels) | !nzchar(module_label_levels)] <- module_levels[is.na(module_label_levels) | !nzchar(module_label_levels)]
plot_df$module_label <- factor(plot_df$module_label, levels = rev(module_label_levels))
plot_df$rank_label <- factor(plot_df$rank_label, levels = rank_levels)
plot_df$sig_mark <- ifelse(
  plot_df$fisher_enrich_fdr_within_module <= 0.05,
  "*",
  ifelse(plot_df$fisher_deplete_fdr_within_module <= 0.05, "-", "")
)

p <- ggplot(plot_df, aes(x = rank_label, y = module_label, fill = log2_obs_exp)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sig_mark), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#d9d9d9",
    mid = "white",
    high = "#663333",
    midpoint = 0,
    name = "log2(obs/exp)"
  ) +
  labs(
    title = "Module Phylorank Enrichment",
    subtitle = "* enriched, - depleted; within-module Fisher BH FDR<=0.05",
    x = "Phylorank",
    y = "Module"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(file.path(out_fig, "module_phylorank_enrichment_heatmap.pdf"), p, width = 8.5, height = 7.5)

cat("Wrote:\n")
cat(" - ", file.path(out_tab, "module_phylorank_global_association.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "module_phylorank_enrichment_stats.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "module_phylorank_enrichment_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_fig, "module_phylorank_enrichment_heatmap.pdf"), "\n", sep = "")
