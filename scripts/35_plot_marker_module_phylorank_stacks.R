#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
out_root <- paths$output_root
out_fig <- file.path(out_root, "figures")
out_tab <- file.path(out_root, "tables")
ensure_dir(out_fig)
ensure_dir(out_tab)

age_file <- paths$gene_age_tsv
marker_file <- file.path(out_root, "markers", "marker_table.tsv")
module_file <- file.path(out_root, "modules", "module_membership_labeled.tsv")
if (!file.exists(module_file)) module_file <- file.path(out_root, "modules", "module_membership.tsv")

if (!file.exists(age_file)) stop("Missing gene_age_tsv: ", age_file, call. = FALSE)
if (!file.exists(marker_file)) stop("Missing marker table: ", marker_file, call. = FALSE)
if (!file.exists(module_file)) stop("Missing module membership table: ", module_file, call. = FALSE)

# ---- Gene age / phylorank ----
ages <- read.delim(age_file, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ages) < 3) stop("gene_age_tsv must have at least 3 columns", call. = FALSE)
colnames(ages)[1:3] <- c("transcript_id", "age_name", "rank")
ages$gene_id <- sub("\\.t[0-9]+$", "", as.character(ages$transcript_id))
ages$rank <- suppressWarnings(as.integer(ages$rank))
ages <- ages[!is.na(ages$gene_id) & nzchar(ages$gene_id) & is.finite(ages$rank), c("gene_id", "rank", "age_name"), drop = FALSE]
ages <- unique(ages)
rank_vals <- sort(unique(ages$rank))
rank_name_map <- tapply(ages$age_name, ages$rank, function(x) names(sort(table(x), decreasing = TRUE))[1])
rank_label_map <- setNames(paste0("PS", rank_vals, ": ", unname(rank_name_map[as.character(rank_vals)])), rank_vals)
ages$rank_label <- unname(rank_label_map[as.character(ages$rank)])
rank_levels <- unname(rank_label_map[as.character(rank_vals)])
# True grouped qualitative mapping by phylorank blocks: 1-2, 3-5, 6-8
rank_nums <- suppressWarnings(as.integer(sub("^PS([0-9]+):.*$", "\\1", rank_levels)))
pal_rank <- setNames(rep("#999999", length(rank_levels)), rank_levels)
# 1-2: light cream shades (user-specified)
pal_rank[rank_nums == 1] <- "#f7efe5"
pal_rank[rank_nums == 2] <- "#e7d5bc"
# 3-5: teal shades
pal_rank[rank_nums == 3] <- "#8CC8C1"
pal_rank[rank_nums == 4] <- "#4FA79D"
pal_rank[rank_nums == 5] <- "#2F7F78"
# 6-8: maroon shades
pal_rank[rank_nums == 6] <- "#B35E5E"
pal_rank[rank_nums == 7] <- "#8F3A3A"
pal_rank[rank_nums == 8] <- "#662424"

# ---- Genome reference ----
genome <- as.data.frame(table(ages$rank_label), stringsAsFactors = FALSE)
colnames(genome) <- c("rank_label", "value")
genome$set_id <- "Genome_average"

# ---- Marker sets ----
markers <- read.delim(marker_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
flag_col <- if ("marker_flag" %in% colnames(markers)) "marker_flag" else if ("strict_marker_flag" %in% colnames(markers)) "strict_marker_flag" else NA_character_
if (is.na(flag_col)) stop("No marker flag column found in marker_table.tsv", call. = FALSE)
mark <- markers[as.logical(markers[[flag_col]]) & !is.na(markers$marker_group) & nzchar(markers$marker_group), c("gene_id", "marker_group"), drop = FALSE]
mark <- unique(mark)
mark <- merge(mark, ages, by = "gene_id", all.x = FALSE, all.y = FALSE)
mark_tab <- as.data.frame(table(mark$marker_group, mark$rank_label), stringsAsFactors = FALSE)
colnames(mark_tab) <- c("set_id", "rank_label", "value")
mark_tab <- mark_tab[mark_tab$value > 0, , drop = FALSE]
mark_tab <- rbind(mark_tab, genome)
mark_tab$set_id <- as.character(mark_tab$set_id)
mark_tab$rank_label <- factor(mark_tab$rank_label, levels = rank_levels)
stage_order_pref <- c(
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
marker_sets <- unique(setdiff(mark_tab$set_id, "Genome_average"))
mark_levels <- c(intersect(stage_order_pref, marker_sets), setdiff(marker_sets, stage_order_pref), "Genome_average")
mark_levels <- rev(mark_levels)
mark_tab$set_id <- factor(mark_tab$set_id, levels = mark_levels)

write_tsv(mark_tab, file.path(out_tab, "marker_phylorank_distribution.tsv"))

p_marker <- ggplot(mark_tab, aes(fill = rank_label, y = value, x = set_id)) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
  scale_fill_manual(values = pal_rank, breaks = rank_levels, limits = rank_levels, drop = FALSE) +
  coord_flip() +
  labs(
    title = "Marker Genes: Phylorank Composition",
    x = "Marker group",
    y = "Relative composition",
    fill = "Phylorank"
  ) +
  theme_bw(base_size = 11)

# Omitted non-manuscript PDF export: ggsave(file.path(out_fig, "marker_phylorank_stacked_bar_relative.pdf"), p_marker, width = 9.5, height = 6)

# ---- Module sets ----
mods <- read.delim(module_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# Prefer full renamed module labels when available; fall back to module rank/id.
label_cols <- grep("^module_label($|\\.)", colnames(mods), value = TRUE)
rank_cols <- grep("^module_rank($|\\.)", colnames(mods), value = TRUE)
module_col <- if (length(label_cols) > 0) label_cols[1] else if (length(rank_cols) > 0) rank_cols[1] else "module_id"
mod <- unique(mods[, c("gene_id", module_col), drop = FALSE])
colnames(mod) <- c("gene_id", "set_id")
mod <- merge(mod, ages, by = "gene_id", all.x = FALSE, all.y = FALSE)
mod_tab <- as.data.frame(table(mod$set_id, mod$rank_label), stringsAsFactors = FALSE)
colnames(mod_tab) <- c("set_id", "rank_label", "value")
mod_tab <- mod_tab[mod_tab$value > 0, , drop = FALSE]
mod_tab <- rbind(mod_tab, genome)
mod_tab$set_id <- as.character(mod_tab$set_id)
mod_tab$rank_label <- factor(mod_tab$rank_label, levels = rank_levels)

mids <- unique(setdiff(mod_tab$set_id, "Genome_average"))
mid_rank <- suppressWarnings(as.integer(sub("^M([0-9]+).*$", "\\1", mids)))
if (all(is.finite(mid_rank))) {
  mids <- mids[order(mid_rank)]
} else {
  mids <- sort(mids)
}
mod_tab$set_id <- factor(mod_tab$set_id, levels = c(mids, "Genome_average"))
mod_tab$set_id <- factor(mod_tab$set_id, levels = rev(c(mids, "Genome_average")))

write_tsv(mod_tab, file.path(out_tab, "module_phylorank_distribution.tsv"))

p_module <- ggplot(mod_tab, aes(fill = rank_label, y = value, x = set_id)) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
  scale_fill_manual(values = pal_rank, breaks = rank_levels, limits = rank_levels, drop = FALSE) +
  coord_flip() +
  labs(
    title = "Modules: Phylorank Composition",
    x = "Module",
    y = "Relative composition",
    fill = "Phylorank"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(out_fig, "module_phylorank_stacked_bar_relative.pdf"), p_module, width = 10, height = 7.5)

cat("Wrote:\n")
cat(" - ", file.path(out_tab, "marker_phylorank_distribution.tsv"), "\n", sep = "")
cat(" - ", file.path(out_tab, "module_phylorank_distribution.tsv"), "\n", sep = "")
cat(" - ", file.path(out_fig, "marker_phylorank_stacked_bar_relative.pdf"), "\n", sep = "")
cat(" - ", file.path(out_fig, "module_phylorank_stacked_bar_relative.pdf"), "\n", sep = "")
