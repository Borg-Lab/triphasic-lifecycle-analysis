#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/io.R")
source("R/qc.R")
source("R/de.R")
source("R/enrichment.R")

pairwise_overlap_stats <- function(set_map, universe) {
  nm <- sort(names(set_map))
  out <- data.frame(
    set_a = character(), set_b = character(),
    n_a = integer(), n_b = integer(),
    n_intersection = integer(), jaccard = numeric(),
    fisher_or = numeric(), fisher_p = numeric(),
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
      out <- rbind(
        out,
        data.frame(
          set_a = nm[i], set_b = nm[j],
          n_a = length(a), n_b = length(b),
          n_intersection = length(ia),
          jaccard = if ((length(a) + length(b) - length(ia)) > 0) length(ia) / (length(a) + length(b) - length(ia)) else 0,
          fisher_or = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
          fisher_p = if (!is.null(ft)) as.numeric(ft$p.value) else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  out
}

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))
contrasts <- read_contrast_table(Sys.getenv("CONTRASTS_CONFIG", unset = "config/contrasts.tsv"))

expr <- prepare_expression_inputs(paths, design, analysis)
means <- compute_group_tpm_means(expr$tpm, design, group_col = "stage")
means$gene_id <- normalize_gene_id(means$gene_id)

de_file <- file.path(paths$output_root, "de", "de_results.tsv")
if (!file.exists(de_file)) stop("Missing DE results: ", de_file, call. = FALSE)
de <- read.delim(de_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
de$gene_id <- normalize_gene_id(de$gene_id)
de$log2fc <- suppressWarnings(as.numeric(de$log2fc))
de$fdr <- suppressWarnings(as.numeric(de$fdr))

lfc_thr <- 1
fdr_thr <- 0.05

de_sig <- de %>%
  filter(is.finite(log2fc), is.finite(fdr), fdr <= fdr_thr, abs(log2fc) >= lfc_thr)

# --------- 1) Generation-vs-stage DE classes ----------
stage_gen <- design %>% distinct(stage, generation) %>% rename(stage_id = stage)
gen_map <- setNames(stage_gen$generation, stage_gen$stage_id)

contrast_meta <- contrasts %>%
  mutate(
    numerator = as.character(numerator),
    denominator = as.character(denominator),
    gen_num = unname(gen_map[numerator]),
    gen_den = unname(gen_map[denominator]),
    cross_generation = !is.na(gen_num) & !is.na(gen_den) & gen_num != gen_den,
    within_generation = !is.na(gen_num) & !is.na(gen_den) & gen_num == gen_den
  )

de_sig <- de_sig %>%
  left_join(
    contrast_meta[, c("contrast_id", "numerator", "denominator", "gen_num", "gen_den", "cross_generation", "within_generation")],
    by = "contrast_id"
  )

gene_bias <- de_sig %>%
  group_by(gene_id) %>%
  summarise(
    generation_biased = any(cross_generation, na.rm = TRUE),
    stage_biased = any(within_generation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bias_class = case_when(
      generation_biased & stage_biased ~ "generation+stage_biased",
      generation_biased & !stage_biased ~ "generation_biased_only",
      !generation_biased & stage_biased ~ "stage_biased_only",
      TRUE ~ "unclassified"
    )
  )

stage_to_simple_group <- function(x) {
  dplyr::case_when(
    x %in% c("Bm3235m_vegtissue", "Bm3235m_spermatangia", "Bm3235m_spermatia", "Bm3235f_vegtissue", "Bm3235f_procarps") ~ "GAM",
    x %in% c("Bm3235t_vegtissue", "Bm3235t_stichidia") ~ "TETRA",
    x %in% c("Bm3235c_carposporophyte") ~ "CSP",
    TRUE ~ NA_character_
  )
}

stage_to_core_group <- function(x) {
  dplyr::case_when(
    x %in% c("Bm3235m_vegtissue", "Bm3235m_spermatangia", "Bm3235f_vegtissue", "Bm3235f_procarps") ~ "GAM_core",
    x %in% c("Bm3235t_vegtissue", "Bm3235t_stichidia") ~ "TETRA",
    x %in% c("Bm3235c_carposporophyte") ~ "CSP",
    TRUE ~ NA_character_
  )
}

build_triphasy_sets <- function(grouped_design, group_levels, label_prefix = "") {
  ctr <- do.call(rbind, lapply(group_levels, function(g) {
    others <- setdiff(group_levels, g)
    data.frame(
      contrast_id = paste0(g, "_vs_", others),
      factor = "triphasy_group",
      numerator = g,
      denominator = others,
      lfc_threshold = lfc_thr,
      fdr_threshold = fdr_thr,
      stringsAsFactors = FALSE
    )
  }))
  rownames(ctr) <- NULL

  de_res <- run_differential_expression(expr$counts, grouped_design, ctr)
  de_res$log2fc <- suppressWarnings(as.numeric(de_res$log2fc))
  de_res$fdr <- suppressWarnings(as.numeric(de_res$fdr))
  sig_up <- de_res %>%
    filter(is.finite(log2fc), is.finite(fdr), fdr <= fdr_thr, log2fc >= lfc_thr)
  sig_down <- de_res %>%
    filter(is.finite(log2fc), is.finite(fdr), fdr <= fdr_thr, log2fc <= -lfc_thr)

  up_sets <- lapply(group_levels, function(g) {
    others <- setdiff(group_levels, g)
    Reduce(intersect, lapply(others, function(o) {
      sig_up$gene_id[sig_up$contrast_id == paste0(g, "_vs_", o)]
    }))
  })
  names(up_sets) <- paste0(label_prefix, group_levels, "_biased")

  down_sets <- lapply(group_levels, function(g) {
    others <- setdiff(group_levels, g)
    Reduce(intersect, lapply(others, function(o) {
      sig_down$gene_id[sig_down$contrast_id == paste0(g, "_vs_", o)]
    }))
  })
  names(down_sets) <- paste0(label_prefix, group_levels, "_downregulated")

  list(
    contrasts = ctr,
    de = de_res,
    sets = up_sets,
    down_sets = down_sets,
    membership = bind_rows(lapply(names(up_sets), function(nm) {
      data.frame(gene_id = up_sets[[nm]], set = nm, direction = "upregulated", stringsAsFactors = FALSE)
    })) %>% distinct(gene_id, set),
    down_membership = bind_rows(lapply(names(down_sets), function(nm) {
      data.frame(gene_id = down_sets[[nm]], set = nm, direction = "downregulated", stringsAsFactors = FALSE)
    })) %>% distinct(gene_id, set),
    sizes = bind_rows(lapply(names(up_sets), function(nm) {
      data.frame(set = nm, direction = "upregulated", n_genes = length(unique(up_sets[[nm]])), stringsAsFactors = FALSE)
    })) %>% arrange(desc(n_genes), set),
    down_sizes = bind_rows(lapply(names(down_sets), function(nm) {
      data.frame(set = nm, direction = "downregulated", n_genes = length(unique(down_sets[[nm]])), stringsAsFactors = FALSE)
    })) %>% arrange(desc(n_genes), set)
  )
}

triphasy_design <- design %>%
  mutate(triphasy_group = stage_to_simple_group(stage)) %>%
  filter(!is.na(triphasy_group)) %>%
  mutate(triphasy_group = factor(triphasy_group, levels = c("GAM", "CSP", "TETRA")))

triphasy_core_design <- design %>%
  mutate(triphasy_group = stage_to_core_group(stage)) %>%
  filter(!is.na(triphasy_group)) %>%
  mutate(triphasy_group = factor(triphasy_group, levels = c("GAM_core", "CSP", "TETRA")))

triphasy_all <- build_triphasy_sets(triphasy_design, c("GAM", "CSP", "TETRA"))
triphasy_core <- build_triphasy_sets(triphasy_core_design, c("GAM_core", "CSP", "TETRA"))

biased_membership <- triphasy_core$membership
biased_set_sizes <- triphasy_core$sizes
biased_gene_sets <- triphasy_core$sets
down_membership <- triphasy_core$down_membership
down_set_sizes <- triphasy_core$down_sizes
down_gene_sets <- triphasy_core$down_sets

core_stage_groups <- list(
  GAM_core_biased = c("Bm3235m_vegtissue", "Bm3235m_spermatangia", "Bm3235f_vegtissue", "Bm3235f_procarps"),
  TETRA_biased = c("Bm3235t_vegtissue", "Bm3235t_stichidia"),
  CSP_biased = c("Bm3235c_carposporophyte")
)

mean_tpm_cols <- setNames(
  paste0("mean_tpm_", make.names(unique(unlist(core_stage_groups)))),
  unique(unlist(core_stage_groups))
)

classify_target_specificity <- function(membership_df, group_map, mean_df, tpm_floor) {
  if (nrow(membership_df) == 0) {
    return(data.frame(
      gene_id = character(), set = character(), target_group = character(),
      target_max_tpm = numeric(), opposing_max_tpm = numeric(),
      target_detected = logical(), opposing_detected = logical(),
      specificity_class = character(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(seq_len(nrow(membership_df)), function(i) {
    gid <- membership_df$gene_id[i]
    set_id <- membership_df$set[i]
    target_stages <- group_map[[set_id]]
    other_stages <- setdiff(unique(unlist(group_map)), target_stages)
    target_cols <- intersect(unname(mean_tpm_cols[target_stages]), colnames(mean_df))
    other_cols <- intersect(unname(mean_tpm_cols[other_stages]), colnames(mean_df))
    idx <- match(gid, mean_df$gene_id)
    target_vals <- if (!is.na(idx) && length(target_cols) > 0) as.numeric(mean_df[idx, target_cols, drop = TRUE]) else numeric()
    other_vals <- if (!is.na(idx) && length(other_cols) > 0) as.numeric(mean_df[idx, other_cols, drop = TRUE]) else numeric()
    target_max <- if (length(target_vals) > 0 && any(is.finite(target_vals))) max(target_vals, na.rm = TRUE) else NA_real_
    other_max <- if (length(other_vals) > 0 && any(is.finite(other_vals))) max(other_vals, na.rm = TRUE) else NA_real_
    target_detected <- is.finite(target_max) && target_max >= tpm_floor
    opposing_detected <- is.finite(other_max) && other_max >= tpm_floor
    data.frame(
      gene_id = gid,
      set = set_id,
      target_group = sub("_biased$", "", set_id),
      target_max_tpm = target_max,
      opposing_max_tpm = other_max,
      target_detected = target_detected,
      opposing_detected = opposing_detected,
      specificity_class = dplyr::case_when(
        target_detected & !opposing_detected ~ "target_detectable_only",
        target_detected & opposing_detected ~ "target_enriched_but_detected_elsewhere",
        !target_detected & opposing_detected ~ "opposing_detected_only",
        TRUE ~ "not_detected_at_floor"
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

up_specificity <- classify_target_specificity(
  biased_membership,
  core_stage_groups,
  means,
  as.numeric(analysis$tpm_detect_threshold)
)

up_specificity_summary <- up_specificity %>%
  count(set, specificity_class, name = "n_genes") %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()

classify_downregulated_target_off <- function(membership_df, stage_groups, means_df, threshold) {
  mean_cols <- grep("^mean_tpm_", colnames(means_df), value = TRUE)
  rows <- lapply(names(stage_groups), function(up_set_id) {
    down_set_id <- sub("_biased$", "_downregulated", up_set_id)
    target_cols <- paste0("mean_tpm_", stage_groups[[up_set_id]])
    target_cols <- intersect(target_cols, colnames(means_df))
    opposing_cols <- setdiff(mean_cols, target_cols)
    if (length(target_cols) == 0 || length(opposing_cols) == 0) {
      return(data.frame())
    }
    dat <- membership_df %>%
      filter(set == down_set_id) %>%
      inner_join(means_df[, c("gene_id", mean_cols), drop = FALSE], by = "gene_id")
    if (nrow(dat) == 0) return(data.frame())
    target_max <- apply(dat[, target_cols, drop = FALSE], 1, max, na.rm = TRUE)
    opposing_max <- apply(dat[, opposing_cols, drop = FALSE], 1, max, na.rm = TRUE)
    all_max <- apply(dat[, mean_cols, drop = FALSE], 1, max, na.rm = TRUE)
    data.frame(
      gene_id = dat$gene_id,
      set = down_set_id,
      target_group = sub("_biased$", "", up_set_id),
      target_max_tpm = target_max,
      opposing_max_tpm = opposing_max,
      all_max_tpm = all_max,
      target_generation_off = target_max < threshold,
      target_off_opposing_detected = target_max < threshold & opposing_max >= threshold,
      all_stages_off = all_max < threshold,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

down_off_status <- classify_downregulated_target_off(
  down_membership,
  core_stage_groups,
  means,
  as.numeric(analysis$tpm_detect_threshold)
)

down_off_summary <- down_off_status %>%
  count(set, target_generation_off, target_off_opposing_detected, all_stages_off, name = "n_genes") %>%
  group_by(set) %>%
  mutate(total_genes = sum(n_genes), fraction = n_genes / total_genes) %>%
  ungroup()

down_target_off_membership <- down_off_status %>%
  filter(target_off_opposing_detected) %>%
  transmute(
    gene_id,
    set = paste0(target_group, "_target_off_downregulated")
  ) %>%
  distinct(gene_id, set)

down_target_off_set_sizes <- down_target_off_membership %>%
  count(set, name = "n_genes") %>%
  arrange(desc(n_genes))

direction_sizes <- bind_rows(
  biased_set_sizes,
  down_set_sizes
) %>%
  mutate(
    generation = sub("_biased$|_downregulated$", "", set),
    generation = recode(generation, GAM_core = "GAM_core", TETRA = "TETRA", CSP = "CSP"),
    direction = factor(direction, levels = c("upregulated", "downregulated"))
  )

gam_sensitivity <- data.frame(
  metric = c("GAM_biased_with_spermatia", "GAM_core_biased_without_spermatia", "overlap", "only_with_spermatia", "only_without_spermatia"),
  n_genes = c(
    length(unique(triphasy_all$sets$GAM_biased)),
    length(unique(triphasy_core$sets$GAM_core_biased)),
    length(intersect(triphasy_all$sets$GAM_biased, triphasy_core$sets$GAM_core_biased)),
    length(setdiff(triphasy_all$sets$GAM_biased, triphasy_core$sets$GAM_core_biased)),
    length(setdiff(triphasy_core$sets$GAM_core_biased, triphasy_all$sets$GAM_biased))
  ),
  stringsAsFactors = FALSE
)

# --------- 3) Cystocarp overlap with tetrasporo/spore programmes ----------
cystocarp_contrasts <- c(
  "Bm3235f_cystocarps_imm_vs_Bm3235f_vegtissue",
  "Bm3235f_cystocarps_mat_vs_Bm3235f_vegtissue",
  "Bm3235f_cystocarps_imm_vs_Bm3235f_procarps",
  "Bm3235f_cystocarps_mat_vs_Bm3235f_procarps"
)

cystocarp_upregulated <- unique(de$gene_id[
  de$contrast_id %in% cystocarp_contrasts &
    is.finite(de$log2fc) & is.finite(de$fdr) &
    de$log2fc >= lfc_thr & de$fdr <= fdr_thr
])

# Preserve spore programme as a separate DE-up set; tetraspores and carpospores are haploid spore products, not TETRA body.
spore_contrasts <- c(
  "Bm3235t_tetraspores_vs_Bm3235t_stichidia",
  "Bm3235t_tetraspores_vs_Bm3235t_vegtissue",
  "Bm3235t_carpospores_vs_Bm3235f_cystocarps_mat",
  "Bm3235t_carpospores_vs_Bm3235f_cystocarps_imm"
)
spore_upregulated <- unique(de$gene_id[
  de$contrast_id %in% spore_contrasts &
    is.finite(de$log2fc) & is.finite(de$fdr) &
    de$log2fc >= lfc_thr & de$fdr <= fdr_thr
])

membership_overlap <- bind_rows(
  biased_membership,
  data.frame(gene_id = cystocarp_upregulated, set = "cystocarp_upregulated", stringsAsFactors = FALSE),
  data.frame(gene_id = spore_upregulated, set = "spore_upregulated", stringsAsFactors = FALSE)
) %>% distinct(gene_id, set)

universe <- unique(means$gene_id)
protein_coding_universe <- universe
set_map <- split(membership_overlap$gene_id, membership_overlap$set)
overlap_stats <- pairwise_overlap_stats(set_map, universe) %>%
  mutate(fisher_fdr = p.adjust(fisher_p, method = "BH")) %>%
  arrange(fisher_fdr, desc(jaccard))

focus_pairs <- overlap_stats %>%
  filter((set_a == "cystocarp_upregulated" & set_b %in% c("GAM_core_biased", "TETRA_biased", "CSP_biased", "spore_upregulated")) |
           (set_b == "cystocarp_upregulated" & set_a %in% c("GAM_core_biased", "TETRA_biased", "CSP_biased", "spore_upregulated")))

# --------- 4) Functional enrichment and phylorank distribution ----------
anno <- read_functional_annotation(paths)
anno_keep <- intersect(colnames(anno), c("gene_id", "GOs", "KEGG_ko", "Description", "Preferred_name"))
anno_use <- anno[, anno_keep, drop = FALSE]
background_anno <- merge(
  data.frame(gene_id = protein_coding_universe, stringsAsFactors = FALSE),
  anno_use,
  by = "gene_id", all.x = TRUE, sort = FALSE
)

make_enrichment_input <- function(membership_df) {
  out <- background_anno
  out$marker_group <- membership_df$set[match(out$gene_id, membership_df$gene_id)]
  out$marker_flag <- !is.na(out$marker_group)
  out$marker_group[!out$marker_flag] <- "background"
  out
}

enrich_df <- make_enrichment_input(biased_membership)
enrich <- run_marker_functional_enrichment(
  marker_df = enrich_df,
  term_specs = list(
    list(source = "GO", col = "GOs", split = ","),
    list(source = "KEGG_ko", col = "KEGG_ko", split = ",")
  ),
  min_term_genes = 5L,
  marker_classes = NULL,
  universe_genes = protein_coding_universe
)

if (nrow(enrich) > 0) {
  enrich$global_fdr <- enrich$fdr
  enrich <- do.call(rbind, lapply(split(enrich, enrich$marker_group), function(d) {
    d$fdr <- stats::p.adjust(d$p_value, method = "BH")
    d
  }))
  rownames(enrich) <- NULL
}

enrich_top <- enrich
if (nrow(enrich_top) > 0) {
  enrich_top <- enrich_top[order(enrich_top$fdr, enrich_top$p_value), , drop = FALSE]
  enrich_top <- do.call(rbind, lapply(split(enrich_top, enrich_top$marker_group), function(d) head(d, 15)))
  rownames(enrich_top) <- NULL
}

down_enrich_df <- make_enrichment_input(down_membership)
down_enrich <- run_marker_functional_enrichment(
  marker_df = down_enrich_df,
  term_specs = list(
    list(source = "GO", col = "GOs", split = ","),
    list(source = "KEGG_ko", col = "KEGG_ko", split = ",")
  ),
  min_term_genes = 5L,
  marker_classes = NULL,
  universe_genes = protein_coding_universe
)

if (nrow(down_enrich) > 0) {
  down_enrich$global_fdr <- down_enrich$fdr
  down_enrich <- do.call(rbind, lapply(split(down_enrich, down_enrich$marker_group), function(d) {
    d$fdr <- stats::p.adjust(d$p_value, method = "BH")
    d
  }))
  rownames(down_enrich) <- NULL
}

down_enrich_top <- down_enrich
if (nrow(down_enrich_top) > 0) {
  down_enrich_top <- down_enrich_top[order(down_enrich_top$fdr, down_enrich_top$p_value), , drop = FALSE]
  down_enrich_top <- do.call(rbind, lapply(split(down_enrich_top, down_enrich_top$marker_group), function(d) head(d, 15)))
  rownames(down_enrich_top) <- NULL
}

down_target_off_enrich_df <- make_enrichment_input(down_target_off_membership)
down_target_off_enrich <- run_marker_functional_enrichment(
  marker_df = down_target_off_enrich_df,
  term_specs = list(
    list(source = "GO", col = "GOs", split = ","),
    list(source = "KEGG_ko", col = "KEGG_ko", split = ",")
  ),
  min_term_genes = 5L,
  marker_classes = NULL,
  universe_genes = protein_coding_universe
)

if (nrow(down_target_off_enrich) > 0) {
  down_target_off_enrich$global_fdr <- down_target_off_enrich$fdr
  down_target_off_enrich <- do.call(rbind, lapply(split(down_target_off_enrich, down_target_off_enrich$marker_group), function(d) {
    d$fdr <- stats::p.adjust(d$p_value, method = "BH")
    d
  }))
  rownames(down_target_off_enrich) <- NULL
}

down_target_off_enrich_top <- down_target_off_enrich
if (nrow(down_target_off_enrich_top) > 0) {
  down_target_off_enrich_top <- down_target_off_enrich_top[order(down_target_off_enrich_top$fdr, down_target_off_enrich_top$p_value), , drop = FALSE]
  down_target_off_enrich_top <- do.call(rbind, lapply(split(down_target_off_enrich_top, down_target_off_enrich_top$marker_group), function(d) head(d, 15)))
  rownames(down_target_off_enrich_top) <- NULL
}

age_file <- paths$gene_age_tsv
if (!file.exists(age_file)) stop("Missing gene_age_tsv: ", age_file, call. = FALSE)
ages <- read.delim(age_file, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ages) < 3) stop("gene_age_tsv must have at least 3 columns", call. = FALSE)
colnames(ages)[1:3] <- c("transcript_id", "age_name", "rank")
ages$gene_id <- normalize_gene_id(ages$transcript_id)
ages$rank <- suppressWarnings(as.integer(ages$rank))
ages <- ages[!is.na(ages$gene_id) & nzchar(ages$gene_id) & is.finite(ages$rank), c("gene_id", "rank", "age_name"), drop = FALSE]
ages <- unique(ages)
rank_vals <- sort(unique(ages$rank))
rank_name_map <- tapply(ages$age_name, ages$rank, function(x) names(sort(table(x), decreasing = TRUE))[1])
rank_label_map <- setNames(paste0("PS", rank_vals, ": ", unname(rank_name_map[as.character(rank_vals)])), rank_vals)
ages$rank_label <- unname(rank_label_map[as.character(ages$rank)])
rank_levels <- unname(rank_label_map[as.character(rank_vals)])
rank_nums <- suppressWarnings(as.integer(sub("^PS([0-9]+):.*$", "\\1", rank_levels)))
pal_rank <- setNames(rep("#999999", length(rank_levels)), rank_levels)
pal_rank[grep("^PS1:", rank_levels)] <- "#f7efe5"
pal_rank[grep("^PS2:", rank_levels)] <- "#e7d5bc"
pal_rank[grep("^PS3:", rank_levels)] <- "#8CC8C1"
pal_rank[grep("^PS4:", rank_levels)] <- "#4FA79D"
pal_rank[grep("^PS5:", rank_levels)] <- "#2F7F78"
pal_rank[grep("^PS6:", rank_levels)] <- "#B35E5E"
pal_rank[grep("^PS7:", rank_levels)] <- "#8F3A3A"
pal_rank[grep("^PS8:", rank_levels)] <- "#662424"

genome_rank <- as.data.frame(table(ages$rank_label), stringsAsFactors = FALSE)
colnames(genome_rank) <- c("rank_label", "value")
genome_rank$set <- "Genome_average"

rank_tab <- merge(biased_membership, ages[, c("gene_id", "rank_label"), drop = FALSE], by = "gene_id", all = FALSE)
rank_tab <- as.data.frame(table(rank_tab$set, rank_tab$rank_label), stringsAsFactors = FALSE)
colnames(rank_tab) <- c("set", "rank_label", "value")
rank_tab <- rank_tab[rank_tab$value > 0, , drop = FALSE]
rank_tab <- rbind(rank_tab, genome_rank)
rank_tab$rank_label <- factor(rank_tab$rank_label, levels = rank_levels)
set_order <- c("GAM_core_biased", "TETRA_biased", "CSP_biased", "Genome_average")
rank_tab$set <- factor(rank_tab$set, levels = rev(intersect(set_order, unique(rank_tab$set))))

rank_set_labels <- c(
  GAM_core_biased = "GAM biased",
  GAM_core_downregulated = "GAM repressed",
  CSP_biased = "Carposporophyte biased",
  CSP_downregulated = "Carposporophyte repressed",
  TETRA_biased = "Tetrasporophyte biased",
  TETRA_downregulated = "Tetrasporophyte repressed",
  Genome_average = "Genome average"
)
rank_set_order <- c(
  "GAM_core_biased",
  "GAM_core_downregulated",
  "CSP_biased",
  "CSP_downregulated",
  "TETRA_biased",
  "TETRA_downregulated",
  "Genome_average"
)
rank_membership_up_down <- bind_rows(
  biased_membership %>% mutate(direction = "biased"),
  down_membership %>% mutate(direction = "repressed")
) %>%
  distinct(gene_id, set, direction)
rank_tab_up_down <- merge(rank_membership_up_down, ages[, c("gene_id", "rank_label"), drop = FALSE], by = "gene_id", all = FALSE)
rank_tab_up_down <- as.data.frame(table(rank_tab_up_down$set, rank_tab_up_down$rank_label), stringsAsFactors = FALSE)
colnames(rank_tab_up_down) <- c("set", "rank_label", "value")
rank_tab_up_down <- rank_tab_up_down[rank_tab_up_down$value > 0, , drop = FALSE]
rank_tab_up_down <- rbind(rank_tab_up_down, genome_rank)
rank_tab_up_down$rank_label <- factor(rank_tab_up_down$rank_label, levels = rank_levels)
rank_tab_up_down$set <- factor(rank_tab_up_down$set, levels = rev(intersect(rank_set_order, unique(rank_tab_up_down$set))))
rank_tab_up_down$set_label <- unname(rank_set_labels[as.character(rank_tab_up_down$set)])
rank_tab_up_down$set_label <- factor(rank_tab_up_down$set_label, levels = rev(unname(rank_set_labels[intersect(rank_set_order, unique(as.character(rank_tab_up_down$set)))])))

# --------- outputs ----------
out_dir <- file.path(paths$output_root, "generation")
fig_dir <- file.path(paths$output_root, "figures")
ensure_dir(out_dir)
ensure_dir(fig_dir)

write_tsv(gene_bias, file.path(out_dir, "generation_stage_bias_gene_classification.tsv"))
write_tsv(
  gene_bias %>% count(bias_class, name = "n_genes") %>% arrange(desc(n_genes)),
  file.path(out_dir, "generation_stage_bias_class_summary.tsv")
)

write_tsv(triphasy_core$de, file.path(out_dir, "triphasy_core_group_de_results.tsv"))
write_tsv(triphasy_all$de, file.path(out_dir, "triphasy_with_spermatia_group_de_results.tsv"))
write_tsv(biased_membership, file.path(out_dir, "triphasy_biased_gene_set_membership.tsv"))
write_tsv(biased_set_sizes, file.path(out_dir, "triphasy_biased_gene_set_sizes.tsv"))
write_tsv(down_membership, file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv"))
write_tsv(down_set_sizes, file.path(out_dir, "triphasy_downregulated_gene_set_sizes.tsv"))
write_tsv(direction_sizes, file.path(out_dir, "triphasy_up_down_gene_set_sizes.tsv"))
write_tsv(up_specificity, file.path(out_dir, "triphasy_biased_target_specificity.tsv"))
write_tsv(up_specificity_summary, file.path(out_dir, "triphasy_biased_target_specificity_summary.tsv"))
write_tsv(down_off_status, file.path(out_dir, "triphasy_downregulated_target_off_status.tsv"))
write_tsv(down_off_summary, file.path(out_dir, "triphasy_downregulated_target_off_summary.tsv"))
write_tsv(down_target_off_membership, file.path(out_dir, "triphasy_downregulated_target_off_gene_set_membership.tsv"))
write_tsv(down_target_off_set_sizes, file.path(out_dir, "triphasy_downregulated_target_off_gene_set_sizes.tsv"))
strict_biased_union <- unique(biased_membership$gene_id)
strict_repressed_union <- unique(down_membership$gene_id)
strict_up_down_status_sizes <- data.frame(
  set = c("biased_only", "repressed_only", "biased_and_repressed", "other_protein_coding"),
  n_genes = c(
    length(setdiff(strict_biased_union, strict_repressed_union)),
    length(setdiff(strict_repressed_union, strict_biased_union)),
    length(intersect(strict_biased_union, strict_repressed_union)),
    length(setdiff(protein_coding_universe, union(strict_biased_union, strict_repressed_union)))
  ),
  stringsAsFactors = FALSE
)
strict_up_down_set_sizes_with_other <- bind_rows(
  direction_sizes %>%
    transmute(set, generation, direction, n_genes),
  data.frame(
    set = "other_protein_coding",
    generation = "other",
    direction = "other",
    n_genes = length(setdiff(protein_coding_universe, union(strict_biased_union, strict_repressed_union))),
    stringsAsFactors = FALSE
  )
)
strict_biased_with_other_sizes <- bind_rows(
  biased_set_sizes,
  data.frame(
    set = "other_protein_coding",
    n_genes = length(setdiff(protein_coding_universe, strict_biased_union)),
    stringsAsFactors = FALSE
  )
)
write_tsv(strict_biased_with_other_sizes, file.path(out_dir, "triphasy_biased_gene_set_sizes_with_other.tsv"))
write_tsv(strict_up_down_status_sizes, file.path(out_dir, "triphasy_up_down_gene_status_sizes_with_other.tsv"))
write_tsv(strict_up_down_set_sizes_with_other, file.path(out_dir, "triphasy_up_down_gene_set_sizes_with_other.tsv"))
write_tsv(gam_sensitivity, file.path(out_dir, "triphasy_gam_spermatia_sensitivity.tsv"))

write_tsv(membership_overlap, file.path(out_dir, "triphasy_overlap_gene_set_membership.tsv"))
write_tsv(
  membership_overlap %>% count(set, name = "n_genes") %>% arrange(desc(n_genes)),
  file.path(out_dir, "triphasy_overlap_gene_set_sizes.tsv")
)
write_tsv(overlap_stats, file.path(out_dir, "triphasy_overlap_pairwise_overlap_stats.tsv"))
write_tsv(focus_pairs, file.path(out_dir, "cystocarp_vs_tetrasporo_spore_overlap.tsv"))
write_tsv(enrich, file.path(out_dir, "triphasy_biased_functional_enrichment.tsv"))
write_tsv(enrich_top, file.path(out_dir, "triphasy_biased_functional_enrichment_top.tsv"))
write_tsv(down_enrich, file.path(out_dir, "triphasy_downregulated_functional_enrichment.tsv"))
write_tsv(down_enrich_top, file.path(out_dir, "triphasy_downregulated_functional_enrichment_top.tsv"))
write_tsv(down_target_off_enrich, file.path(out_dir, "triphasy_downregulated_target_off_functional_enrichment.tsv"))
write_tsv(down_target_off_enrich_top, file.path(out_dir, "triphasy_downregulated_target_off_functional_enrichment_top.tsv"))
write_tsv(rank_tab, file.path(out_dir, "triphasy_biased_phylorank_distribution.tsv"))
write_tsv(rank_tab_up_down, file.path(out_dir, "triphasy_up_down_phylorank_distribution.tsv"))

# --------- plots ----------

# P1: generation vs stage bias class counts
class_summary <- gene_bias %>% count(bias_class, name = "n_genes") %>% arrange(desc(n_genes))
p1 <- ggplot(class_summary, aes(x = reorder(bias_class, n_genes), y = n_genes)) +
  geom_col(fill = "#7A4A45", width = 0.75) +
  coord_flip() +
  labs(title = "Generation-Biased vs Stage-Biased Genes", x = NULL, y = "Genes") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "generation_stage_bias_class_barplot.pdf"), p1, width = 7.2, height = 4.2)

# P2: strict DE-biased triphasy gene sets, with all other protein-coding genes as context.
strict_plot_levels <- c("GAM_core_biased", "TETRA_biased", "CSP_biased", "other_protein_coding")
strict_plot_pal <- c(
  GAM_core_biased = "#8E2F3F",
  TETRA_biased = "#8FA33E",
  CSP_biased = "#6A3D9A",
  other_protein_coding = "#D9D2C8"
)
p2df <- strict_biased_with_other_sizes %>%
  mutate(set = factor(set, levels = strict_plot_levels))

p2 <- ggplot(p2df, aes(x = set, y = n_genes, fill = set)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  scale_fill_manual(values = strict_plot_pal) +
  labs(title = "Triphasy Group-Biased Genes in Protein-Coding Context", x = NULL, y = "Genes") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_gene_set_barplot.pdf"), p2, width = 6.2, height = 4.2)

# P3: donut view of strict DE-biased/repressed triphasy gene-set sizes, with all other protein-coding genes.
donut_plot_levels <- c(
  "GAM_core_biased", "GAM_core_downregulated",
  "CSP_biased", "CSP_downregulated",
  "TETRA_biased", "TETRA_downregulated",
  "other_protein_coding"
)
donut_plot_pal <- c(
  GAM_core_biased = "#8E2F3F",
  GAM_core_downregulated = "#C9969A",
  CSP_biased = "#6A3D9A",
  CSP_downregulated = "#B8A4D0",
  TETRA_biased = "#8FA33E",
  TETRA_downregulated = "#C9D69A",
  other_protein_coding = "#D9D2C8"
)
donut_df <- strict_up_down_set_sizes_with_other %>%
  mutate(set = factor(set, levels = donut_plot_levels)) %>%
  arrange(set) %>%
  mutate(
    fraction = n_genes / sum(n_genes),
    ymax = cumsum(fraction),
    ymin = dplyr::lag(ymax, default = 0),
    label = paste0(set, "\n", n_genes, " (", round(100 * fraction, 1), "%)")
  )

p_donut <- ggplot(donut_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.4, fill = set)) +
  geom_rect(color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  scale_fill_manual(values = donut_plot_pal) +
  labs(title = "Strict Triphasy Biased and Repressed Gene Sets in Protein-Coding Context", fill = NULL) +
  theme_void(base_size = 10) +
  theme(legend.position = "right")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_gene_set_donut.pdf"), p_donut, width = 5.5, height = 4.5)
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_gene_set_donut.pdf"), p_donut, width = 5.8, height = 4.6)

p_dir <- ggplot(direction_sizes, aes(x = generation, y = n_genes, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = n_genes),
    position = position_dodge(width = 0.72),
    vjust = -0.25,
    size = 3
  ) +
  scale_fill_manual(values = c(upregulated = "#663333", downregulated = "#6A7F8F"), name = NULL) +
  labs(
    title = "Triphasy Upregulated and Downregulated Gene Sets",
    x = NULL,
    y = "Genes"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_gene_set_barplot.pdf"), p_dir, width = 6.8, height = 4.4)

specific_fraction_values <- bind_rows(
  up_specificity_summary %>%
    filter(specificity_class == "target_detectable_only") %>%
    transmute(
      generation = recode(sub("_biased$", "", set), GAM_core = "GAM", TETRA = "TETRA", CSP = "CSP"),
      direction = "Upregulated",
      total_genes = total_genes,
      strict_genes = n_genes,
      strict_class = "target-specific",
      strict_fraction = fraction
    ),
  down_off_summary %>%
    filter(target_off_opposing_detected) %>%
    transmute(
      generation = recode(sub("_downregulated$", "", set), GAM_core = "GAM", TETRA = "TETRA", CSP = "CSP"),
      direction = "Downregulated",
      total_genes = total_genes,
      strict_genes = n_genes,
      strict_class = "target-off",
      strict_fraction = fraction
    )
) %>%
  mutate(
    other_genes = total_genes - strict_genes,
    generation = factor(generation, levels = c("GAM", "CSP", "TETRA")),
    direction = factor(direction, levels = c("Upregulated", "Downregulated"))
  ) %>%
  arrange(generation, direction)

specific_fraction_long <- bind_rows(
  specific_fraction_values %>%
    transmute(generation, direction, segment = "specific/off subset", n_genes = strict_genes, total_genes, strict_genes, strict_fraction, strict_class),
  specific_fraction_values %>%
    transmute(generation, direction, segment = "remaining genes", n_genes = other_genes, total_genes, strict_genes, strict_fraction, strict_class)
) %>%
  mutate(segment = factor(segment, levels = c("specific/off subset", "remaining genes")))

write_tsv(specific_fraction_values, file.path(out_dir, "triphasy_up_down_specific_fraction_barplot_values.tsv"))

p_specific_fraction <- ggplot(specific_fraction_long, aes(x = direction, y = n_genes, fill = segment)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(
    data = specific_fraction_values,
    aes(x = direction, y = total_genes, label = paste0(strict_genes, "/", total_genes, "\n", round(100 * strict_fraction), "%")),
    inherit.aes = FALSE,
    vjust = -0.25,
    size = 3
  ) +
  facet_wrap(~ generation, nrow = 1) +
  scale_fill_manual(values = c("specific/off subset" = "#663333", "remaining genes" = "#D9D2C8"), name = NULL) +
  labs(
    title = "Strict Subsets Within Triphasy Up/Down Gene Sets",
    subtitle = "Up: target-specific detectable only; Down: target generation off but detectable elsewhere",
    x = NULL,
    y = "Genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
ggsave(file.path(fig_dir, "triphasy_up_down_specific_fraction_square_facets.pdf"), p_specific_fraction, width = 8.4, height = 4.2)

p_specific <- ggplot(up_specificity_summary, aes(x = set, y = n_genes, fill = specificity_class)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = n_genes),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      target_detectable_only = "#663333",
      target_enriched_but_detected_elsewhere = "#B79AC8",
      opposing_detected_only = "#6A7F8F",
      not_detected_at_floor = "#D9D2C8"
    ),
    name = NULL
  ) +
  labs(
    title = "Detectability Classes Within Upregulated Triphasy Sets",
    subtitle = paste0("Target-only means target generation TPM >= ", analysis$tpm_detect_threshold, " and all opposing generations below threshold"),
    x = NULL,
    y = "Genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_target_specificity_stackedbar.pdf"), p_specific, width = 8.4, height = 4.8)

p_down_off <- down_off_summary %>%
  mutate(
    off_class = if_else(target_off_opposing_detected, "target off; detected elsewhere", "not target-off"),
    set = factor(set, levels = c("GAM_core_downregulated", "TETRA_downregulated", "CSP_downregulated"))
  ) %>%
  group_by(set, off_class) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  ggplot(aes(x = set, y = n_genes, fill = off_class)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = n_genes),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "target off; detected elsewhere" = "#663333",
      "not target-off" = "#D9D2C8"
    ),
    name = NULL
  ) +
  labs(
    title = "Target-Off Genes Within Downregulated Triphasy Sets",
    subtitle = paste0("Target-off means target generation max TPM < ", analysis$tpm_detect_threshold, " and at least one opposing generation detectable"),
    x = NULL,
    y = "Genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_downregulated_target_off_stackedbar.pdf"), p_down_off, width = 8.2, height = 4.8)

# P4: overlap heatmap for DE-biased triphasy, cystocarp-up, and spore-up sets
heat_levels <- intersect(
  c("GAM_core_biased", "TETRA_biased", "CSP_biased", "cystocarp_upregulated", "spore_upregulated"),
  names(set_map)
)
heat_label_map <- c(
  GAM_core_biased = "GAM biased",
  TETRA_biased = "Tetrasporophyte biased",
  CSP_biased = "Carposporophyte biased",
  cystocarp_upregulated = "Whole-cystocarp upregulated",
  spore_upregulated = "Spore upregulated"
)
heat_df <- overlap_stats %>%
  filter(set_a %in% heat_levels, set_b %in% heat_levels) %>%
  transmute(set_a, set_b, jaccard, fisher_or, fisher_fdr) %>%
  bind_rows(
    overlap_stats %>%
      filter(set_a %in% heat_levels, set_b %in% heat_levels) %>%
      select(set_a_orig = set_a, set_b_orig = set_b, jaccard, fisher_or, fisher_fdr) %>%
      transmute(set_a = set_b_orig, set_b = set_a_orig, jaccard, fisher_or, fisher_fdr)
  ) %>%
  bind_rows(data.frame(
    set_a = heat_levels, set_b = heat_levels, jaccard = NA_real_,
    fisher_or = NA_real_, fisher_fdr = NA_real_, stringsAsFactors = FALSE
  )) %>%
  filter(set_a %in% heat_levels, set_b %in% heat_levels) %>%
  mutate(
    set_a_label = factor(unname(heat_label_map[set_a]), levels = unname(heat_label_map[heat_levels])),
    set_b_label = factor(unname(heat_label_map[set_b]), levels = rev(unname(heat_label_map[heat_levels]))),
    log2_fisher_or = log2(fisher_or),
    log2_or_label = if_else(is.na(log2_fisher_or), "", sprintf("%.2f", log2_fisher_or)),
    fdr_sig = is.finite(fisher_fdr) & fisher_fdr <= 0.05,
    jaccard_bin = cut(
      jaccard,
      breaks = c(-Inf, 0, 0.01, 0.03, 0.07, 0.12, Inf),
      labels = c("0", "0-0.01", "0.01-0.03", "0.03-0.07", "0.07-0.12", ">0.12"),
      include.lowest = TRUE,
      right = TRUE
    ),
    jaccard_label = if_else(is.na(jaccard), "", sprintf("%.3f", jaccard))
  )
write_tsv(heat_df, file.path(out_dir, "triphasy_overlap_jaccard_heatmap_values.tsv"))

p4 <- ggplot(heat_df, aes(x = set_a_label, y = set_b_label, fill = jaccard_bin)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = jaccard_label), size = 2.6, color = "black") +
  scale_fill_manual(
    values = c(
      "0" = "#F5F5F5",
      "0-0.01" = "#E7D5BC",
      "0.01-0.03" = "#D3B089",
      "0.03-0.07" = "#BE946C",
      "0.07-0.12" = "#8F3A3A",
      ">0.12" = "#663333"
    ),
    na.value = "#D9D9D9",
    drop = FALSE
  ) +
  labs(title = "Triphasy/Spore Secondary Set Overlap (Jaccard)", x = NULL, y = NULL, fill = "Jaccard") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_overlap_jaccard_heatmap.pdf"), p4, width = 6.5, height = 5.8)

or_limit <- max(abs(heat_df$log2_fisher_or[is.finite(heat_df$log2_fisher_or)]), na.rm = TRUE)
or_limit <- min(or_limit, 8)
p4_or <- ggplot(heat_df, aes(x = set_a_label, y = set_b_label, fill = pmax(pmin(log2_fisher_or, or_limit), -or_limit))) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = log2_or_label, fontface = if_else(fdr_sig, "bold", "plain")), size = 2.6, color = "black") +
  scale_fill_gradientn(
    colors = c("#6B6B6B", "white", "#CC9966", "#663333"),
    values = scales::rescale(c(-or_limit, 0, or_limit * 0.45, or_limit)),
    limits = c(-or_limit, or_limit),
    oob = scales::squish,
    na.value = "#D9D9D9",
    name = "log2 odds\nratio"
  ) +
  labs(
    title = "Triphasy/Spore Secondary Set Overlap Enrichment",
    subtitle = "Positive = enriched overlap; negative = depleted overlap; bold labels indicate Fisher FDR <= 0.05",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_overlap_log2_oddsratio_heatmap.pdf"), p4_or, width = 6.8, height = 5.8)

# P5: functional enrichment dotplot
if (nrow(enrich_top) > 0) {
  enrich_plot <- enrich_top %>%
    filter(is.finite(fdr), fdr <= 0.05) %>%
    mutate(
      term_label = if_else(source == "GO" & !is.na(go_term_name) & nzchar(go_term_name), go_term_name, term),
      term_label = paste0(source, ": ", term_label),
      marker_group = factor(marker_group, levels = set_order)
    ) %>%
    group_by(marker_group) %>%
    slice_min(order_by = fdr, n = 6, with_ties = FALSE) %>%
    ungroup()
  if (nrow(enrich_plot) > 0) {
    enrich_plot$term_label <- factor(enrich_plot$term_label, levels = rev(unique(enrich_plot$term_label)))
    p_enrich <- ggplot(enrich_plot, aes(x = marker_group, y = term_label)) +
      geom_point(aes(size = overlap, color = -log10(pmax(fdr, .Machine$double.xmin)))) +
      scale_color_gradient(low = "#CCCCCC", high = "#663333") +
      labs(title = "Functional Enrichment of Triphasy-Biased Gene Sets", x = NULL, y = NULL, color = "-log10 FDR", size = "Overlap") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.minor = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_functional_enrichment_dotplot.pdf"), p_enrich, width = 9, height = 7)
  }
}

if (nrow(down_target_off_enrich_top) > 0) {
  target_off_enrich_plot <- down_target_off_enrich_top %>%
    filter(is.finite(fdr), fdr <= 0.05) %>%
    mutate(
      term_label = if_else(source == "GO" & !is.na(go_term_name) & nzchar(go_term_name), go_term_name, term),
      term_label = paste0(source, ": ", term_label),
      marker_group = factor(
        marker_group,
        levels = c(
          "GAM_core_target_off_downregulated",
          "TETRA_target_off_downregulated",
          "CSP_target_off_downregulated"
        )
      )
    ) %>%
    group_by(marker_group) %>%
    slice_min(order_by = fdr, n = 8, with_ties = FALSE) %>%
    ungroup()
  if (nrow(target_off_enrich_plot) > 0) {
    target_off_enrich_plot$term_label <- factor(target_off_enrich_plot$term_label, levels = rev(unique(target_off_enrich_plot$term_label)))
    p_target_off_enrich <- ggplot(target_off_enrich_plot, aes(x = marker_group, y = term_label)) +
      geom_point(aes(size = overlap, color = -log10(pmax(fdr, .Machine$double.xmin)))) +
      scale_color_gradient(low = "#CCCCCC", high = "#663333") +
      labs(
        title = "Functional Enrichment of Target-Off Downregulated Gene Sets",
        x = NULL,
        y = NULL,
        color = "-log10 FDR",
        size = "Overlap"
      ) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.minor = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_downregulated_target_off_functional_enrichment_dotplot.pdf"), p_target_off_enrich, width = 9, height = 6.5)
  }
}

# P6: phylorank composition
p_rank <- ggplot(rank_tab, aes(x = set, y = value, fill = rank_label)) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
  scale_fill_manual(values = pal_rank, breaks = rank_levels, limits = rank_levels, drop = FALSE) +
  coord_flip() +
  labs(title = "Triphasy-Biased Gene Sets: Phylorank Composition", x = NULL, y = "Relative composition", fill = "Phylorank") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_biased_phylorank_stacked_bar_relative.pdf"), p_rank, width = 9.2, height = 5.8)

p_rank_up_down <- ggplot(rank_tab_up_down, aes(x = set_label, y = value, fill = rank_label)) +
  geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
  scale_fill_manual(values = pal_rank, breaks = rank_levels, limits = rank_levels, drop = FALSE) +
  coord_flip() +
  labs(title = "Triphasy Biased and Repressed Gene Sets: Phylorank Composition", x = NULL, y = "Relative composition", fill = "Phylorank") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
# Omitted non-manuscript PDF export: ggsave(file.path(fig_dir, "triphasy_up_down_phylorank_stacked_bar_relative.pdf"), p_rank_up_down, width = 9.4, height = 6.4)

cat("Wrote:\n")
cat(" - ", file.path(out_dir, "generation_stage_bias_gene_classification.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "generation_stage_bias_class_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_core_group_de_results.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_with_spermatia_group_de_results.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_gene_set_membership.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_gene_set_sizes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_gene_set_membership.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_gene_set_sizes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_gene_set_sizes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_specific_fraction_barplot_values.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_target_specificity.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_target_specificity_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_target_off_status.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_target_off_summary.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_target_off_gene_set_membership.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_target_off_gene_set_sizes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_gene_set_sizes_with_other.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_gene_status_sizes_with_other.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_gene_set_sizes_with_other.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_gam_spermatia_sensitivity.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_overlap_gene_set_membership.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_overlap_gene_set_sizes.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_overlap_pairwise_overlap_stats.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "cystocarp_vs_tetrasporo_spore_overlap.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_functional_enrichment.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_functional_enrichment_top.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_functional_enrichment.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_functional_enrichment_top.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_target_off_functional_enrichment.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_downregulated_target_off_functional_enrichment_top.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_biased_phylorank_distribution.tsv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "triphasy_up_down_phylorank_distribution.tsv"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "generation_stage_bias_class_barplot.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_gene_set_barplot.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_gene_set_donut.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_gene_set_donut.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_gene_set_barplot.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_specific_fraction_square_facets.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_target_specificity_stackedbar.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_downregulated_target_off_stackedbar.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_overlap_jaccard_heatmap.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_overlap_log2_oddsratio_heatmap.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_downregulated_target_off_functional_enrichment_dotplot.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_biased_phylorank_stacked_bar_relative.pdf"), "\n", sep = "")
cat(" - ", file.path(fig_dir, "triphasy_up_down_phylorank_stacked_bar_relative.pdf"), "\n", sep = "")
