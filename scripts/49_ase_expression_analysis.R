############################################################
# Allele-specific expression (ASE) analysis — R side
# Complements the SNPsplit + STAR + htseq-count bash pipeline
#
# Inputs expected from the pipeline (per-sample):
#   counts/<sample>/<sample>.M.counts.txt
#   counts/<sample>/<sample>.F.counts.txt
#
# Required companion files you provide:
#   sample_metadata.tsv   (see template below)
#   gene annotation (GTF) optional but recommended
#
# Author: you + ChatGPT
# Date: 2025-09-03
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(ggridges)
  library(patchwork)
  library(rlang)
  library(GenomicRanges)
  library(rtracklayer)
})

# ---------------------------
# Helper utilities
# ---------------------------
write_tsv_smart <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (stringr::str_ends(path, "\\.gz$")) {
    con <- gzfile(path, "wb")
    on.exit(close(con), add = TRUE)
    readr::write_tsv(df, con)
  } else {
    readr::write_tsv(df, path)
  }
}

annotate_if_available <- function(df, annotations, key = "gene_id") {
  if (!is.null(annotations) && key %in% names(df)) {
    suppressMessages(dplyr::left_join(df, annotations, by = key))
  } else {
    df
  }
}

# ---------------------------
# Paths & inputs (EDIT THESE)
# ---------------------------
script_path <- normalizePath(sys.frame(1)$ofile %||% "scripts/49_ase_expression_analysis.R", mustWork = FALSE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
paths_yml <- Sys.getenv("CFG_PATHS", file.path(project_root, "config", "paths.yml"))

read_path_value <- function(key, default = NA_character_) {
  if (!file.exists(paths_yml)) return(default)
  lines <- readLines(paths_yml, warn = FALSE)
  hit <- grep(paste0("^\\s*", key, "\\s*:"), lines, value = TRUE)
  if (length(hit) == 0) return(default)
  val <- sub(paste0("^\\s*", key, "\\s*:\\s*"), "", hit[1])
  val <- gsub("^['\"]|['\"]$", "", trimws(val))
  if (!nzchar(val)) default else val
}

# Dedicated ASE output root in this project
ase_root <- Sys.getenv("ASE_ROOT", file.path(project_root, "results", "ase_mf"))
counts_dir <- file.path(ase_root, "counts")
meta_path <- Sys.getenv("ASE_META", file.path(project_root, "config", "ase_sample_metadata.tsv"))
gtf_path <- read_path_value(
  "gff_gtf",
  "data/annotation/genes.gff"
)

# Where to save outputs
plots_dir  <- file.path(ase_root, "plots")
results_dir<- file.path(ase_root, "results")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directories:\n  results: ", results_dir, "\n  plots:   ", plots_dir)

# Optional functional annotation (used to enrich outputs)
anno_path <- read_path_value(
  "annotation_tsv",
  "data/annotation/functional_annotation.tsv"
)
anno_slim <- NULL
if (!is.na(anno_path) && file.exists(anno_path)) {
  anno_base <- read.delim(file = anno_path, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
  if ("query" %in% names(anno_base)) {
    anno_base <- dplyr::mutate(anno_base, gene_id = sub("\\.t\\d+$", "", query))
  } else if ("ID" %in% names(anno_base)) {
    anno_base <- dplyr::mutate(anno_base, gene_id = sub("\\.t\\d+$", "", ID))
  } else {
    first_col <- names(anno_base)[1]
    anno_base <- dplyr::mutate(anno_base, gene_id = sub("\\.t\\d+$", "", .data[[first_col]]))
  }
  keep_cols <- c("gene_id","query","Description","Preferred_name","GOs",
                 "KEGG_ko","KEGG_Pathway","CAZy","TAP_family","Cazy_category")
  anno_slim <- dplyr::select(anno_base, dplyr::any_of(keep_cols)) %>%
    dplyr::distinct()
  message("Loaded functional annotation table: ", basename(anno_path),
          " (", dplyr::n_distinct(anno_slim$gene_id), " unique genes).")
} else {
  message("Functional annotation file not found; gene-level outputs will lack annotations.")
}

# Minimal filtering thresholds (tweak as needed)
min_total_per_sample <- 10     # min M+F counts per gene within a sample to be considered for that sample’s binomial test
min_total_per_group  <- 30     # min pooled M+F counts per gene within a group for group-level tests
min_reps_per_group   <- 2      # min replicates in a group to compute replicate summaries
wis_permutations <- suppressWarnings(as.integer(Sys.getenv("ASE_WIS_PERMUTATIONS", "1000")))
if (is.na(wis_permutations) || wis_permutations < 10) wis_permutations <- 1000

# ------------------------------------
# Sample metadata template (TSV, tabs)
# ------------------------------------
# sample    group       tissue       replicate  ploidy
# Bm3235F_vegtissue_1   veg         vegetative  1         2
# Bm3235F_vegtissue_2   veg         vegetative  2         2
# Bm3235F_cysto_1       cysto       cystocarp   1         2
# ...
meta <- read_delim(
  file = meta_path,
  delim = "\t",
  locale = locale(encoding = "UTF-8"),  # try UTF-16LE below if needed
  trim_ws = TRUE,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "\r", "") |> str_trim()))

stopifnot(all(c("sample","group") %in% names(meta)))  # must have at least these

# Handle ploidy if present
if ("ploidy" %in% names(meta)) {
  meta <- meta %>% mutate(ploidy = as.integer(ploidy))
} else {
  message("Note: 'ploidy' column not found in metadata; WIS group tests will skip ploidy-stratified permutations.")
}

# Keep all stage labels separate (11 groups) for downstream summaries/plots.
meta <- meta %>% mutate(stage11 = group)

# ---------------------------
# PCA-aligned color palettes
# ---------------------------
pca_stage_cols <- c(
  "Bm3235f_cystocarps_imm" = "#B79AC8",
  "Bm3235f_cystocarps_mat" = "#6E4B9E",
  "Bm3235c_carposporophyte" = "#8C6BB1",
  "Bm3235t_carpospores"    = "#4C2A85",
  "Bm3235f_vegtissue"      = "#8E2F3F",
  "Bm3235f_procarps"       = "#5C1F2A",
  "Bm3235f_pericarp"       = "#C96F5B",
  "Bm3235t_vegtissue"      = "#B5C46A",
  "Bm3235t_stichidia"      = "#8FA33E",
  "Bm3235t_tetraspores"    = "#6E8629",
  "Bm3235m_vegtissue"      = "#D6AF3D",
  "Bm3235m_spermatangia"   = "#B8860B",
  "Bm3235m_spermatia"      = "#8C6A00"
)

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

tissue_cols <- c(
  "cystocarps_immature" = "#B79AC8",
  "cystocarps_mature"   = "#6E4B9E",
  "carpospores"         = "#4C2A85",
  "vegtissue"           = "#8E2F3F",
  "procarps"            = "#5C1F2A",
  "pericarp"            = "#C96F5B",
  "carposporophyte"     = "#8C6BB1",
  "stichidia"           = "#8FA33E",
  "tetraspores"         = "#6E8629",
  "spermatangia"        = "#B8860B",
  "spermatia"           = "#8C6A00",
  "thallus"             = "#A67C52"
)

palette_for_values <- function(values, mapping, default_col = "#9E9E9E") {
  lv <- as.character(unique(values))
  out <- setNames(rep(default_col, length(lv)), lv)
  hit <- intersect(lv, names(mapping))
  out[hit] <- mapping[hit]
  out
}

# ---------------------------------
# 1) Read allele counts from HTSeq
# ---------------------------------
message("Scanning counts directory: ", counts_dir)
m_files <- list.files(counts_dir, pattern = "\\.M\\.counts\\.txt$", recursive = TRUE, full.names = TRUE)
f_files <- list.files(counts_dir, pattern = "\\.F\\.counts\\.txt$", recursive = TRUE, full.names = TRUE)

if (length(m_files) == 0 || length(f_files) == 0) {
  stop("No M/F counts found under: ", counts_dir, call. = FALSE)
}

read_counts <- function(path) {
  # htseq-count default format: <gene_id>\t<count>
  # skip potential special rows like "__no_feature" if present
  df <- suppressMessages(read_tsv(path, col_names = c("gene_id","count"), col_types = "ci"))
  df <- df %>% filter(!str_starts(gene_id, "__"))
  df
}

# Helper to extract sample name from path: counts/<sample>/<sample>.M.counts.txt
extract_sample <- function(p) {
  # parent dir name is the sample
  basename(dirname(p))
}

counts_long <- bind_rows(
  tibble(file = m_files, allele = "M"),
  tibble(file = f_files, allele = "F")
) %>%
  mutate(sample = map_chr(file, extract_sample)) %>%
  mutate(data = map(file, read_counts)) %>%
  select(-file) %>%
  unnest(data)

# ---------------------------------
# 2) Tidy matrices & quick QC
# ---------------------------------
# Wide matrix per sample with M/F columns
counts_wide <- counts_long %>%
  pivot_wider(names_from = allele, values_from = count, values_fill = 0) %>%
  mutate(total = M + F,
         AF_M  = if_else(total > 0, M/total, NA_real_),
         log2_aFC = if_else(total > 0, log2((M + 0.5) / (F + 0.5)), NA_real_))  # pseudocount keeps single-allele observations

# Exclude UV sex-chromosome genes from ASE analyses.
# Build a gene -> chromosome lookup from the GTF and drop genes on ChrVm/ChrUf
try({
  gtf_for_chr <- rtracklayer::import(gtf_path)
  genes_gr <- gtf_for_chr[gtf_for_chr$type == "gene"]
  gene_cols <- names(S4Vectors::mcols(genes_gr))
  gene_id_col <- intersect(c("gene_id", "ID"), gene_cols)
  if (length(gene_id_col) > 0) {
    gene_chr_tbl <- tibble(
      gene_id = as.character(S4Vectors::mcols(genes_gr)[[gene_id_col[1]]]),
      chr     = as.character(GenomicRanges::seqnames(genes_gr))
    )
  } else {
    # Fallback: derive per-gene chromosome from exons (pick predominant)
    ex_gr <- gtf_for_chr[gtf_for_chr$type == "exon"]
    if (!is.null(S4Vectors::mcols(ex_gr)$gene_id)) {
      gene_chr_tbl <- tibble(
        gene_id = S4Vectors::mcols(ex_gr)$gene_id,
        chr     = as.character(GenomicRanges::seqnames(ex_gr))
      ) %>%
      count(gene_id, chr, name = "n") %>%
      arrange(gene_id, desc(n)) %>%
      group_by(gene_id) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(gene_id, chr)
    } else {
      gene_chr_tbl <- NULL
      warning("GTF lacks 'gene_id' on genes/exons; cannot filter ChrVm genes.")
    }
  }

  if (!is.null(gene_chr_tbl)) {
    sex_chromosomes <- c("ChrVm", "ChrUf")
    excluded_genes <- gene_chr_tbl %>% filter(chr %in% sex_chromosomes) %>% pull(gene_id) %>% unique()
    n_excl <- length(excluded_genes)
    if (n_excl > 0) {
      counts_wide <- counts_wide %>% filter(!(gene_id %in% excluded_genes))
      message(sprintf("Excluded %d genes on UV sex chromosomes (ChrVm/ChrUf) from analysis.", n_excl))
    } else {
      message("No genes annotated on ChrVm/ChrUf found to exclude.")
    }
  }
}, silent = TRUE)

# Merge metadata
dat <- counts_wide %>%
  inner_join(meta, by = "sample")

# ---------------------------------
# 2c) Biotype breakdown of analyzed genes
# ---------------------------------
extract_biotype_table <- function(gr) {
  if (length(gr) == 0) return(NULL)
  gene_gr <- gr[tolower(as.character(mcols(gr)$type)) == "gene"]
  if (length(gene_gr) == 0) return(NULL)

  cn <- names(mcols(gene_gr))
  biotype_col <- intersect(c("biotype", "gene_biotype", "gene_type"), cn)
  biotype_col <- if (length(biotype_col) > 0) biotype_col[1] else NA_character_

  gene_tbl <- tibble(
    gene_id = as.character(mcols(gene_gr)$gene_id),
    biotype = if (!is.na(biotype_col)) as.character(mcols(gene_gr)[[biotype_col]]) else NA_character_
  ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    mutate(
      biotype = if_else(is.na(biotype) | biotype == "", "unknown", biotype),
      biotype = as.character(biotype)
    ) %>%
    distinct(gene_id, .keep_all = TRUE)
  gene_tbl
}

biotype_breakdown <- NULL
try({
  anno_gr <- rtracklayer::import(gtf_path)
  anno_biotype <- extract_biotype_table(anno_gr)
  if (!is.null(anno_biotype) && nrow(anno_biotype) > 0) {
    ase_genes <- dat %>% distinct(gene_id) %>% mutate(in_ase = TRUE)
    biotype_breakdown <- anno_biotype %>%
      left_join(ase_genes, by = "gene_id") %>%
      mutate(in_ase = if_else(is.na(in_ase), FALSE, in_ase)) %>%
      group_by(biotype) %>%
      summarise(
        n_genes_annotation = n(),
        n_genes_in_ase = sum(in_ase),
        frac_in_ase = if_else(n_genes_annotation > 0, n_genes_in_ase / n_genes_annotation, NA_real_),
        .groups = "drop"
      ) %>%
      arrange(desc(n_genes_in_ase), desc(n_genes_annotation), biotype)

    write_tsv_smart(biotype_breakdown, file.path(results_dir, "ase_biotype_breakdown.tsv"))
  } else {
    message("Biotype breakdown skipped: no gene-level biotype information found in annotation.")
  }
}, silent = TRUE)

# Per-sample QC
qc_sample <- dat %>%
  group_by(sample, group) %>%
  summarise(total_M = sum(M), total_F = sum(F), total = sum(total),
            mean_AF_M = mean(AF_M, na.rm = TRUE),
            genes_w_total_ge_min = sum(total >= min_total_per_sample),
            .groups = "drop") %>%
  arrange(desc(total))
write_tsv(qc_sample, file.path(results_dir, "ase_qc_per_sample.tsv"))

# Quick plots
p1 <- qc_sample %>% 
  ggplot(aes(x = reorder(sample, -total), y = total, fill = group)) +
  geom_col() +
  scale_fill_manual(values = palette_for_values(qc_sample$group, pca_stage_cols), drop = FALSE) +
  coord_flip() + labs(title = "Total allele counts per sample", x = NULL, y = "Total M+F") +
  theme_bw()

# Ridgeline distribution of log2 aFC by group (replaces histogram)
dat_ridge <- dat %>%
  dplyr::filter(total >= min_total_per_sample, is.finite(log2_aFC)) %>%
  dplyr::mutate(group = factor(group, levels = rev(c(intersect(stage_order, unique(group)), setdiff(unique(group), stage_order)))))

p2 <- ggplot(dat_ridge, aes(
        x = log2_aFC,
        y = group,
        fill = group
      )) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 1.1, rel_min_height = 0.001, color = NA) +
  scale_fill_manual(values = palette_for_values(dat_ridge$group, pca_stage_cols), drop = FALSE) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Distribution of log2 aFC by group (ridgeline)",
       x = "log2(M/F)", y = NULL)

# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "01_qc_totals.pdf"), p1, width = 8, height = 6)

## Ridgeline distribution of log2 aFC by stage (11 groups)
if ("stage11" %in% names(dat)) {
  dat_ridge_tissue <- dat %>%
    dplyr::filter(total >= min_total_per_sample, is.finite(log2_aFC)) %>%
    dplyr::mutate(stage11 = factor(stage11, levels = rev(c(intersect(stage_order, unique(stage11)), setdiff(unique(stage11), stage_order)))))

  p2_tissue <- ggplot(dat_ridge_tissue, aes(
          x = log2_aFC,
          y = stage11,
          fill = stage11
        )) +
    ggridges::geom_density_ridges(alpha = 0.7, scale = 1.1, rel_min_height = 0.001, color = NA) +
    scale_fill_manual(values = palette_for_values(dat_ridge_tissue$stage11, pca_stage_cols), drop = FALSE) +
    geom_vline(xintercept = 0, linetype = 2) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Distribution of log2 aFC by stage (ridgeline)",
         x = "log2(M/F)", y = NULL)

  ggsave(file.path(plots_dir, "02b_ridge_log2aFC_by_tissue.pdf"), p2_tissue, width = 9, height = 6)
} else {
  message("Ridge by stage skipped: 'stage11' column missing in metadata.")
}

# ---------------------------------
# 2b) Genes with informative SNPs (from pipeline SNP table)
#     Intersect fixed M/F SNPs with exons and summarise coverage
# ---------------------------------
{
  snp_table_path <- file.path(ase_root, "snps", "MF.fixed_diffs.table.tsv")
  if (!file.exists(snp_table_path)) {
    warning("Informative SNP table not found: ", snp_table_path,
            ". Skipping informative SNP coverage section.")
  } else {
    # Load informative SNPs (fixed M vs F differences)
    snps_dt <- data.table::fread(
      snp_table_path,
      col.names = c("chr","pos","ref","alt","M_allele","F_allele")
    )

    # Build GRanges for SNPs (1 bp)
    snp_gr <- GenomicRanges::GRanges(
      seqnames = snps_dt$chr,
      ranges   = IRanges::IRanges(start = snps_dt$pos, width = 1)
    )

    # Load GTF and select exons (must carry gene_id)
    gtf_gr <- rtracklayer::import(gtf_path)
    exons  <- gtf_gr[gtf_gr$type == "exon"]
    if (is.null(S4Vectors::mcols(exons)$gene_id)) {
      warning("GTF lacks 'gene_id' in exon attributes. Skipping informative SNP coverage.")
    } else {
      exons <- exons[!is.na(exons$gene_id)]
      excluded_gene_ids <- if (exists("excluded_genes", inherits = FALSE)) excluded_genes else NULL

      summarise_snp_overlap <- function(feature_gr, label, col_suffix, excluded_ids = NULL) {
        col_name <- paste0("n_informative_snps_", col_suffix)
        if (length(feature_gr) == 0) {
          message(sprintf("No %s annotations found; skipping %s overlap counts.", label))
          empty_tbl <- tibble::tibble(gene_id = character()) %>%
            dplyr::mutate(!!rlang::sym(col_name) := integer())
          return(list(
            label = label,
            col_name = col_name,
            counts = empty_tbl,
            n_snps = 0L,
            n_genes = 0L
          ))
        }
        gene_ids <- S4Vectors::mcols(feature_gr)$gene_id
        if (is.null(gene_ids)) {
          warning(sprintf("GTF lacks 'gene_id' for %s entries; skipping %s overlap.", label, label))
          return(NULL)
        }
        keep <- !is.na(gene_ids)
        if (!any(keep)) {
          warning(sprintf("All %s entries lack gene_id; skipping %s overlap.", label, label))
          return(NULL)
        }
        feature_gr <- feature_gr[keep]
        gene_ids <- gene_ids[keep]

        hits_all <- GenomicRanges::findOverlaps(snp_gr, feature_gr, ignore.strand = TRUE)
        if (!is.null(excluded_ids) && length(hits_all) > 0) {
          keep_idx <- !(gene_ids[S4Vectors::subjectHits(hits_all)] %in% excluded_ids)
          hits <- hits_all[keep_idx]
        } else {
          hits <- hits_all
        }

        snp_idx <- if (length(hits) > 0) unique(S4Vectors::queryHits(hits)) else integer()
        genes_hit <- if (length(hits) > 0) gene_ids[S4Vectors::subjectHits(hits)] else character()

        if (length(genes_hit) > 0) {
          gene_counts <- as.data.frame(table(genes_hit), stringsAsFactors = FALSE) %>%
            dplyr::rename(gene_id = genes_hit, !!rlang::sym(col_name) := Freq) %>%
            dplyr::arrange(dplyr::desc(.data[[col_name]]))
        } else {
          gene_counts <- tibble::tibble(gene_id = character()) %>%
            dplyr::mutate(!!rlang::sym(col_name) := integer())
        }

        n_genes_hit <- if (length(genes_hit) > 0) length(unique(genes_hit)) else 0L
        list(
          label = label,
          col_name = col_name,
          counts = tibble::as_tibble(gene_counts),
          n_snps = length(snp_idx),
          n_genes = n_genes_hit
        )
      }

      exon_summary <- summarise_snp_overlap(
        feature_gr = exons,
        label = "exons",
        col_suffix = "exonic",
        excluded_ids = excluded_gene_ids
      )

      type_lower <- tolower(as.character(gtf_gr$type))
      utr_gr <- gtf_gr[type_lower %in% c("five_prime_utr", "three_prime_utr")]
      utr_summary <- summarise_snp_overlap(
        feature_gr = utr_gr,
        label = "UTRs (5' + 3')",
        col_suffix = "utr",
        excluded_ids = excluded_gene_ids
      )

      feature_tables <- purrr::compact(list(
        if (!is.null(exon_summary)) exon_summary$counts else NULL,
        if (!is.null(utr_summary)) utr_summary$counts else NULL
      ))

      if (length(feature_tables) == 0) {
        warning("No informative SNP overlaps could be summarised; skipping per-gene/coverage outputs.")
      } else {
        gene_snp_counts <- purrr::reduce(feature_tables, dplyr::full_join, by = "gene_id") %>%
          dplyr::mutate(across(-gene_id, ~ tidyr::replace_na(.x, 0L)))

        feature_count_cols <- setdiff(names(gene_snp_counts), "gene_id")
        total_informative_per_gene <- if (length(feature_count_cols) > 0) {
          rowSums(gene_snp_counts[feature_count_cols], na.rm = TRUE)
        } else {
          integer(nrow(gene_snp_counts))
        }

        # High-level summaries (now that totals per gene are available)
        n_total_snps      <- nrow(snps_dt)
        n_exonic_snps     <- if (!is.null(exon_summary)) exon_summary$n_snps else 0L
        n_genes_exonic    <- if (!is.null(exon_summary)) exon_summary$n_genes else 0L
        n_utr_snps        <- if (!is.null(utr_summary)) utr_summary$n_snps else 0L
        n_genes_utr       <- if (!is.null(utr_summary)) utr_summary$n_genes else 0L
        n_genes_any       <- if (length(total_informative_per_gene) > 0) sum(total_informative_per_gene > 0) else 0L
        gtf_gene_ids      <- unique(S4Vectors::mcols(gtf_gr[gtf_gr$type == "gene"])$gene_id)
        n_total_genes_gtf <- length(gtf_gene_ids)

        summary_lines <- c(
          paste("Total informative SNPs (fixed M↔F):", n_total_snps),
          paste("Informative SNPs inside exons:",     n_exonic_snps),
          paste("Genes with ≥1 exonic informative SNP:", n_genes_exonic),
          paste("Informative SNPs inside UTRs:", n_utr_snps),
          paste("Genes with ≥1 UTR informative SNP:", n_genes_utr),
          paste("Genes with ≥1 informative SNP (any feature):", n_genes_any),
          paste("Total genes in GTF:",                n_total_genes_gtf),
          sprintf("Fraction genes w/ exonic informative SNP: %.3f",
                  ifelse(n_total_genes_gtf > 0, n_genes_exonic / n_total_genes_gtf, NA_real_)),
          sprintf("Fraction genes w/ informative SNP anywhere: %.3f",
                  ifelse(n_total_genes_gtf > 0, n_genes_any / n_total_genes_gtf, NA_real_))
        )

        readr::write_lines(summary_lines, file.path(results_dir, "informative_snp_summary.txt"))

        write_tsv_smart(
          annotate_if_available(gene_snp_counts, anno_slim),
          file.path(results_dir, "informative_snps_per_gene.tsv")
        )

        genes_with_snp <- gene_snp_counts %>%
          dplyr::mutate(total_informative = total_informative_per_gene) %>%
          dplyr::filter(total_informative > 0) %>%
          dplyr::pull(gene_id)

        # Per-sample coverage among expressed genes (same expr threshold as above)
        per_sample_coverage <- dat %>%
          dplyr::mutate(has_informative_snp = gene_id %in% genes_with_snp,
                        expressed = total >= min_total_per_sample) %>%
          dplyr::group_by(sample, group) %>%
          dplyr::summarise(
            genes_expressed        = sum(expressed, na.rm = TRUE),
            genes_expr_with_infSNP = sum(expressed & has_informative_snp, na.rm = TRUE),
            frac_expr_with_infSNP  = dplyr::if_else(genes_expressed > 0,
                                                    genes_expr_with_infSNP / genes_expressed, NA_real_),
            .groups = "drop"
          ) %>%
          dplyr::arrange(dplyr::desc(frac_expr_with_infSNP))

        readr::write_tsv(per_sample_coverage, file.path(results_dir, "per_sample_informative_snp_coverage.tsv"))

        ase_snp_summary <- tibble::tibble(
          step = factor(
            c("Fixed M/F SNPs", "SNPs in exons", "SNPs in UTRs", "Genes with exon/UTR SNP", "Genes with ASE reads", "Genes testable", "Genes in ridge plot"),
            levels = c("Fixed M/F SNPs", "SNPs in exons", "SNPs in UTRs", "Genes with exon/UTR SNP", "Genes with ASE reads", "Genes testable", "Genes in ridge plot")
          ),
          count = c(
            n_total_snps,
            n_exonic_snps,
            n_utr_snps,
            n_genes_any,
            dplyr::n_distinct(dat$gene_id[dat$total > 0]),
            dplyr::n_distinct(dat$gene_id[dat$total >= min_total_per_sample]),
            dplyr::n_distinct(dat$gene_id[dat$total >= min_total_per_sample & is.finite(dat$log2_aFC)])
          ),
          unit = c("SNPs", "SNPs", "SNPs", "Genes", "Genes", "Genes", "Genes"),
          note = c(
            "Fixed M/F informative SNPs from MF.fixed_diffs.table.tsv",
            "Informative SNPs overlapping exon features in the analysis GTF",
            "Informative SNPs overlapping 5' or 3' UTR features in the analysis GTF",
            "Genes with >=1 informative SNP in exon or UTR",
            "Genes with total M+F allele-assigned reads >0 in >=1 sample",
            paste0("Genes with total M+F allele-assigned reads >=", min_total_per_sample, " in >=1 sample"),
            paste0("Genes plotted in ridgeline after total >=", min_total_per_sample, " and finite log2(M/F)")
          )
        )
        readr::write_tsv(ase_snp_summary, file.path(results_dir, "ase_informative_snp_barplot_summary.tsv"))

        p_snp_summary <- ase_snp_summary %>%
          ggplot(aes(x = step, y = count, fill = unit)) +
          geom_col(width = 0.72, color = NA) +
          geom_text(aes(label = scales::comma(count)), hjust = -0.08, size = 3) +
          coord_flip(clip = "off") +
          scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.14))) +
          scale_fill_manual(values = c("SNPs" = "#BFA78D", "Genes" = "#663333")) +
          theme_bw(base_size = 11) +
          theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = margin(5.5, 28, 5.5, 5.5)
          ) +
          labs(
            title = "Informative SNP attrition for ASE analysis",
            x = NULL,
            y = "Count"
          )
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "03_informative_snp_attrition_barplot.pdf"), p_snp_summary, width = 7.2, height = 4.2)

        message("[Informative SNP coverage written:]\n  ",
                file.path(results_dir, "informative_snp_summary.txt"), "\n  ",
                file.path(results_dir, "informative_snps_per_gene.tsv"), "\n  ",
                file.path(results_dir, "per_sample_informative_snp_coverage.tsv"), "\n  ",
                file.path(results_dir, "ase_informative_snp_barplot_summary.tsv"), "\n  ",
                file.path(plots_dir, "03_informative_snp_attrition_barplot.pdf"))
      }
    }
  }
}

# ---------------------------------
# 3) Per-sample ASE (binomial test)
# ---------------------------------
binom_one_sample <- function(M, total, p0 = 0.5) {
  if (is.na(total) || total <= 0 || is.na(M)) return(NA_real_)
  suppressWarnings(binom.test(M, total, p = p0)$p.value)
}

per_sample_tests <- dat %>%
  filter(total >= min_total_per_sample) %>%
  mutate(p_binom = map2_dbl(M, total, binom_one_sample),
         padj    = p.adjust(p_binom, method = "BH"))

# Tissue-level imprinted gene summary
imprinted_by_tissue <- per_sample_tests %>%
  filter(
    total >= min_total_per_sample,
    !is.na(padj) & padj < 0.05,
    abs(log2_aFC) >= 1,
    AF_M >= 0.9 | AF_M <= 0.1
  ) %>%
  group_by(group, gene_id) %>%
  summarise(
    n_samples = n(),
    n_reps = n_distinct(sample),
    mean_AF_M = mean(AF_M, na.rm = TRUE),
    direction = if_else(mean_AF_M > 0.5, "M_dominant", "F_dominant"),
    all_reps_biased = n_samples == n_reps,
    .groups = "drop"
  ) %>%
  filter(all_reps_biased) %>%
  arrange(group, direction, desc(abs(mean_AF_M - 0.5)))

write_tsv_smart(
  annotate_if_available(imprinted_by_tissue, anno_slim),
  file.path(results_dir, "imprinted_gene_candidates_by_tissue.tsv")
)

# ---------------------------------
# 4) Replicate-level summaries within groups
#    (depth-weighted AF, pooled-binomial, quasi-binomial)
# ---------------------------------

# 4a) Depth-weighted AF per (group, gene)
group_pooled <- dat %>%
  group_by(group, gene_id) %>%
  summarise(M = sum(M), F = sum(F), total = sum(total), n_rep = n_distinct(sample), .groups = "drop") %>%
  mutate(AF_M = if_else(total > 0, M/total, NA_real_),
         log2_aFC = if_else(total > 0, log2((M + 0.5) / (F + 0.5)), NA_real_))

# 4b) Binomial test on pooled counts (simple but can be overconfident)
group_pooled <- group_pooled %>%
  mutate(p_binom_pool = if_else(total >= min_total_per_group,
                                map2_dbl(M, total, binom_one_sample),
                                NA_real_))
group_pooled <- group_pooled %>%
  group_by(group) %>%
  mutate(padj_pool = p.adjust(p_binom_pool, method = "BH")) %>%
  ungroup()

# 4c) Quasi-binomial model to absorb overdispersion across replicates
#     For each (group, gene), fit: cbind(M, F) ~ 1 with family = quasibinomial(link="logit")
#     Wald test: intercept vs logit(0.5) = 0  (so test coef ~= 0)
#
#     NB: This requires >= 2 replicates per group-gene with total > 0.
quasibinom_test <- function(df_one_gene_group) {
  df <- df_one_gene_group %>% filter(total > 0)
  if (nrow(df) < min_reps_per_group) return(tibble(est_logit = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_))
  # model on proportions with weights = total
  fit <- suppressWarnings(glm(cbind(M, F) ~ 1, family = quasibinomial(link = "logit"), data = df))
  co <- coef(summary(fit))
  # intercept estimate on logit scale; test vs 0
  est <- unname(co[1, "Estimate"])
  se  <- unname(co[1, "Std. Error"])
  z   <- ifelse(is.finite(se) && se > 0, est / se, NA_real_)
  p   <- ifelse(is.finite(z), 2 * pnorm(-abs(z)), NA_real_)
  tibble(est_logit = est, se = se, z = z, p = p)
}

group_quasi <- dat %>%
  group_by(group, gene_id) %>%
  group_modify(~ quasibinom_test(.x)) %>%
  ungroup() %>%
  mutate(padj_quasi = p.adjust(p, method = "BH"))

# 4d) Merge pooled & quasi results
group_results <- group_pooled %>%
  left_join(group_quasi, by = c("group","gene_id")) %>%
  arrange(group, padj_quasi, padj_pool)

# Pre-compute WIS per sample (reused by group/tissue summaries)
compute_wis_tbl <- function(dat, meta, min_total) {
  dat %>%
    dplyr::filter(total >= min_total) %>%
    dplyr::group_by(sample, group) %>%
    dplyr::summarise(
      WIS   = sum(abs((M/total) - 0.5) * total) / sum(total),
      depth = sum(total),
      genes = n(),
      .groups = "drop"
    ) %>%
    dplyr::left_join(meta %>% dplyr::select(any_of(c("sample","ploidy"))), by = "sample")
}

wis_tbl <- compute_wis_tbl(dat, meta, min_total_per_sample)
write_tsv_smart(wis_tbl, file.path(results_dir, "wis_per_sample.tsv"))

# ---------------------------------
# 4f) Stage-level ASE burden: WIS (diploids), with permutation test
# ---------------------------------
{
  if (!"stage11" %in% names(meta)) {
    message("Stage-level WIS skipped: 'stage11' column missing in metadata.")
  } else {
    # Add stage to the existing per-sample WIS table
    wis_tbl_t <- wis_tbl %>% left_join(meta %>% select(sample, stage11, ploidy), by = "sample")
    readr::write_tsv(wis_tbl_t, file.path(results_dir, "wis_per_sample_with_tissue.tsv"))

    # Focus on diploids for inference
    dip_t <- if ("ploidy" %in% names(wis_tbl_t)) dplyr::filter(wis_tbl_t, ploidy == 2) else wis_tbl_t

    # Plot WIS by stage (diploids)
    if (nrow(dip_t) > 0) {
      p_wis_t <- dip_t %>%
        ggplot(aes(x = stage11, y = WIS, fill = stage11)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
        scale_fill_manual(values = palette_for_values(dip_t$stage11, pca_stage_cols), drop = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Weighted Imbalance Score by stage (diploids)", y = "WIS (0–0.5)", x = NULL) +
        coord_cartesian(ylim = c(0, 0.5))
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "06_wis_by_tissue_diploid.pdf"), p_wis_t, width = 9, height = 5)
    }

    # Permutation test across stages: is a stage's mean WIS unusually high?
    if (nrow(dip_t) >= 4 && dplyr::n_distinct(dip_t$stage11) >= 2) {
      set.seed(1)
      B <- wis_permutations

      obs_t <- dip_t %>%
        group_by(stage11) %>%
        summarise(mean_WIS = mean(WIS),
                  n = n(),
                  mean_WIS_weighted = sum(WIS * depth) / sum(depth),
                  .groups = "drop")

      W <- dip_t$WIS
      TISS <- dip_t$stage11

      perm_res_t <- vector("list", B)
      for (b in seq_len(B)) {
        Tp <- sample(TISS, replace = FALSE)
        pm <- tapply(W, Tp, mean)
        perm_res_t[[b]] <- tibble(stage11 = names(pm), mean_WIS = as.numeric(pm))
      }
      perm_long_t <- dplyr::bind_rows(perm_res_t, .id = "perm")

      pvals_t <- perm_long_t %>%
        group_by(stage11) %>%
        summarise(
          p_emp = {
            m_obs <- obs_t$mean_WIS[match(dplyr::first(stage11), obs_t$stage11)]
            if (is.na(m_obs)) NA_real_ else (sum(mean_WIS >= m_obs) + 1) / (n() + 1)
          },
          .groups = "drop"
        )

      wis_tissue_summary <- obs_t %>% left_join(pvals_t, by = "stage11") %>%
        mutate(padj = p.adjust(p_emp, method = "BH")) %>%
        arrange(padj, desc(mean_WIS))

      readr::write_tsv(wis_tissue_summary, file.path(results_dir, "wis_tissue_summary.tsv"))
    } else {
      message("Stage-level WIS permutation test skipped: need >= 4 diploid samples and >= 2 stages.")
    }
  }
}

# 4e) Group-level ASE burden: Weighted Imbalance Score (WIS)
#     WIS = weighted mean |AF_M - 0.5| per sample; then compare groups
{
  # Plot WIS for diploid samples by group
  if ("ploidy" %in% names(wis_tbl)) {
    wis_dip <- wis_tbl %>% filter(ploidy == 2)
  } else {
    wis_dip <- wis_tbl
  }
  
  if (nrow(wis_dip) > 0) {
    p_wis <- wis_dip %>%
      ggplot(aes(x = group, y = WIS, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      scale_fill_manual(values = palette_for_values(wis_dip$group, pca_stage_cols), drop = FALSE) +
      theme_bw() +
      labs(title = "Weighted Imbalance Score (diploids)", y = "WIS (0–0.5)", x = NULL) +
      coord_cartesian(ylim = c(0, 0.5))
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "05_wis_by_group_diploid.pdf"), p_wis, width = 8, height = 5)
  }
  
  # Permutation test: is a group's mean WIS higher than expected by chance?
  # Stratify by ploidy (diploid only); shuffle group labels among diploids
  if ("ploidy" %in% names(wis_tbl)) {
    dip <- wis_tbl %>% filter(ploidy == 2)
  } else {
    dip <- wis_tbl
  }
  
  if (nrow(dip) >= 4 && dplyr::n_distinct(dip$group) >= 2) {
    set.seed(1)
    B <- wis_permutations
    
    obs <- dip %>%
      group_by(group) %>%
      summarise(mean_WIS = mean(WIS),
                n = n(),
                mean_WIS_weighted = sum(WIS * depth) / sum(depth),
                .groups = "drop")
    
    # Prepare vectors for fast permutation
    W <- dip$WIS
    G <- dip$group
    
    # Permute group labels and compute mean WIS per permuted group
    perm_res <- vector("list", B)
    for (b in seq_len(B)) {
      Gp <- sample(G, replace = FALSE)
      pm <- tapply(W, Gp, mean)
      perm_res[[b]] <- tibble(group = names(pm), mean_WIS = as.numeric(pm))
    }
    perm_long <- dplyr::bind_rows(perm_res, .id = "perm")
    
    # Empirical one-sided p-value: P(mean_WIS_perm >= mean_WIS_obs)
    pvals <- perm_long %>%
      group_by(group) %>%
      summarise(
        p_emp = {
          m_obs <- obs$mean_WIS[match(dplyr::first(group), obs$group)]
          if (is.na(m_obs)) NA_real_ else (sum(mean_WIS >= m_obs) + 1) / (n() + 1)
        },
        .groups = "drop"
      )
    
    wis_group_summary <- obs %>% left_join(pvals, by = "group") %>%
      mutate(padj = p.adjust(p_emp, method = "BH")) %>%
      arrange(padj, desc(mean_WIS))
    
    readr::write_tsv(wis_group_summary, file.path(results_dir, "wis_group_summary.tsv"))
  } else {
    message("WIS permutation test skipped: need >= 4 diploid samples and >= 2 groups.")
  }
}

# ---------------------------------
# 4g) ASE-positive genes per sample/stage (diploids): counts, proportions, tests, plots
# ---------------------------------
{
  alpha_ase <- 0.05

  ase_sample <- per_sample_tests %>%
    filter(total >= min_total_per_sample) %>%
    mutate(hit = !is.na(padj) & padj < alpha_ase) %>%
    group_by(sample, group) %>%
    summarise(
      genes_tested = n(),
      ase_hits = sum(hit, na.rm = TRUE),
      prop = if_else(genes_tested > 0, ase_hits / genes_tested, NA_real_),
      .groups = "drop"
    ) %>%
    left_join(meta %>% select(any_of(c("sample","stage11","ploidy"))), by = "sample")

  readr::write_tsv(ase_sample, file.path(results_dir, "ase_hits_per_sample.tsv"))

  # Diploid-only
  ase_dip <- if ("ploidy" %in% names(ase_sample)) dplyr::filter(ase_sample, ploidy == 2) else ase_sample
  ase_dip <- ase_dip %>% filter(!is.na(prop))

  # Plots: stage and group (diploids)
  if (nrow(ase_dip) > 0 && "stage11" %in% names(ase_dip)) {
    p_prop_tissue <- ase_dip %>%
      ggplot(aes(x = stage11, y = prop, fill = stage11)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      scale_fill_manual(values = palette_for_values(ase_dip$stage11, pca_stage_cols), drop = FALSE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Fraction of ASE-positive genes by stage (diploids)",
           y = "Proportion ASE-positive (FDR < 0.05)", x = NULL) +
      coord_cartesian(ylim = c(0, 1))
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "07_ase_prop_by_tissue_diploid.pdf"), p_prop_tissue, width = 9, height = 5)
  }

  if (nrow(ase_dip) > 0) {
    p_prop_group <- ase_dip %>%
      ggplot(aes(x = group, y = prop, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      scale_fill_manual(values = palette_for_values(ase_dip$group, pca_stage_cols), drop = FALSE) +
      theme_bw() +
      labs(title = "Fraction of ASE-positive genes by group (diploids)",
           y = "Proportion ASE-positive (FDR < 0.05)", x = NULL) +
      coord_cartesian(ylim = c(0, 1))
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "08_ase_prop_by_group_diploid.pdf"), p_prop_group, width = 8, height = 5)
  }

  # Permutation test across stages on mean proportion (diploids)
  if (nrow(ase_dip) >= 4 && "stage11" %in% names(ase_dip) && dplyr::n_distinct(ase_dip$stage11) >= 2) {
    set.seed(1)
    B <- wis_permutations

    obs_pt <- ase_dip %>%
      group_by(stage11) %>%
      summarise(mean_prop = mean(prop), n = n(),
                mean_prop_weighted = sum(prop * genes_tested) / sum(genes_tested),
                .groups = "drop")

    P <- ase_dip$prop
    TISS <- ase_dip$stage11

    perm_res <- vector("list", B)
    for (b in seq_len(B)) {
      Tp <- sample(TISS, replace = FALSE)
      pm <- tapply(P, Tp, mean)
      perm_res[[b]] <- tibble(stage11 = names(pm), mean_prop = as.numeric(pm))
    }
    perm_long <- dplyr::bind_rows(perm_res, .id = "perm")

    pvals <- perm_long %>%
      group_by(stage11) %>%
      summarise(
        p_emp = {
          m_obs <- obs_pt$mean_prop[match(dplyr::first(stage11), obs_pt$stage11)]
          if (is.na(m_obs)) NA_real_ else (sum(mean_prop >= m_obs) + 1) / (n() + 1)
        },
        .groups = "drop"
      )

    ase_prop_tissue_summary <- obs_pt %>% left_join(pvals, by = "stage11") %>%
      mutate(padj = p.adjust(p_emp, method = "BH")) %>%
      arrange(padj, desc(mean_prop))

    readr::write_tsv(ase_prop_tissue_summary, file.path(results_dir, "ase_prop_tissue_summary.tsv"))
  } else {
    message("ASE proportion permutation test (stage) skipped: need >= 4 diploid samples and >= 2 stages.")
  }
}

# Quick volcano for a group (edit group name below)
one_group <- unique(meta$group)[1]
v1 <- group_results %>%
  filter(group == one_group, total >= min_total_per_group) %>%
  mutate(signif = padj_quasi < 0.05) %>%
  ggplot(aes(x = log2_aFC, y = -log10(padj_quasi), colour = signif)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_colour_manual(values = c(`TRUE` = "#8E2F3F", `FALSE` = "#BDBDBD")) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() + labs(title = paste0("ASE (quasi-binomial) — group: ", one_group),
                    x = "log2 aFC (M/F)", y = "-log10(FDR)")
# Omitted non-manuscript PDF export: 03_volcano_quasi_<group>.pdf

# ---------------------------------
# 5) Differential ASE between groups (optional)
#    Test whether AF differs between (exactly) two groups.
#    Model: cbind(M, F) ~ group  (quasi-binomial), Wald test for group effect.
# ---------------------------------
groups <- unique(meta$group)
if (length(groups) == 2) {
  gA <- groups[1]; gB <- groups[2]
  dat2 <- dat %>% filter(group %in% c(gA, gB))
  diffase <- dat2 %>%
    group_by(gene_id) %>%
    group_modify(function(df, key) {
      df <- df %>% filter(total > 0)
      # Need at least one replicate per group
      if (length(unique(df$group)) < 2 ||
          min(table(df$group)) < 1) {
        return(tibble(beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_))
      }
      fit <- suppressWarnings(glm(cbind(M, F) ~ group, family = quasibinomial(link = "logit"), data = df))
      co <- coef(summary(fit))
      # group effect is the coefficient for the second level
      # relevel so gA is reference
      df$group <- relevel(factor(df$group), ref = gA)
      fit <- suppressWarnings(glm(cbind(M, F) ~ group, family = quasibinomial(link = "logit"), data = df))
      co <- coef(summary(fit))
      if (!("group" %in% rownames(co))) {
        # glm may code as groupgB; get the second row programmatically
        idx <- setdiff(rownames(co), "(Intercept)")[1]
      } else {
        idx <- "group"
      }
      beta <- unname(co[idx, "Estimate"])
      se   <- unname(co[idx, "Std. Error"])
      z    <- ifelse(is.finite(se) && se > 0, beta / se, NA_real_)
      p    <- ifelse(is.finite(z), 2 * pnorm(-abs(z)), NA_real_)
      tibble(beta = beta, se = se, z = z, p = p)
    }) %>%
    ungroup() %>%
    mutate(padj = p.adjust(p, method = "BH"),
           groupA = gA, groupB = gB)
  
  pv <- diffase %>%
    ggplot(aes(x = beta, y = -log10(padj))) +
    geom_point(alpha = 0.7, size = 1.2, color = "#5C1F2A") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme_bw() + labs(title = paste0("Differential ASE (quasi-binomial): ", gA, " vs ", gB),
                      x = "Group effect (logit AF_M)", y = "-log10(FDR)")
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "04_diffASE_quasi.pdf"), pv, width = 7, height = 5)
}

# ---------------------------------
# 5b) Extended ASE integrations
#     - Cystocarp partitioning (maternal-only vs diploid-balanced signal)
#     - Marker/hub gene ASE tests
#     - Parent-specific TAI
#     - Overlay parental bias onto WGCNA modules
#     - Immature vs mature cystocarp parent-contribution comparison
# ---------------------------------
{
  message("Running extended ASE integrations...")

  module_membership_path <- file.path(project_root, "results", "modules", "module_membership_labeled_annotated.tsv")
  marker_table_path <- file.path(project_root, "results", "markers", "marker_table.tsv")

  read_tsv_maybe <- function(path) {
    if (!file.exists(path)) return(NULL)
    suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  }

  module_tbl <- read_tsv_maybe(module_membership_path)
  marker_tbl <- read_tsv_maybe(marker_table_path)

  module_gene_tbl <- NULL
  if (!is.null(module_tbl) && all(c("gene_id", "module_label") %in% names(module_tbl))) {
    module_gene_tbl <- module_tbl %>%
      mutate(gene_id = as.character(gene_id),
             module_label = as.character(module_label))
  }

  # ---------------------------
  # 5b-1) Cystocarp read partition
  # ---------------------------
  cyst_groups <- c("Bm3235f_cystocarps_imm", "Bm3235f_cystocarps_mat")
  cyst_dat <- dat %>%
    filter(group %in% cyst_groups, total > 0) %>%
    mutate(
      maternal_only_reads = pmax(M - F, 0),
      diploid_balanced_reads = 2 * pmin(M, F),
      F_excess_reads = pmax(F - M, 0),
      maternal_only_frac = maternal_only_reads / total,
      diploid_balanced_frac = diploid_balanced_reads / total
    )

  if (nrow(cyst_dat) > 0) {
    cyst_gene_partition <- cyst_dat %>%
      group_by(sample, group, gene_id) %>%
      summarise(
        total_reads = sum(total),
        maternal_only_reads = sum(maternal_only_reads),
        diploid_balanced_reads = sum(diploid_balanced_reads),
        F_excess_reads = sum(F_excess_reads),
        maternal_only_frac = if_else(total_reads > 0, maternal_only_reads / total_reads, NA_real_),
        diploid_balanced_frac = if_else(total_reads > 0, diploid_balanced_reads / total_reads, NA_real_),
        .groups = "drop"
      )

    cyst_sample_partition <- cyst_dat %>%
      group_by(sample, group) %>%
      summarise(
        total_reads = sum(total),
        maternal_only_reads = sum(maternal_only_reads),
        diploid_balanced_reads = sum(diploid_balanced_reads),
        F_excess_reads = sum(F_excess_reads),
        maternal_only_frac = if_else(total_reads > 0, maternal_only_reads / total_reads, NA_real_),
        diploid_balanced_frac = if_else(total_reads > 0, diploid_balanced_reads / total_reads, NA_real_),
        .groups = "drop"
      )

    cyst_group_partition <- cyst_sample_partition %>%
      group_by(group) %>%
      summarise(
        n_samples = n(),
        mean_maternal_only_frac = mean(maternal_only_frac, na.rm = TRUE),
        mean_diploid_balanced_frac = mean(diploid_balanced_frac, na.rm = TRUE),
        sd_maternal_only_frac = sd(maternal_only_frac, na.rm = TRUE),
        sd_diploid_balanced_frac = sd(diploid_balanced_frac, na.rm = TRUE),
        .groups = "drop"
      )

    # Immature vs mature test on sample-level maternal-only fraction
    cyst_imm <- cyst_sample_partition %>% filter(group == "Bm3235f_cystocarps_imm") %>% pull(maternal_only_frac)
    cyst_mat <- cyst_sample_partition %>% filter(group == "Bm3235f_cystocarps_mat") %>% pull(maternal_only_frac)
    cyst_wilcox <- tibble(
      metric = "maternal_only_frac",
      group_a = "Bm3235f_cystocarps_imm",
      group_b = "Bm3235f_cystocarps_mat",
      n_a = length(cyst_imm),
      n_b = length(cyst_mat),
      median_a = ifelse(length(cyst_imm) > 0, median(cyst_imm, na.rm = TRUE), NA_real_),
      median_b = ifelse(length(cyst_mat) > 0, median(cyst_mat, na.rm = TRUE), NA_real_),
      p_value = ifelse(length(cyst_imm) >= 2 && length(cyst_mat) >= 2,
                       suppressWarnings(wilcox.test(cyst_imm, cyst_mat, exact = FALSE)$p.value),
                       NA_real_)
    )

    write_tsv_smart(cyst_gene_partition,  file.path(results_dir, "cystocarp_partition_gene_level.tsv.gz"))
    write_tsv_smart(cyst_sample_partition, file.path(results_dir, "cystocarp_partition_sample_level.tsv"))
    write_tsv_smart(cyst_group_partition,  file.path(results_dir, "cystocarp_partition_group_summary.tsv"))
    write_tsv_smart(cyst_wilcox,           file.path(results_dir, "cystocarp_partition_imm_vs_mat_stats.tsv"))

    p_cyst_stack <- cyst_sample_partition %>%
      select(sample, group, maternal_only_reads, diploid_balanced_reads, F_excess_reads) %>%
      pivot_longer(cols = c(maternal_only_reads, diploid_balanced_reads, F_excess_reads),
                   names_to = "component", values_to = "reads") %>%
      mutate(component = factor(component,
                                levels = c("maternal_only_reads", "diploid_balanced_reads", "F_excess_reads"),
                                labels = c("Maternal-only", "Diploid-balanced", "F-excess"))) %>%
      ggplot(aes(x = sample, y = reads, fill = component)) +
      geom_col(position = "fill") +
      facet_wrap(~group, scales = "free_x", nrow = 1) +
      scale_fill_manual(values = c("Maternal-only" = "#8E2F3F",
                                   "Diploid-balanced" = "#B79AC8",
                                   "F-excess" = "#C9C9C9")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = "Cystocarp read partition by sample",
           x = NULL, y = "Fraction of reads", fill = NULL)
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "09_cystocarp_partition_stacked_fraction.pdf"), p_cyst_stack, width = 10, height = 4.5)

    p_cyst_box <- cyst_sample_partition %>%
      ggplot(aes(x = group, y = maternal_only_frac, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.12, alpha = 0.8, size = 2) +
      scale_fill_manual(values = palette_for_values(cyst_sample_partition$group, pca_stage_cols), drop = FALSE) +
      theme_bw() +
      labs(title = "Maternal-only signal in cystocarp samples",
           x = NULL, y = "Maternal-only fraction of reads")
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "10_cystocarp_maternal_fraction_imm_vs_mat.pdf"), p_cyst_box, width = 6.5, height = 4.5)
  } else {
    message("Cystocarp partition skipped: no cystocarp rows found in ASE data.")
  }

  # ---------------------------
  # 5b-2) ASE tests for marker genes and module hub genes
  # ---------------------------
  ase_positive <- group_results %>%
    mutate(ase_pos = !is.na(padj_quasi) & padj_quasi < 0.05 & !is.na(log2_aFC) & abs(log2_aFC) >= 1)

  run_set_enrichment <- function(set_tbl, set_name) {
    if (is.null(set_tbl) || nrow(set_tbl) == 0) return(NULL)
    out <- ase_positive %>%
      mutate(in_set = gene_id %in% set_tbl$gene_id) %>%
      group_by(group) %>%
      group_modify(~{
        d <- .x
        a <- sum(d$ase_pos & d$in_set, na.rm = TRUE)     # set + ASE+
        b <- sum((!d$ase_pos) & d$in_set, na.rm = TRUE)  # set + ASE-
        c <- sum(d$ase_pos & (!d$in_set), na.rm = TRUE)  # non-set + ASE+
        dd <- sum((!d$ase_pos) & (!d$in_set), na.rm = TRUE)
        mat <- matrix(c(a, b, c, dd), nrow = 2, byrow = TRUE)
        p <- suppressWarnings(fisher.test(mat, alternative = "greater"))
        tibble(
          gene_set = set_name,
          n_set = a + b,
          n_background = c + dd,
          n_set_ase_pos = a,
          frac_set_ase_pos = if_else((a + b) > 0, a / (a + b), NA_real_),
          odds_ratio = unname(p$estimate),
          p_value = p$p.value
        )
      }) %>%
      ungroup() %>%
      mutate(padj = p.adjust(p_value, method = "BH"))
    out
  }

  marker_set <- NULL
  if (!is.null(marker_tbl) && "gene_id" %in% names(marker_tbl)) {
    marker_set <- marker_tbl %>%
      mutate(gene_id = as.character(gene_id),
             marker_flag = tolower(as.character(marker_flag)) %in% c("true","t","1","yes","y")) %>%
      filter(marker_flag) %>%
      distinct(gene_id)
  }

  hub_set <- NULL
  if (!is.null(module_tbl) && all(c("gene_id", "hub_flag") %in% names(module_tbl))) {
    hub_set <- module_tbl %>%
      mutate(gene_id = as.character(gene_id),
             hub_flag = tolower(as.character(hub_flag)) %in% c("true","t","1","yes","y")) %>%
      filter(hub_flag) %>%
      distinct(gene_id)
  }

  marker_enrich <- run_set_enrichment(marker_set, "marker_genes")
  hub_enrich <- run_set_enrichment(hub_set, "module_hub_genes")
  marker_hub_enrich <- bind_rows(marker_enrich, hub_enrich)
  if (nrow(marker_hub_enrich) > 0) {
    write_tsv_smart(marker_hub_enrich, file.path(results_dir, "ase_marker_hub_enrichment_by_group.tsv"))
  }

  marker_group_test <- NULL
  if (!is.null(marker_set) && nrow(marker_set) > 0) {
    marker_group_test <- ase_positive %>%
      filter(gene_id %in% marker_set$gene_id) %>%
      select(group, gene_id, M, F, total, AF_M, log2_aFC, p_binom_pool, padj_pool, p, padj_quasi) %>%
      arrange(group, padj_quasi, padj_pool)
    write_tsv_smart(marker_group_test, file.path(results_dir, "ase_marker_genes_group_tests.tsv.gz"))
  }

  hub_group_test <- NULL
  if (!is.null(hub_set) && nrow(hub_set) > 0) {
    hub_group_test <- ase_positive %>%
      filter(gene_id %in% hub_set$gene_id) %>%
      left_join(module_tbl %>% select(gene_id, module_label, module_rank, hub_flag) %>% distinct(), by = "gene_id") %>%
      select(group, module_label, module_rank, gene_id, M, F, total, AF_M, log2_aFC, p_binom_pool, padj_pool, p, padj_quasi) %>%
      arrange(group, module_rank, padj_quasi, padj_pool)
    write_tsv_smart(hub_group_test, file.path(results_dir, "ase_module_hub_genes_group_tests.tsv.gz"))
  }

  # ---------------------------
  # 5b-3) Parent-specific TAI from ASE counts + per-gene phylorank
  # ---------------------------
  phyl_tbl <- NULL
  gene_age_path <- read_path_value(
    "gene_age_tsv",
    "data/annotation/gene_ages.tsv"
  )
  if (!is.na(gene_age_path) && file.exists(gene_age_path)) {
    ga <- suppressMessages(readr::read_tsv(
      gene_age_path,
      col_names = FALSE,
      show_col_types = FALSE,
      progress = FALSE
    ))
    # Expected schema (no header): transcript_id, phylostrata_name, rank, support
    if (ncol(ga) >= 3) {
      phyl_tbl <- ga %>%
        transmute(
          gene_id = sub("\\.t\\d+$", "", as.character(ga[[1]])),
          phylorank = suppressWarnings(as.numeric(ga[[3]]))
        ) %>%
        filter(!is.na(gene_id), gene_id != "", !is.na(phylorank)) %>%
        distinct(gene_id, .keep_all = TRUE)
      write_tsv_smart(
        phyl_tbl %>% count(phylorank, name = "n_genes") %>% arrange(phylorank),
        file.path(results_dir, "ase_parental_tai_phylorank_source_summary.tsv")
      )
    }
  }

  if (!is.null(phyl_tbl) && nrow(phyl_tbl) > 0) {
    parent_tai_by_sample <- dat %>%
      select(sample, group, gene_id, M, F, total) %>%
      left_join(phyl_tbl, by = "gene_id") %>%
      filter(!is.na(phylorank)) %>%
      group_by(sample, group) %>%
      summarise(
        n_genes = n_distinct(gene_id),
        M_total = sum(M, na.rm = TRUE),
        F_total = sum(F, na.rm = TRUE),
        TAI_M = if_else(M_total > 0, sum(phylorank * M, na.rm = TRUE) / M_total, NA_real_),
        TAI_F = if_else(F_total > 0, sum(phylorank * F, na.rm = TRUE) / F_total, NA_real_),
        TAI_delta_M_minus_F = TAI_M - TAI_F,
        .groups = "drop"
      )

    parent_tai_by_group <- parent_tai_by_sample %>%
      group_by(group) %>%
      summarise(
        n_samples = n(),
        mean_TAI_M = mean(TAI_M, na.rm = TRUE),
        mean_TAI_F = mean(TAI_F, na.rm = TRUE),
        mean_delta_M_minus_F = mean(TAI_delta_M_minus_F, na.rm = TRUE),
        .groups = "drop"
      )

    write_tsv_smart(parent_tai_by_sample, file.path(results_dir, "ase_parental_tai_by_sample.tsv"))
    write_tsv_smart(parent_tai_by_group,  file.path(results_dir, "ase_parental_tai_by_group.tsv"))

    p_tai_parent <- parent_tai_by_sample %>%
      pivot_longer(cols = c(TAI_M, TAI_F), names_to = "parent", values_to = "TAI") %>%
      mutate(parent = recode(parent, TAI_M = "Maternal (M)", TAI_F = "Paternal/Female (F)")) %>%
      ggplot(aes(x = group, y = TAI, colour = parent, group = interaction(sample, parent))) +
      geom_point(size = 2.2, alpha = 0.9, position = position_dodge(width = 0.3)) +
      geom_line(alpha = 0.35, position = position_dodge(width = 0.3)) +
      scale_color_manual(values = c("Maternal (M)" = "#8E2F3F", "Paternal/Female (F)" = "#4C2A85")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Parent-specific TAI by sample and group", x = NULL, y = "TAI (ASE-weighted)")
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "11_parental_tai_by_group.pdf"), p_tai_parent, width = 11, height = 4.8)
  } else {
    message("Parent-specific TAI skipped: no per-gene phylorank column ('Number') found.")
  }

  # ---------------------------
  # 5b-4) Overlay parental bias onto WGCNA modules
  # ---------------------------
  if (!is.null(module_gene_tbl)) {
    module_bias <- ase_positive %>%
      left_join(module_gene_tbl %>% select(gene_id, module_label, module_rank, hub_flag) %>% distinct(),
                by = "gene_id") %>%
      filter(!is.na(module_label)) %>%
      mutate(ase_pos = !is.na(padj_quasi) & padj_quasi < 0.05 & !is.na(log2_aFC) & abs(log2_aFC) >= 1) %>%
      group_by(module_label, module_rank, group) %>%
      summarise(
        n_genes = n(),
        mean_AF_M = mean(AF_M, na.rm = TRUE),
        mean_log2_aFC = mean(log2_aFC, na.rm = TRUE),
        mean_abs_log2_aFC = mean(abs(log2_aFC), na.rm = TRUE),
        frac_ase_pos = mean(ase_pos, na.rm = TRUE),
        .groups = "drop"
      )

    write_tsv_smart(module_bias, file.path(results_dir, "ase_module_parental_bias_overlay.tsv.gz"))

    p_module_bias <- module_bias %>%
      mutate(
        module_label = factor(module_label, levels = unique(module_label[order(module_rank, module_label)])),
        group = factor(group, levels = unique(meta$group))
      ) %>%
      ggplot(aes(x = group, y = module_label, fill = mean_log2_aFC)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_gradient2(low = "#4C2A85", mid = "#F7F7F7", high = "#8E2F3F", midpoint = 0) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Parental bias overlay across WGCNA modules",
           x = NULL, y = "Module", fill = "Mean log2 aFC\n(M/F)")
# Omitted non-manuscript PDF export: ggsave(file.path(plots_dir, "12_module_parental_bias_heatmap.pdf"), p_module_bias, width = 11, height = 7.5)
  } else {
    message("Module parental-bias overlay skipped: module membership table not available.")
  }
}

# ---------------------------------
# 6) Optional: join functional annotations
# ---------------------------------
if (!is.null(anno_slim)) {
  matched <- length(intersect(unique(per_sample_tests$gene_id), unique(anno_slim$gene_id)))
  total   <- length(unique(per_sample_tests$gene_id))
  message(sprintf("Functional annotations available for %d / %d unique genes (%.1f%%).",
                  matched, total, ifelse(total > 0, 100 * matched / total, NaN)))
} else {
  message("Functional annotations unavailable for this run; gene-level tables will omit annotation fields.")
}

# ---------------------------------
# 7) Write tidy master tables for downstream work
# ---------------------------------
write_tsv_smart(annotate_if_available(dat, anno_slim),             file.path(results_dir, "ase_counts_long.tsv.gz"))
write_tsv_smart(annotate_if_available(per_sample_tests, anno_slim), file.path(results_dir, "ase_per_sample_binomial.tsv.gz"))
write_tsv_smart(annotate_if_available(group_results, anno_slim),   file.path(results_dir, "ase_group_results.tsv.gz"))
if (exists("diffase")) write_tsv_smart(annotate_if_available(diffase, anno_slim), file.path(results_dir, "ase_diff_between_groups.tsv.gz"))

message("Done. Outputs written next to ASE counts: \n  results: ", results_dir, "\n  plots:   ", plots_dir)
