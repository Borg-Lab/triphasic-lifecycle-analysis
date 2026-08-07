required_path_keys <- c(
  "counts_tsv", "tpm_tsv", "annotation_tsv", "gff_gtf",
  "gene_age_tsv", "orthogroup_tsv", "output_root"
)

required_analysis_keys <- c(
  "min_rowmax_counts", "tpm_detect_threshold", "marker_tau_threshold",
  "marker_min_log2fc", "marker_max_fdr", "wgcna_enable",
  "chrom_cluster_enable", "chrom_perm_n", "chrom_seed",
  "chrom_distance_metric"
)

required_design_cols <- c(
  "sample_id", "generation", "stage", "tissue", "sex", "replicate", "batch", "include"
)

required_contrast_cols <- c(
  "contrast_id", "factor", "numerator", "denominator", "lfc_threshold", "fdr_threshold"
)

require_columns <- function(df, cols, name) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(sprintf("Missing required columns in %s: %s", name, paste(miss, collapse = ", ")), call. = FALSE)
  }
}

read_yaml_checked <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with install.packages('yaml').", call. = FALSE)
  }
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path), call. = FALSE)
  yaml::read_yaml(path)
}

check_keys <- function(x, keys, name) {
  miss <- setdiff(keys, names(x))
  if (length(miss) > 0) {
    stop(sprintf("Missing required keys in %s: %s", name, paste(miss, collapse = ", ")), call. = FALSE)
  }
}

read_paths_config <- function(path = "config/paths.yml") {
  cfg <- read_yaml_checked(path)
  check_keys(cfg, required_path_keys, "config/paths.yml")
  cfg
}

read_analysis_config <- function(path = "config/analysis.yml") {
  cfg <- read_yaml_checked(path)
  check_keys(cfg, required_analysis_keys, "config/analysis.yml")
  cfg
}

read_design_table <- function(path = "config/design.tsv", include_all = FALSE) {
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path), call. = FALSE)
  d <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(d, required_design_cols, "config/design.tsv")
  if (!isTRUE(include_all)) {
    d <- d[as.logical(d$include), , drop = FALSE]
  }
  d
}

read_contrast_table <- function(path = "config/contrasts.tsv") {
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path), call. = FALSE)
  ctb <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  require_columns(ctb, required_contrast_cols, "config/contrasts.tsv")
  ctb
}

validate_configs <- function(paths_cfg, analysis_cfg, design_df, contrast_df) {
  check_keys(paths_cfg, required_path_keys, "paths config")
  check_keys(analysis_cfg, required_analysis_keys, "analysis config")
  require_columns(design_df, required_design_cols, "design table")
  require_columns(contrast_df, required_contrast_cols, "contrast table")
  TRUE
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_read_table <- function(path, id_col = NULL) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(data.frame())
  x <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(id_col) && id_col %in% colnames(x)) {
    x[[id_col]] <- as.character(x[[id_col]])
  }
  x
}

normalize_gene_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

.protein_coding_cache <- new.env(parent = emptyenv())

extract_gene_ids_from_gr <- function(gr) {
  meta <- S4Vectors::mcols(gr)
  gid <- rep(NA_character_, length(gr))
  candidates <- c("gene_id", "ID", "Name", "gene", "locus_tag", "Parent", "transcript_id")
  for (nm in candidates) {
    if (!nm %in% colnames(meta)) next
    vals <- as.character(meta[[nm]])
    vals <- sub(",.*$", "", vals)
    vals <- sub("^[^:]+:", "", vals)
    vals <- sub("\\.t[0-9]+$", "", vals)
    fill <- (is.na(gid) | !nzchar(gid)) & !is.na(vals) & nzchar(vals)
    gid[fill] <- vals[fill]
  }
  normalize_gene_id(gid)
}

get_protein_coding_gene_ids <- function(paths_cfg) {
  gff <- paths_cfg$gff_gtf
  if (is.null(gff) || !nzchar(gff) || !file.exists(gff)) {
    stop("Protein-coding filtering requires a valid gff_gtf path in config/paths.yml.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required for protein-coding filtering from GFF/GTF.", call. = FALSE)
  }
  if (!exists(gff, envir = .protein_coding_cache, inherits = FALSE)) {
    gr <- rtracklayer::import(gff)
    meta <- S4Vectors::mcols(gr)
    type_vec <- if ("type" %in% colnames(meta)) as.character(meta[["type"]]) else rep(NA_character_, length(gr))
    idx <- which(type_vec == "gene")
    if (length(idx) == 0) {
      idx <- which(type_vec %in% c("mRNA", "transcript"))
    }
    if (length(idx) == 0) {
      stop("Could not find gene/transcript features in gff_gtf for protein-coding filtering.", call. = FALSE)
    }
    g <- gr[idx]
    gmeta <- S4Vectors::mcols(g)
    biotype_cols <- c("biotype", "gene_biotype", "gene_type", "transcript_biotype", "transcript_type")
    bt_col <- biotype_cols[biotype_cols %in% colnames(gmeta)][1]
    if (is.na(bt_col) || !nzchar(bt_col)) {
      stop("No biotype column found in gff_gtf (expected one of: biotype, gene_biotype, gene_type, transcript_biotype, transcript_type).", call. = FALSE)
    }
    bt <- tolower(as.character(gmeta[[bt_col]]))
    keep <- !is.na(bt) & bt == "protein_coding"
    gid <- extract_gene_ids_from_gr(g)
    gids <- unique(gid[keep & !is.na(gid) & nzchar(gid)])
    if (length(gids) == 0) {
      stop("Protein-coding filtering found zero genes. Check gff_gtf biotype annotations.", call. = FALSE)
    }
    assign(gff, gids, envir = .protein_coding_cache)
  }
  get(gff, envir = .protein_coding_cache, inherits = FALSE)
}

read_functional_annotation <- function(paths_cfg) {
  anno <- safe_read_table(paths_cfg$annotation_tsv)
  if (nrow(anno) == 0) return(data.frame(gene_id = character(), stringsAsFactors = FALSE))

  id_col <- NULL
  for (cand in c("gene_id", "query", "transcript_id", "ID", "Name")) {
    if (cand %in% colnames(anno)) {
      id_col <- cand
      break
    }
  }
  if (is.null(id_col)) return(data.frame(gene_id = character(), stringsAsFactors = FALSE))

  anno$gene_id <- normalize_gene_id(anno[[id_col]])
  anno <- anno[!is.na(anno$gene_id) & nzchar(anno$gene_id), , drop = FALSE]
  pc_genes <- get_protein_coding_gene_ids(paths_cfg)
  anno <- anno[anno$gene_id %in% pc_genes, , drop = FALSE]
  anno <- anno[!duplicated(anno$gene_id), , drop = FALSE]
  rownames(anno) <- NULL
  anno
}

write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  path
}

read_coordinates <- function(paths_cfg) {
  coord_cols <- c("gene_id", "chr", "start", "end", "strand")
  anno <- safe_read_table(paths_cfg$annotation_tsv)
  if (nrow(anno) > 0 && all(coord_cols %in% colnames(anno))) {
    out <- anno[, coord_cols]
  } else if (!is.null(paths_cfg$gff_gtf) && nzchar(paths_cfg$gff_gtf) && file.exists(paths_cfg$gff_gtf) && requireNamespace("rtracklayer", quietly = TRUE)) {
    gr <- rtracklayer::import(paths_cfg$gff_gtf)
    meta <- S4Vectors::mcols(gr)
    type_vec <- if ("type" %in% colnames(meta)) as.character(meta[["type"]]) else rep(NA_character_, length(gr))
    idx <- which(type_vec == "gene")
    if (length(idx) == 0) {
      idx <- which(type_vec %in% c("mRNA", "transcript"))
    }
    g <- gr[idx]
    gid <- extract_gene_ids_from_gr(g)

    out <- data.frame(
      gene_id = gid,
      chr = as.character(GenomeInfoDb::seqnames(g)),
      start = as.integer(IRanges::start(g)),
      end = as.integer(IRanges::end(g)),
      strand = as.character(BiocGenerics::strand(g)),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(gene_id = character(), chr = character(), start = integer(), end = integer(), strand = character())
  }

  out <- out[!is.na(out$gene_id) & nzchar(out$gene_id), ]
  out <- out[!is.na(out$chr) & !is.na(out$start) & !is.na(out$end), ]
  out <- out[out$end >= out$start, ]
  pc_genes <- get_protein_coding_gene_ids(paths_cfg)
  out <- out[out$gene_id %in% pc_genes, , drop = FALSE]
  out <- out[!duplicated(out$gene_id), ]
  rownames(out) <- NULL
  out
}

filter_coordinates <- function(coords_df, analysis_cfg) {
  if (nrow(coords_df) == 0) return(coords_df)
  exclude_chr <- analysis_cfg$exclude_chromosomes
  if (is.null(exclude_chr) || length(exclude_chr) == 0) {
    exclude_chr <- c("ChrPt", "ChrMt")
  }
  exclude_chr <- as.character(unlist(exclude_chr))
  coords_df <- coords_df[!(coords_df$chr %in% exclude_chr), , drop = FALSE]
  rownames(coords_df) <- NULL
  coords_df
}

prepare_expression_inputs <- function(paths_cfg, design_df, analysis_cfg) {
  counts <- safe_read_table(paths_cfg$counts_tsv)
  tpm <- safe_read_table(paths_cfg$tpm_tsv)
  pc_genes <- get_protein_coding_gene_ids(paths_cfg)

  if (nrow(counts) > 0) {
    id_col <- colnames(counts)[1]
    colnames(counts)[1] <- "gene_id"
    if (id_col != "gene_id") counts$gene_id <- as.character(counts$gene_id)
    counts$gene_id <- normalize_gene_id(counts$gene_id)
    counts <- counts[!is.na(counts$gene_id) & nzchar(counts$gene_id), , drop = FALSE]
    counts <- counts[counts$gene_id %in% pc_genes, , drop = FALSE]
    counts <- counts[!duplicated(counts$gene_id), , drop = FALSE]
  } else {
    counts <- data.frame(gene_id = character())
  }

  if (nrow(tpm) > 0) {
    id_col <- colnames(tpm)[1]
    colnames(tpm)[1] <- "gene_id"
    if (id_col != "gene_id") tpm$gene_id <- as.character(tpm$gene_id)
    tpm$gene_id <- normalize_gene_id(tpm$gene_id)
    tpm <- tpm[!is.na(tpm$gene_id) & nzchar(tpm$gene_id), , drop = FALSE]
    tpm <- tpm[tpm$gene_id %in% pc_genes, , drop = FALSE]
    tpm <- tpm[!duplicated(tpm$gene_id), , drop = FALSE]
  } else {
    tpm <- data.frame(gene_id = character())
  }

  keep_samples <- design_df$sample_id[as.logical(design_df$include)]
  counts_cols <- intersect(keep_samples, colnames(counts))
  tpm_cols <- intersect(keep_samples, colnames(tpm))

  if (length(counts_cols) > 0) {
    counts <- counts[, c("gene_id", counts_cols), drop = FALSE]
    mat <- as.matrix(counts[, -1, drop = FALSE])
    storage.mode(mat) <- "numeric"
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required for count filtering via filterByExpr.", call. = FALSE)
    }
    group <- as.factor(as.character(design_df$stage[match(colnames(mat), design_df$sample_id)]))
    if (any(is.na(group))) {
      stop("Some count columns are missing stage assignments in design.tsv.", call. = FALSE)
    }
    dge <- edgeR::DGEList(counts = mat, group = group)
    keep <- edgeR::filterByExpr(dge, group = group)
    if (!any(keep)) {
      stop("All genes were filtered by edgeR::filterByExpr; check counts/design settings.", call. = FALSE)
    }
    counts <- counts[keep, , drop = FALSE]
  }

  if (length(tpm_cols) > 0) {
    tpm <- tpm[, c("gene_id", tpm_cols), drop = FALSE]
  }

  list(counts = counts, tpm = tpm)
}
