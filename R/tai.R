normalize_gene_id <- function(x) {
  x <- as.character(x)
  x <- sub(",.*$", "", x)
  x <- sub("^[^:]+:", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

order_stage_columns <- function(stage_cols) {
  preferred_suffix <- c(
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
  preferred <- paste0("mean_tpm_", preferred_suffix)
  ordered <- intersect(preferred, stage_cols)
  c(ordered, setdiff(stage_cols, ordered))
}

read_gene_age_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(data.frame())

  age <- tryCatch(
    read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  has_expected_header <- any(c("gene_id", "ID", "id", "rank", "phylostratum", "age", "age_rank") %in% colnames(age))
  if (!has_expected_header) {
    age <- tryCatch(
      read.delim(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) data.frame()
    )
    if (ncol(age) >= 1) colnames(age)[1] <- "gene_id"
    if (ncol(age) >= 2) colnames(age)[2] <- "taxon"
    if (ncol(age) >= 3) colnames(age)[3] <- "rank"
    if (ncol(age) >= 4) colnames(age)[4] <- "support"
  }

  age
}

prepare_mytai_input <- function(group_tpm_df, age_df) {
  if (nrow(group_tpm_df) == 0 || nrow(age_df) == 0) return(data.frame())

  age_id_col <- NULL
  for (cand in c("gene_id", "ID", "id")) {
    if (cand %in% colnames(age_df)) {
      age_id_col <- cand
      break
    }
  }
  if (is.null(age_id_col)) {
    age_id_col <- colnames(age_df)[1]
  }

  age_col <- NULL
  for (cand in c("rank", "phylostratum", "age", "age_rank")) {
    if (cand %in% colnames(age_df)) {
      age_col <- cand
      break
    }
  }
  if (is.null(age_col)) return(data.frame())

  age_df$gene_id <- normalize_gene_id(age_df[[age_id_col]])
  age_df$rank_num <- suppressWarnings(as.numeric(as.character(age_df[[age_col]])))
  age_df <- age_df[!is.na(age_df$gene_id) & nzchar(age_df$gene_id) & is.finite(age_df$rank_num), c("gene_id", "rank_num"), drop = FALSE]
  age_df <- age_df[!duplicated(age_df$gene_id), , drop = FALSE]

  group_tpm_df$gene_id <- normalize_gene_id(group_tpm_df$gene_id)
  merged <- merge(age_df, group_tpm_df, by = "gene_id", all = FALSE)
  if (nrow(merged) == 0) return(data.frame())

  stage_cols_all <- setdiff(colnames(group_tpm_df), "gene_id")
  stage_keep <- vapply(stage_cols_all, function(cc) {
    x <- suppressWarnings(as.numeric(merged[[cc]]))
    any(is.finite(x) & x > 0)
  }, logical(1))
  stage_cols <- order_stage_columns(stage_cols_all[stage_keep])
  if (length(stage_cols) == 0) return(data.frame())

  for (cc in stage_cols) {
    x <- suppressWarnings(as.numeric(merged[[cc]]))
    x[!is.finite(x)] <- 0
    merged[[cc]] <- x
  }

  out <- data.frame(
    rank = merged$rank_num,
    ID = merged$gene_id,
    merged[, stage_cols, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out
}

build_tai_context <- function(group_tpm_df, paths_cfg) {
  schema <- list(
    phyex_df = data.frame(),
    phyex_set = NULL,
    stage_cols = character(),
    expr_mat = matrix(numeric(0), nrow = 0, ncol = 0)
  )
  if (nrow(group_tpm_df) == 0) return(schema)
  if (!requireNamespace("myTAI", quietly = TRUE)) {
    stop("Package 'myTAI' is required for TAI analysis but is not installed.", call. = FALSE)
  }

  age <- read_gene_age_table(paths_cfg$gene_age_tsv)
  if (nrow(age) == 0) return(schema)

  phyex_df <- prepare_mytai_input(group_tpm_df, age)
  if (nrow(phyex_df) == 0) return(schema)

  phyex_set <- myTAI::as_BulkPhyloExpressionSet(phyex_df, name = "triphasic_phyex_set")
  stage_cols <- setdiff(colnames(phyex_df), c("rank", "ID"))
  expr_mat <- as.matrix(phyex_df[, stage_cols, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"

  list(
    phyex_df = phyex_df,
    phyex_set = phyex_set,
    stage_cols = stage_cols,
    expr_mat = expr_mat
  )
}

run_tai_analysis <- function(group_tpm_df, paths_cfg) {
  schema <- data.frame(stage = character(), TAI = numeric(), n_genes = integer(), stringsAsFactors = FALSE)
  ctx <- build_tai_context(group_tpm_df, paths_cfg)
  if (is.null(ctx$phyex_set) || length(ctx$stage_cols) == 0) return(schema)

  tai_vals <- myTAI::TAI(ctx$phyex_set)
  n_genes <- colSums(is.finite(ctx$expr_mat) & ctx$expr_mat > 0, na.rm = TRUE)

  tai_num <- as.numeric(tai_vals)
  names(tai_num) <- names(tai_vals)
  out <- data.frame(
    stage = ctx$stage_cols,
    TAI = tai_num[ctx$stage_cols],
    n_genes = as.integer(n_genes[ctx$stage_cols]),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$TAI) & out$n_genes > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

safe_pdf_plot <- function(path, plotting_fn) {
  grDevices::pdf(path, width = 8, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plotting_fn()
  if (inherits(p, "ggplot")) {
    print(p)
  } else if (inherits(p, "grob") || inherits(p, "gTree")) {
    if (requireNamespace("grid", quietly = TRUE)) grid::grid.draw(p)
  }
  path
}


transformed_expression_matrix <- function(x, stage_cols) {
  expr <- tryCatch(attributes(x)$.expression, error = function(e) NULL)
  if (is.null(expr)) expr <- tryCatch(x$expression, error = function(e) NULL)
  if (is.null(expr)) return(NULL)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  keep_cols <- intersect(stage_cols, colnames(expr))
  if (length(keep_cols) == 0) return(NULL)
  expr[, keep_cols, drop = FALSE]
}

valid_tai_expression_weights <- function(x, stage_cols) {
  expr <- transformed_expression_matrix(x, stage_cols)
  if (is.null(expr)) return(FALSE)
  all(is.finite(expr)) && all(expr >= 0)
}

plot_transformed_tai_table <- function(trans_tbl) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot.new()
    text(0.5, 0.5, "ggplot2 unavailable")
    return(invisible(NULL))
  }
  ggplot2::ggplot(trans_tbl, ggplot2::aes(x = stage, y = value, color = transformation, group = transformation)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = "TAI", color = "Transformation")
}

run_tai_diagnostics <- function(group_tpm_df, paths_cfg, output_dir) {
  ensure_dir(output_dir)
  out_files <- character()
  ctx <- build_tai_context(group_tpm_df, paths_cfg)
  if (is.null(ctx$phyex_set) || length(ctx$stage_cols) == 0) return(out_files)

  phyex_set <- ctx$phyex_set
  stage_cols <- ctx$stage_cols

  # Export plotting-ready numeric tables for all TAI curve variants.
  # 1) Raw signature
  raw_tai <- myTAI::TAI(phyex_set)
  raw_tbl <- data.frame(
    curve = "raw",
    stage = stage_cols,
    value = as.numeric(raw_tai[stage_cols]),
    stringsAsFactors = FALSE
  )
  raw_file <- file.path(output_dir, "tai_signature_values.tsv")
  write_tsv(raw_tbl, raw_file)
  out_files <- c(out_files, raw_file)

  # 2) Count-transformed signatures (same variants used by plot_signature_transformed)
  trans_tbl <- data.frame(
    transformation = character(),
    stage = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  transforms <- NULL
  transforms <- tryCatch(myTAI:::COUNT_TRANSFORMS, error = function(e) NULL)
  if (!is.null(transforms) && length(transforms) > 0) {
    for (nm in names(transforms)) {
      set_tr <- tryCatch(
        myTAI::transform_counts(phyex_set, FUN = transforms[[nm]], FUN_name = nm),
        error = function(e) NULL
      )
      if (is.null(set_tr)) next
      if (!valid_tai_expression_weights(set_tr, stage_cols)) {
        warning(sprintf("Skipping TAI transformation '%s' because it produced negative or non-finite expression weights.", nm), call. = FALSE)
        next
      }
      tvals <- tryCatch(myTAI::TAI(set_tr), error = function(e) NULL)
      if (is.null(tvals)) next
      trans_tbl <- rbind(
        trans_tbl,
        data.frame(
          transformation = nm,
          stage = stage_cols,
          value = as.numeric(tvals[stage_cols]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  if (nrow(trans_tbl) > 0) {
    trans_file <- file.path(output_dir, "tai_signature_transformed_values.tsv")
    write_tsv(trans_tbl, trans_file)
    out_files <- c(out_files, trans_file)
  }

  # 3) Rank-transformed signature (same variant used in tai_signature_rank_transformed.pdf)
  rank_tbl <- data.frame(
    curve = character(),
    stage = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  phyex_rank <- tryCatch(
    myTAI::tf(phyex_set, FUN = function(x) apply(x, 2, base::rank)),
    error = function(e) NULL
  )
  if (!is.null(phyex_rank)) {
    rvals <- tryCatch(myTAI::TAI(phyex_rank), error = function(e) NULL)
    if (!is.null(rvals)) {
      rank_tbl <- data.frame(
        curve = "rank_transformed",
        stage = stage_cols,
        value = as.numeric(rvals[stage_cols]),
        stringsAsFactors = FALSE
      )
      rank_file <- file.path(output_dir, "tai_signature_rank_transformed_values.tsv")
      write_tsv(rank_tbl, rank_file)
      out_files <- c(out_files, rank_file)
    }
  }

  # 4) Gene-quantile signatures (same variants used by plot_signature_gene_quantiles)
  q_tbl <- data.frame(
    quantile = numeric(),
    stage = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  quantiles <- c(1, 0.99, 0.95, 0.9, 0.8)
  for (q in quantiles) {
    genes_q <- tryCatch(myTAI::genes_top_mean(phyex_set, p = q), error = function(e) NULL)
    if (is.null(genes_q)) next
    set_q <- tryCatch(myTAI::remove_genes(phyex_set, genes = genes_q), error = function(e) NULL)
    if (is.null(set_q)) next
    qvals <- tryCatch(myTAI::TAI(set_q), error = function(e) NULL)
    if (is.null(qvals)) next
    q_tbl <- rbind(
      q_tbl,
      data.frame(
        quantile = q,
        stage = stage_cols,
        value = as.numeric(qvals[stage_cols]),
        stringsAsFactors = FALSE
      )
    )
  }
  if (nrow(q_tbl) > 0) {
    q_file <- file.path(output_dir, "tai_signature_gene_quantiles_values.tsv")
    write_tsv(q_tbl, q_file)
    out_files <- c(out_files, q_file)
  }

  flatline_ok <- TRUE
  flt <- tryCatch(
    myTAI::stat_flatline_test(phyex_set, plot_result = FALSE),
    error = function(e) {
      flatline_ok <<- FALSE
      NULL
    }
  )

  if (flatline_ok && !is.null(flt)) {
    test_tbl <- data.frame(
      method = "flatline_test",
      method_name = flt@method_name,
      test_stat = as.numeric(flt@test_stat),
      p_label = flt@p_label,
      p_value = as.numeric(flt@p_value),
      alternative = flt@alternative,
      null_n = length(flt@null_sample),
      stringsAsFactors = FALSE
    )
    test_file <- file.path(output_dir, "tai_flatline_test.tsv")
    write_tsv(test_tbl, test_file)
    out_files <- c(out_files, test_file)
  }

  sig_file <- file.path(output_dir, "tai_signature.pdf")
  safe_pdf_plot(sig_file, function() myTAI::plot_signature(phyex_set, show_p_val = TRUE))
  out_files <- c(out_files, sig_file)

  sig_tf_file <- file.path(output_dir, "tai_signature_transformed.pdf")
  if (nrow(trans_tbl) > 0) {
    safe_pdf_plot(sig_tf_file, function() plot_transformed_tai_table(trans_tbl))
    out_files <- c(out_files, sig_tf_file)
  }

  sig_rank_file <- file.path(output_dir, "tai_signature_rank_transformed.pdf")
  safe_pdf_plot(sig_rank_file, function() {
    phyex_rank_plot <- myTAI::tf(phyex_set, FUN = function(x) apply(x, 2, base::rank))
    myTAI::plot_signature(phyex_rank_plot)
  })
  out_files <- c(out_files, sig_rank_file)

  sig_q_file <- file.path(output_dir, "tai_signature_gene_quantiles.pdf")
  safe_pdf_plot(sig_q_file, function() myTAI::plot_signature_gene_quantiles(phyex_set))
  out_files <- c(out_files, sig_q_file)

  if (flatline_ok && !is.null(flt)) {
    cf_file <- file.path(output_dir, "tai_flatline_cullen_frey.pdf")
    safe_pdf_plot(cf_file, function() myTAI::plot_cullen_frey(flt))
    out_files <- c(out_files, cf_file)

    null_file <- file.path(output_dir, "tai_flatline_null_txi_sample.pdf")
    safe_pdf_plot(null_file, function() myTAI::plot_null_txi_sample(flt))
    out_files <- c(out_files, null_file)
  }

  out_files[file.exists(out_files)]
}
