compute_midpoint <- function(start, end) {
  (as.numeric(start) + as.numeric(end)) / 2
}

identify_proximity_blocks <- function(set_id, gene_ids, coords_df, gap_bp = 50000L, min_block_genes = 2L) {
  sub <- coords_df[coords_df$gene_id %in% gene_ids, , drop = FALSE]
  if (nrow(sub) == 0) {
    return(data.frame(
      set_id = character(),
      chr = character(),
      block_id = integer(),
      block_start = integer(),
      block_end = integer(),
      span_bp = numeric(),
      n_genes = integer(),
      gene_ids = character(),
      stringsAsFactors = FALSE
    ))
  }

  sub <- sub[order(sub$chr, sub$midpoint), , drop = FALSE]
  out <- list()
  out_idx <- 0L
  block_id <- 0L

  for (ch in unique(sub$chr)) {
    chr_df <- sub[sub$chr == ch, , drop = FALSE]
    if (nrow(chr_df) < min_block_genes) next

    current_idx <- 1L
    for (i in 2:nrow(chr_df)) {
      gap <- chr_df$midpoint[i] - chr_df$midpoint[i - 1]
      if (is.na(gap) || gap > gap_bp) {
        block_idx <- current_idx:(i - 1L)
        if (length(block_idx) >= min_block_genes) {
          block_id <- block_id + 1L
          block_start <- min(chr_df$start[block_idx], na.rm = TRUE)
          block_end <- max(chr_df$end[block_idx], na.rm = TRUE)
          span_bp <- block_end - block_start + 1
          out_idx <- out_idx + 1L
          out[[out_idx]] <- data.frame(
            set_id = set_id,
            chr = ch,
            block_id = block_id,
            block_start = as.integer(block_start),
            block_end = as.integer(block_end),
            span_bp = as.numeric(span_bp),
            n_genes = as.integer(length(block_idx)),
            gene_ids = paste(chr_df$gene_id[block_idx], collapse = ";"),
            stringsAsFactors = FALSE
          )
        }
        current_idx <- i
      }
    }

    # Final block on this chromosome
    block_idx <- current_idx:nrow(chr_df)
    if (length(block_idx) >= min_block_genes) {
      block_id <- block_id + 1L
      block_start <- min(chr_df$start[block_idx], na.rm = TRUE)
      block_end <- max(chr_df$end[block_idx], na.rm = TRUE)
      span_bp <- block_end - block_start + 1
      out_idx <- out_idx + 1L
      out[[out_idx]] <- data.frame(
        set_id = set_id,
        chr = ch,
        block_id = block_id,
        block_start = as.integer(block_start),
        block_end = as.integer(block_end),
        span_bp = as.numeric(span_bp),
        n_genes = as.integer(length(block_idx)),
        gene_ids = paste(chr_df$gene_id[block_idx], collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) {
    return(data.frame(
      set_id = character(),
      chr = character(),
      block_id = integer(),
      block_start = integer(),
      block_end = integer(),
      span_bp = numeric(),
      n_genes = integer(),
      gene_ids = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

summarize_proximity_blocks <- function(block_df, total_set_genes) {
  if (nrow(block_df) == 0) {
    return(
      list(
        n_blocks = 0L,
        n_genes_in_blocks = 0L,
        pct_genes_in_proximity_blocks = 0,
        avg_block_span_bp = NA_real_,
        median_block_span_bp = NA_real_
      )
    )
  }

  block_spans <- block_df$span_bp
  block_genes <- unique(unlist(strsplit(block_df$gene_ids, ";", fixed = TRUE)))
  n_blocks <- nrow(block_df)
  n_genes_in_blocks <- length(unique(block_genes))
  pct <- if (total_set_genes > 0) 100 * n_genes_in_blocks / total_set_genes else 0

  list(
    n_blocks = as.integer(n_blocks),
    n_genes_in_blocks = as.integer(n_genes_in_blocks),
    pct_genes_in_proximity_blocks = pct,
    avg_block_span_bp = if (n_blocks > 0) mean(block_spans) else NA_real_,
    median_block_span_bp = if (n_blocks > 0) stats::median(block_spans) else NA_real_
  )
}

nearest_neighbor_distances <- function(coords_df) {
  if (nrow(coords_df) < 2) return(numeric())
  coords_df <- coords_df[order(coords_df$chr, coords_df$midpoint), , drop = FALSE]

  nn <- numeric(0)
  for (ch in unique(coords_df$chr)) {
    x <- coords_df[coords_df$chr == ch, "midpoint"]
    x <- sort(x)
    if (length(x) < 2) next
    d <- sapply(seq_along(x), function(i) {
      others <- x[-i]
      min(abs(x[i] - others))
    })
    nn <- c(nn, d)
  }
  nn
}

compute_set_stat <- function(gene_ids, coords_df) {
  sub <- coords_df[coords_df$gene_id %in% gene_ids, , drop = FALSE]
  if (nrow(sub) < 2) return(NA_real_)
  nn <- nearest_neighbor_distances(sub)
  if (length(nn) == 0) return(NA_real_)
  stats::median(nn, na.rm = TRUE)
}

permute_set_stat <- function(gene_ids, coords_df, n_perm = 1000, seed = 42) {
  set.seed(seed)
  sub <- coords_df[coords_df$gene_id %in% gene_ids, , drop = FALSE]
  if (nrow(sub) < 2) return(rep(NA_real_, n_perm))

  chr_counts <- table(sub$chr)
  pool_by_chr <- split(coords_df$gene_id, coords_df$chr)

  stats_vec <- rep(NA_real_, n_perm)
  for (i in seq_len(n_perm)) {
    pick <- character(0)
    valid <- TRUE
    for (ch in names(chr_counts)) {
      n <- as.integer(chr_counts[[ch]])
      pool <- pool_by_chr[[ch]]
      if (length(pool) < n) {
        valid <- FALSE
        break
      }
      pick <- c(pick, sample(pool, size = n, replace = FALSE))
    }
    if (!valid) next
    stats_vec[i] <- compute_set_stat(pick, coords_df)
  }
  stats_vec
}

run_chromosomal_clustering <- function(module_df, coords_df, analysis_cfg, extra_sets = NULL) {
  summary_schema <- data.frame(
    set_id = character(), chr_scope = character(), n_genes = integer(),
    obs_median_nn_dist = numeric(), perm_median = numeric(), z_score = numeric(),
    empirical_p = numeric(), fdr = numeric(), cluster_flag = logical(),
    n_proximity_blocks = integer(), n_genes_in_proximity_blocks = integer(),
    pct_genes_in_proximity_blocks = numeric(),
    avg_block_span_bp = numeric(), median_block_span_bp = numeric(),
    stringsAsFactors = FALSE
  )

  members_schema <- data.frame(set_id = character(), gene_id = character(), stringsAsFactors = FALSE)
  blocks_schema <- data.frame(
    set_id = character(), chr = character(), block_id = integer(),
    block_start = integer(), block_end = integer(), span_bp = numeric(),
    n_genes = integer(), gene_ids = character(), stringsAsFactors = FALSE
  )
  nulls <- list()

  if (!isTRUE(analysis_cfg$chrom_cluster_enable)) {
    return(list(summary = summary_schema, members = members_schema, blocks = blocks_schema, nulls = nulls))
  }

  if (nrow(coords_df) == 0) {
    return(list(summary = summary_schema, members = members_schema, blocks = blocks_schema, nulls = nulls))
  }

  coords_df$midpoint <- compute_midpoint(coords_df$start, coords_df$end)
  coords_df <- coords_df[!is.na(coords_df$midpoint), , drop = FALSE]

  sets <- list()
  if (nrow(module_df) > 0 && all(c("gene_id", "module_id") %in% colnames(module_df))) {
    sets <- split(module_df$gene_id, module_df$module_id)
  }
  if (!is.null(extra_sets) && length(extra_sets) > 0) {
    sets <- c(sets, extra_sets)
  }
  if (length(sets) == 0) {
    return(list(summary = summary_schema, members = members_schema, blocks = blocks_schema, nulls = nulls))
  }

  n_perm <- as.integer(analysis_cfg$chrom_perm_n)
  seed <- as.integer(analysis_cfg$chrom_seed)
  min_set <- if (!is.null(analysis_cfg$chrom_min_set_size)) as.integer(analysis_cfg$chrom_min_set_size) else 5L
  alpha <- if (!is.null(analysis_cfg$alpha)) as.numeric(analysis_cfg$alpha) else 0.05
  block_gap_bp <- if (!is.null(analysis_cfg$chrom_block_gap_bp)) as.integer(analysis_cfg$chrom_block_gap_bp) else 50000L
  block_min_genes <- if (!is.null(analysis_cfg$chrom_block_min_genes)) as.integer(analysis_cfg$chrom_block_min_genes) else 2L

  rows <- list()
  mrows <- list()
  brows <- list()
  idx <- 0
  for (sid in names(sets)) {
    genes <- unique(sets[[sid]])
    genes <- intersect(genes, coords_df$gene_id)
    if (length(genes) < min_set) next

    obs <- compute_set_stat(genes, coords_df)
    perm <- permute_set_stat(genes, coords_df, n_perm = n_perm, seed = seed + idx)
    perm_ok <- perm[!is.na(perm)]
    if (is.na(obs) || length(perm_ok) == 0) next

    emp <- (sum(perm_ok <= obs) + 1) / (length(perm_ok) + 1)
    z <- ifelse(stats::sd(perm_ok) > 0, (obs - mean(perm_ok)) / stats::sd(perm_ok), NA_real_)
    chr_scope <- paste(sort(unique(coords_df$chr[coords_df$gene_id %in% genes])), collapse = ",")
    block_df <- identify_proximity_blocks(
      set_id = sid,
      gene_ids = genes,
      coords_df = coords_df,
      gap_bp = block_gap_bp,
      min_block_genes = block_min_genes
    )
    block_stats <- summarize_proximity_blocks(block_df, total_set_genes = length(genes))

    idx <- idx + 1
    rows[[idx]] <- data.frame(
      set_id = sid,
      chr_scope = chr_scope,
      n_genes = length(genes),
      obs_median_nn_dist = obs,
      perm_median = stats::median(perm_ok),
      z_score = z,
      empirical_p = emp,
      n_proximity_blocks = block_stats$n_blocks,
      n_genes_in_proximity_blocks = block_stats$n_genes_in_blocks,
      pct_genes_in_proximity_blocks = block_stats$pct_genes_in_proximity_blocks,
      avg_block_span_bp = block_stats$avg_block_span_bp,
      median_block_span_bp = block_stats$median_block_span_bp,
      stringsAsFactors = FALSE
    )
    mrows[[idx]] <- data.frame(set_id = sid, gene_id = genes, stringsAsFactors = FALSE)
    brows[[idx]] <- block_df
    nulls[[sid]] <- perm_ok
  }

  if (length(rows) == 0) {
    return(list(summary = summary_schema, members = members_schema, blocks = blocks_schema, nulls = nulls))
  }

  summary_df <- do.call(rbind, rows)
  summary_df$fdr <- stats::p.adjust(summary_df$empirical_p, method = "BH")
  summary_df$cluster_flag <- with(summary_df, fdr <= alpha & z_score < 0)

  list(
    summary = summary_df,
    members = do.call(rbind, mrows),
    blocks = do.call(rbind, brows),
    nulls = nulls
  )
}
