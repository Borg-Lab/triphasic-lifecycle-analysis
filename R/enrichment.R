run_marker_enrichment <- function(cluster_summary_df, cluster_members_df, marker_df, universe_genes) {
  schema <- data.frame(
    cluster_id = character(), marker_group = character(), overlap = integer(),
    overlap_gene_ids = character(),
    odds_ratio = numeric(), p_value = numeric(), fdr = numeric(), test_method = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(cluster_summary_df) == 0 || nrow(cluster_members_df) == 0 || nrow(marker_df) == 0) return(schema)

  clustered <- cluster_summary_df$set_id[cluster_summary_df$cluster_flag]
  if (length(clustered) == 0) return(schema)

  marker_df <- marker_df[isTRUE(marker_df$marker_flag) | marker_df$marker_flag == TRUE, , drop = FALSE]
  if (nrow(marker_df) == 0) return(schema)

  marker_sets <- split(marker_df$gene_id, marker_df$marker_group)
  universe <- unique(universe_genes)

  rows <- list()
  idx <- 0
  for (cid in clustered) {
    cgenes <- unique(cluster_members_df$gene_id[cluster_members_df$set_id == cid])
    cgenes <- intersect(cgenes, universe)
    for (mg in names(marker_sets)) {
      mgenes <- intersect(unique(marker_sets[[mg]]), universe)
      overlap_genes <- sort(intersect(cgenes, mgenes))
      a <- length(overlap_genes)
      b <- length(setdiff(cgenes, mgenes))
      c <- length(setdiff(mgenes, cgenes))
      d <- length(setdiff(universe, union(cgenes, mgenes)))
      tab <- matrix(c(a, b, c, d), nrow = 2)
      ft <- stats::fisher.test(tab, alternative = "greater")

      idx <- idx + 1
      rows[[idx]] <- data.frame(
        cluster_id = cid,
        marker_group = mg,
        overlap = a,
        overlap_gene_ids = paste(overlap_genes, collapse = ","),
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value,
        test_method = "fisher_greater",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) return(schema)
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[, c("cluster_id", "marker_group", "overlap", "odds_ratio", "p_value", "fdr", "test_method", "overlap_gene_ids")]
  out
}

split_annotation_terms <- function(x, pattern = ",") {
  if (length(x) == 0) return(character())
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "-" & x != "NA"]
  if (length(x) == 0) return(character())
  terms <- unlist(strsplit(x, split = pattern))
  terms <- trimws(terms)
  terms <- terms[nzchar(terms) & terms != "-" & terms != "NA"]
  unique(terms)
}

build_term_gene_map <- function(marker_df, term_col, split_pattern = ",") {
  if (!term_col %in% colnames(marker_df)) return(list())
  out <- list()
  for (i in seq_len(nrow(marker_df))) {
    gid <- marker_df$gene_id[i]
    tvec <- split_annotation_terms(marker_df[[term_col]][i], pattern = split_pattern)
    if (length(tvec) == 0) next
    for (tt in tvec) {
      out[[tt]] <- c(out[[tt]], gid)
    }
  }
  out <- lapply(out, unique)
  out
}

annotate_go_term_names <- function(enrich_df) {
  if (nrow(enrich_df) == 0) return(enrich_df)
  enrich_df$go_term_name <- NA_character_
  go_idx <- which(enrich_df$source == "GO" & !is.na(enrich_df$term) & nzchar(enrich_df$term))
  if (length(go_idx) == 0) return(enrich_df)
  go_terms <- unique(enrich_df$term[go_idx])

  # 1) Direct lookup in local GO.db (current IDs present in installed release).
  if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    valid <- intersect(go_terms, AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"))
    if (length(valid) > 0) {
      gtab <- AnnotationDbi::select(
        GO.db::GO.db,
        keys = valid,
        columns = c("TERM"),
        keytype = "GOID"
      )
      gtab <- unique(gtab[, c("GOID", "TERM"), drop = FALSE])
      go_map <- setNames(gtab$TERM, gtab$GOID)
      enrich_df$go_term_name[go_idx] <- unname(go_map[enrich_df$term[go_idx]])
    }
  }

  # 2) Alternate/merged GO IDs frequently seen in older eggNOG annotations.
  alt_to_current <- c(
    "GO:0015037" = "GO:0015035", # alt ID in AmiGO
    "GO:0000780" = "GO:0000779"  # alt ID in AmiGO
  )
  missing_idx <- go_idx[is.na(enrich_df$go_term_name[go_idx]) | !nzchar(enrich_df$go_term_name[go_idx])]
  if (length(missing_idx) > 0 && requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    alt_terms <- enrich_df$term[missing_idx]
    alt_curr <- unname(alt_to_current[alt_terms])
    ok <- !is.na(alt_curr)
    if (any(ok)) {
      valid2 <- intersect(unique(alt_curr[ok]), AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"))
      if (length(valid2) > 0) {
        gtab2 <- AnnotationDbi::select(
          GO.db::GO.db,
          keys = valid2,
          columns = c("TERM"),
          keytype = "GOID"
        )
        gtab2 <- unique(gtab2[, c("GOID", "TERM"), drop = FALSE])
        go_map2 <- setNames(gtab2$TERM, gtab2$GOID)
        fill_vals <- go_map2[alt_curr]
        fill_idx <- missing_idx[ok & !is.na(fill_vals) & nzchar(fill_vals)]
        fill_names <- as.character(fill_vals[ok & !is.na(fill_vals) & nzchar(fill_vals)])
        if (length(fill_idx) > 0) {
          enrich_df$go_term_name[fill_idx] <- fill_names
        }
      }
    }
  }

  # 3) Curated fallback labels for legacy obsolete IDs missing in current GO.db.
  fallback_names <- c(
    "GO:0044424" = "obsolete intracellular part",
    "GO:0044464" = "obsolete cell part",
    "GO:0044444" = "obsolete cytoplasmic part",
    "GO:0044422" = "obsolete organelle part",
    "GO:0043987" = "obsolete histone H3-S10 phosphorylation",
    "GO:0016572" = "obsolete histone phosphorylation",
    "GO:0035404" = "histone-serine phosphorylation"
  )
  missing_idx2 <- go_idx[is.na(enrich_df$go_term_name[go_idx]) | !nzchar(enrich_df$go_term_name[go_idx])]
  if (length(missing_idx2) > 0) {
    fill2 <- unname(fallback_names[enrich_df$term[missing_idx2]])
    ok2 <- !is.na(fill2) & nzchar(fill2)
    if (any(ok2)) enrich_df$go_term_name[missing_idx2[ok2]] <- fill2[ok2]
  }

  # Final fallback: keep ID as label so plots/tables are never blank.
  still_missing <- go_idx[is.na(enrich_df$go_term_name[go_idx]) | !nzchar(enrich_df$go_term_name[go_idx])]
  if (length(still_missing) > 0) {
    enrich_df$go_term_name[still_missing] <- enrich_df$term[still_missing]
  }
  enrich_df
}

drop_obsolete_go_terms <- function(enrich_df) {
  if (nrow(enrich_df) == 0) return(enrich_df)
  if (!all(c("source", "go_term_name") %in% colnames(enrich_df))) return(enrich_df)
  is_obs <- enrich_df$source == "GO" &
    !is.na(enrich_df$go_term_name) &
    grepl("^obsolete\\b", tolower(trimws(enrich_df$go_term_name)))
  enrich_df[!is_obs, , drop = FALSE]
}

run_marker_functional_enrichment <- function(marker_df,
                                             term_specs = list(
                                               list(source = "GO", col = "GOs", split = ","),
                                               list(source = "KEGG_ko", col = "KEGG_ko", split = ","),
                                               list(source = "PFAM", col = "PFAMs", split = ","),
                                               list(source = "COG", col = "COG_category", split = ",")
                                             ),
                                             min_term_genes = 5L,
                                             marker_classes = NULL,
                                             universe_genes = NULL) {
  schema <- data.frame(
    marker_group = character(), marker_class = character(), source = character(),
    term = character(), set_size = integer(), term_genes = integer(), overlap = integer(),
    overlap_gene_ids = character(),
    odds_ratio = numeric(), p_value = numeric(), fdr = numeric(), test_method = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(marker_df) == 0) return(schema)

  marker_df <- marker_df[!is.na(marker_df$gene_id) & nzchar(marker_df$gene_id), , drop = FALSE]
  if (!"marker_group" %in% colnames(marker_df)) return(schema)

  universe <- if (is.null(universe_genes)) unique(marker_df$gene_id) else unique(universe_genes)
  universe <- universe[!is.na(universe) & nzchar(universe)]
  marker_df <- marker_df[marker_df$gene_id %in% universe, , drop = FALSE]
  rows <- list()
  idx <- 0L

  if (is.null(marker_classes)) {
    sub_marker <- marker_df[marker_df$marker_flag == TRUE, , drop = FALSE]
    if (nrow(sub_marker) == 0) return(schema)
    class_label <- "marker_flag_true"
    sets <- split(sub_marker$gene_id, sub_marker$marker_group)
  } else {
    sets <- list()
    class_label <- NA_character_
  }

  class_iter <- if (is.null(marker_classes)) class_label else marker_classes
  for (mc in class_iter) {
    if (is.null(marker_classes)) {
      sub_class <- sub_marker
    } else {
      if (!"marker_class" %in% colnames(marker_df)) next
      sub_class <- marker_df[marker_df$marker_class == mc, , drop = FALSE]
      if (nrow(sub_class) == 0) next
      sets <- split(sub_class$gene_id, sub_class$marker_group)
    }

    for (spec in term_specs) {
      term_map <- build_term_gene_map(marker_df, term_col = spec$col, split_pattern = spec$split)
      if (length(term_map) == 0) next

      keep_terms <- names(term_map)[vapply(term_map, function(g) length(intersect(g, universe)) >= min_term_genes, logical(1))]
      if (length(keep_terms) == 0) next

      for (mg in names(sets)) {
        set_genes <- intersect(unique(sets[[mg]]), universe)
        if (length(set_genes) == 0) next

        for (tt in keep_terms) {
          term_genes <- intersect(term_map[[tt]], universe)
          a <- length(intersect(set_genes, term_genes))
          if (a == 0) next
          overlap_genes <- sort(intersect(set_genes, term_genes))
          b <- length(setdiff(set_genes, term_genes))
          c <- length(setdiff(term_genes, set_genes))
          d <- length(setdiff(universe, union(set_genes, term_genes)))
          tab <- matrix(c(a, b, c, d), nrow = 2)
          ft <- stats::fisher.test(tab, alternative = "greater")

          idx <- idx + 1L
          rows[[idx]] <- data.frame(
            marker_group = mg,
            marker_class = mc,
            source = spec$source,
            term = tt,
            set_size = length(set_genes),
            term_genes = length(term_genes),
            overlap = a,
            overlap_gene_ids = paste(overlap_genes, collapse = ","),
            odds_ratio = unname(ft$estimate),
            p_value = ft$p.value,
            test_method = "fisher_greater",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(rows) == 0) return(schema)
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out <- annotate_go_term_names(out)
  out <- drop_obsolete_go_terms(out)
  out <- out[, c("marker_group", "marker_class", "source", "term", "go_term_name", "set_size", "term_genes", "overlap", "odds_ratio", "p_value", "fdr", "test_method", "overlap_gene_ids")]
  out
}

run_module_functional_enrichment <- function(module_df,
                                             annotation_df,
                                             term_specs = list(
                                               list(source = "GO", col = "GOs", split = ","),
                                               list(source = "KEGG_ko", col = "KEGG_ko", split = ","),
                                               list(source = "PFAM", col = "PFAMs", split = ","),
                                               list(source = "COG", col = "COG_category", split = ",")
                                             ),
                                             min_term_genes = 5L,
                                             min_set_size = 10L) {
  schema <- data.frame(
    module_id = character(), source = character(), term = character(), go_term_name = character(),
    set_size = integer(), term_genes = integer(), overlap = integer(),
    odds_ratio = numeric(), p_value = numeric(), fdr = numeric(), test_method = character(),
    overlap_gene_ids = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(module_df) == 0 || nrow(annotation_df) == 0) return(schema)
  if (!all(c("gene_id", "module_id") %in% colnames(module_df))) return(schema)
  if (!"gene_id" %in% colnames(annotation_df)) return(schema)

  annotation_df <- annotation_df[!is.na(annotation_df$gene_id) & nzchar(annotation_df$gene_id), , drop = FALSE]
  annotation_df <- annotation_df[!duplicated(annotation_df$gene_id), , drop = FALSE]
  universe <- intersect(unique(module_df$gene_id), unique(annotation_df$gene_id))
  if (length(universe) == 0) return(schema)

  modules <- split(module_df$gene_id, module_df$module_id)
  rows <- list()
  idx <- 0L

  for (spec in term_specs) {
    term_map <- build_term_gene_map(annotation_df, term_col = spec$col, split_pattern = spec$split)
    if (length(term_map) == 0) next
    keep_terms <- names(term_map)[vapply(term_map, function(g) length(intersect(g, universe)) >= min_term_genes, logical(1))]
    if (length(keep_terms) == 0) next

    for (mid in names(modules)) {
      set_genes <- intersect(unique(modules[[mid]]), universe)
      if (length(set_genes) < min_set_size) next

      for (tt in keep_terms) {
        term_genes <- intersect(term_map[[tt]], universe)
        overlap_genes <- sort(intersect(set_genes, term_genes))
        a <- length(overlap_genes)
        if (a == 0) next
        b <- length(setdiff(set_genes, term_genes))
        c <- length(setdiff(term_genes, set_genes))
        d <- length(setdiff(universe, union(set_genes, term_genes)))
        tab <- matrix(c(a, b, c, d), nrow = 2)
        ft <- stats::fisher.test(tab, alternative = "greater")

        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          module_id = mid,
          source = spec$source,
          term = tt,
          set_size = length(set_genes),
          term_genes = length(term_genes),
          overlap = a,
          odds_ratio = unname(ft$estimate),
          p_value = ft$p.value,
          test_method = "fisher_greater",
          overlap_gene_ids = paste(overlap_genes, collapse = ","),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) return(schema)
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out <- annotate_go_term_names(out)
  out <- drop_obsolete_go_terms(out)
  out <- out[, c("module_id", "source", "term", "go_term_name", "set_size", "term_genes", "overlap", "odds_ratio", "p_value", "fdr", "test_method", "overlap_gene_ids")]
  out
}

run_marker_enrichment_permutation <- function(cluster_members_df, marker_df, coords_df, n_perm = 1000, seed = 42) {
  set.seed(seed)
  if (nrow(cluster_members_df) == 0 || nrow(marker_df) == 0 || nrow(coords_df) == 0) return(data.frame())

  marker_df <- marker_df[marker_df$marker_flag == TRUE, , drop = FALSE]
  marker_sets <- split(marker_df$gene_id, marker_df$marker_group)
  coords_pool <- split(coords_df$gene_id, coords_df$chr)

  rows <- list()
  idx <- 0
  for (cid in unique(cluster_members_df$set_id)) {
    genes <- unique(cluster_members_df$gene_id[cluster_members_df$set_id == cid])
    sub <- coords_df[coords_df$gene_id %in% genes, , drop = FALSE]
    if (nrow(sub) < 2) next
    chr_counts <- table(sub$chr)

    for (mg in names(marker_sets)) {
      obs <- length(intersect(genes, marker_sets[[mg]]))
      overlap_genes <- sort(intersect(genes, marker_sets[[mg]]))
      perm_hits <- numeric(n_perm)
      for (i in seq_len(n_perm)) {
        pick <- character(0)
        valid <- TRUE
        for (ch in names(chr_counts)) {
          n <- as.integer(chr_counts[[ch]])
          pool <- coords_pool[[ch]]
          if (length(pool) < n) {
            valid <- FALSE
            break
          }
          pick <- c(pick, sample(pool, n, replace = FALSE))
        }
        perm_hits[i] <- if (valid) length(intersect(pick, marker_sets[[mg]])) else NA_real_
      }
      perm_hits <- perm_hits[!is.na(perm_hits)]
      if (length(perm_hits) == 0) next
      p <- (sum(perm_hits >= obs) + 1) / (length(perm_hits) + 1)

      idx <- idx + 1
      rows[[idx]] <- data.frame(
        cluster_id = cid,
        marker_group = mg,
        observed_overlap = obs,
        overlap_gene_ids = paste(overlap_genes, collapse = ","),
        empirical_p = p,
        test_method = "perm_chr_stratified_overlap",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$empirical_p, method = "BH")
  out <- out[, c("cluster_id", "marker_group", "observed_overlap", "empirical_p", "test_method", "fdr", "overlap_gene_ids")]
  out
}
