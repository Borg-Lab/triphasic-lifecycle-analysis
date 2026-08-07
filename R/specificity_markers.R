calc_tau <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || max(x, na.rm = TRUE) <= 0) return(NA_real_)
  x <- x / max(x, na.rm = TRUE)
  sum(1 - x, na.rm = TRUE) / (length(x) - 1)
}

cfg_num <- function(cfg, key, default) {
  if (is.null(cfg[[key]]) || !is.finite(as.numeric(cfg[[key]]))) return(as.numeric(default))
  as.numeric(cfg[[key]])
}

cfg_chr <- function(cfg, key, default) {
  if (is.null(cfg[[key]])) return(as.character(default))
  x <- as.character(cfg[[key]])
  if (!nzchar(x)) return(as.character(default))
  x
}

normalize_marker_de_mode <- function(mode) {
  m <- tolower(trimws(as.character(mode)))
  if (!m %in% c("any", "relevant")) return("any")
  m
}

build_relevant_contrast_lookup <- function(contrast_df) {
  if (is.null(contrast_df) || nrow(contrast_df) == 0) return(list())
  req <- c("contrast_id", "numerator", "denominator")
  if (!all(req %in% colnames(contrast_df))) return(list())

  rows <- contrast_df[, req, drop = FALSE]
  rows$numerator <- as.character(rows$numerator)
  rows$denominator <- as.character(rows$denominator)
  rows$contrast_id <- as.character(rows$contrast_id)
  rows <- rows[!is.na(rows$contrast_id) & nzchar(rows$contrast_id), , drop = FALSE]
  if (nrow(rows) == 0) return(list())

  groups <- sort(unique(c(rows$numerator, rows$denominator)))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  out <- setNames(vector("list", length(groups)), groups)
  for (g in groups) {
    sub_num <- rows[rows$numerator == g, c("contrast_id"), drop = FALSE]
    if (nrow(sub_num) > 0) sub_num$sign_to_group <- 1
    sub_den <- rows[rows$denominator == g, c("contrast_id"), drop = FALSE]
    if (nrow(sub_den) > 0) sub_den$sign_to_group <- -1
    sub <- rbind(sub_num, sub_den)
    if (nrow(sub) == 0) {
      out[[g]] <- data.frame(contrast_id = character(), sign_to_group = numeric(), stringsAsFactors = FALSE)
    } else {
      sub <- sub[!duplicated(sub$contrast_id), , drop = FALSE]
      out[[g]] <- sub
    }
  }
  out
}

build_marker_table <- function(group_tpm_df, de_df, analysis_cfg, func_anno_df = NULL, contrast_df = NULL) {
  schema <- data.frame(
    gene_id = character(), marker_group = character(), tau = numeric(),
    best_log2fc = numeric(), best_fdr = numeric(),
    de_available = logical(), de_testable = logical(),
    marker_specificity_flag = logical(), marker_flag = logical(),
    marker_class = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(group_tpm_df) == 0) return(schema)

  gcols <- setdiff(colnames(group_tpm_df), "gene_id")
  mat <- as.matrix(group_tpm_df[, gcols, drop = FALSE])
  storage.mode(mat) <- "numeric"

  tau <- apply(mat, 1, calc_tau)
  max_idx <- apply(mat, 1, function(z) {
    if (all(is.na(z))) return(NA_integer_)
    which.max(z)
  })
  marker_group <- ifelse(is.na(max_idx), NA_character_, gsub("^mean_tpm_", "", gcols[max_idx]))

  best <- data.frame(
    gene_id = group_tpm_df$gene_id,
    best_log2fc = NA_real_,
    best_fdr = NA_real_,
    de_available = FALSE,
    de_testable = FALSE,
    stringsAsFactors = FALSE
  )

  marker_de_mode <- normalize_marker_de_mode(cfg_chr(analysis_cfg, "marker_de_mode", "any"))
  rel_lookup <- if (marker_de_mode == "relevant") build_relevant_contrast_lookup(contrast_df) else list()

  if (nrow(de_df) > 0) {
    de_split_all <- split(de_df, de_df$gene_id)
    de_valid <- de_df[is.finite(de_df$log2fc) & is.finite(de_df$fdr), , drop = FALSE]
    de_split_valid <- split(de_valid, de_valid$gene_id)

    for (i in seq_len(nrow(best))) {
      gi <- best$gene_id[i]
      mg <- marker_group[i]

      if (!gi %in% names(de_split_all)) next
      d_all <- de_split_all[[gi]]
      d_val <- if (gi %in% names(de_split_valid)) de_split_valid[[gi]] else d_all[0, , drop = FALSE]

      if (marker_de_mode == "relevant" && !is.null(rel_lookup[[mg]])) {
        rel <- rel_lookup[[mg]]
        d_all_r <- d_all[d_all$contrast_id %in% rel$contrast_id, , drop = FALSE]
        d_val_r <- d_val[d_val$contrast_id %in% rel$contrast_id, , drop = FALSE]

        best$de_available[i] <- nrow(d_all_r) > 0
        best$de_testable[i] <- nrow(d_val_r) > 0
        if (nrow(d_val_r) == 0) next

        m <- match(d_val_r$contrast_id, rel$contrast_id)
        sign_to_group <- rel$sign_to_group[m]
        effect_to_group <- d_val_r$log2fc * sign_to_group
        keep <- is.finite(effect_to_group) & is.finite(d_val_r$fdr)
        if (!any(keep)) next
        d_val_r <- d_val_r[keep, , drop = FALSE]
        effect_to_group <- effect_to_group[keep]
        j <- which.max(effect_to_group)
        best$best_log2fc[i] <- effect_to_group[j]
        best$best_fdr[i] <- d_val_r$fdr[j]
      } else {
        best$de_available[i] <- TRUE
        if (nrow(d_val) == 0) next
        best$de_testable[i] <- TRUE
        j <- which.max(abs(d_val$log2fc))
        best$best_log2fc[i] <- d_val$log2fc[j]
        best$best_fdr[i] <- d_val$fdr[j]
      }
    }
  }

  out <- data.frame(
    gene_id = group_tpm_df$gene_id,
    marker_group = marker_group,
    tau = tau,
    best_log2fc = best$best_log2fc,
    best_fdr = best$best_fdr,
    de_available = best$de_available,
    de_testable = best$de_testable,
    stringsAsFactors = FALSE
  )

  out$marker_specificity_flag <- with(out, !is.na(tau) & tau >= as.numeric(analysis_cfg$marker_tau_threshold))
  out$marker_flag <- with(
    out,
    marker_specificity_flag &
      de_testable &
      !is.na(best_log2fc) & abs(best_log2fc) >= as.numeric(analysis_cfg$marker_min_log2fc) &
      !is.na(best_fdr) & best_fdr <= as.numeric(analysis_cfg$marker_max_fdr)
  )

  out$marker_class <- "not_marker"
  out$marker_class[out$marker_specificity_flag & !out$de_available] <- "specificity_only_no_de"
  out$marker_class[out$marker_specificity_flag & out$de_available & !out$de_testable] <- "specificity_only_de_na"
  out$marker_class[out$marker_specificity_flag & out$de_testable & !out$marker_flag] <- "specificity_only_de_not_sig"
  out$marker_class[out$marker_flag] <- "de_supported_marker"

  if (!is.null(func_anno_df) && nrow(func_anno_df) > 0 && "gene_id" %in% colnames(func_anno_df)) {
    keep_cols <- setdiff(colnames(func_anno_df), c("query"))
    keep_cols <- c("gene_id", setdiff(keep_cols, "gene_id"))
    func_anno_df <- func_anno_df[, keep_cols, drop = FALSE]
    m <- match(out$gene_id, func_anno_df$gene_id)
    extra <- func_anno_df[m, setdiff(colnames(func_anno_df), "gene_id"), drop = FALSE]
    out <- cbind(out, extra)
  }

  out
}

build_marker_table_strict <- function(group_tpm_df, de_df, analysis_cfg, func_anno_df = NULL, contrast_df = NULL) {
  out <- build_marker_table(group_tpm_df, de_df, analysis_cfg, func_anno_df = func_anno_df, contrast_df = contrast_df)
  if (nrow(out) == 0) return(out)

  gcols <- grep("^mean_tpm_", colnames(group_tpm_df), value = TRUE)
  mat <- as.matrix(group_tpm_df[, gcols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- group_tpm_df$gene_id

  tpm_floor <- cfg_num(analysis_cfg, "marker_strict_tpm_floor", 1)
  min_max_tpm <- cfg_num(analysis_cfg, "marker_strict_min_max_tpm", 5)
  max_detected <- cfg_num(analysis_cfg, "marker_strict_max_detected_stages", 4)

  mat_floor <- mat
  mat_floor[!is.finite(mat_floor) | mat_floor < tpm_floor] <- 0

  tau_floor <- apply(mat_floor, 1, calc_tau)
  detected_n <- rowSums(mat_floor >= tpm_floor, na.rm = TRUE)
  max_tpm <- apply(mat_floor, 1, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE))
  max_idx <- apply(mat_floor, 1, function(z) {
    if (all(is.na(z)) || max(z, na.rm = TRUE) <= 0) return(NA_integer_)
    which.max(z)
  })
  marker_group_floor <- ifelse(is.na(max_idx), NA_character_, gsub("^mean_tpm_", "", gcols[max_idx]))

  m <- match(out$gene_id, rownames(mat_floor))
  out$strict_tau <- tau_floor[m]
  out$strict_detected_stages <- as.integer(detected_n[m])
  out$strict_max_tpm <- as.numeric(max_tpm[m])
  out$strict_marker_group <- marker_group_floor[m]
  out$strict_tpm_floor <- tpm_floor
  out$strict_max_detected_stages_cutoff <- max_detected
  out$strict_min_max_tpm_cutoff <- min_max_tpm

  out$strict_specificity_flag <- with(
    out,
    !is.na(strict_tau) &
      strict_tau >= as.numeric(analysis_cfg$marker_tau_threshold) &
      !is.na(strict_max_tpm) & strict_max_tpm >= min_max_tpm &
      !is.na(strict_detected_stages) & strict_detected_stages <= max_detected
  )

  out$strict_marker_flag <- with(
    out,
    strict_specificity_flag &
      de_testable &
      !is.na(best_log2fc) & abs(best_log2fc) >= as.numeric(analysis_cfg$marker_min_log2fc) &
      !is.na(best_fdr) & best_fdr <= as.numeric(analysis_cfg$marker_max_fdr)
  )

  out$strict_marker_class <- "strict_not_marker"
  out$strict_marker_class[out$strict_specificity_flag & !out$de_available] <- "strict_specificity_only_no_de"
  out$strict_marker_class[out$strict_specificity_flag & out$de_available & !out$de_testable] <- "strict_specificity_only_de_na"
  out$strict_marker_class[out$strict_specificity_flag & out$de_testable & !out$strict_marker_flag] <- "strict_specificity_only_de_not_sig"
  out$strict_marker_class[out$strict_marker_flag] <- "strict_de_supported_marker"

  out
}

apply_strict_marker_mode <- function(marker_df) {
  if (nrow(marker_df) == 0) return(marker_df)
  needed <- c("strict_marker_group", "strict_tau", "strict_specificity_flag", "strict_marker_flag", "strict_marker_class")
  if (!all(needed %in% colnames(marker_df))) return(marker_df)

  marker_df$orig_marker_group <- marker_df$marker_group
  marker_df$orig_tau <- marker_df$tau
  marker_df$orig_marker_specificity_flag <- marker_df$marker_specificity_flag
  marker_df$orig_marker_flag <- marker_df$marker_flag
  marker_df$orig_marker_class <- marker_df$marker_class

  marker_df$marker_group <- marker_df$strict_marker_group
  marker_df$tau <- marker_df$strict_tau
  marker_df$marker_specificity_flag <- marker_df$strict_specificity_flag
  marker_df$marker_flag <- marker_df$strict_marker_flag
  marker_df$marker_class <- marker_df$strict_marker_class
  marker_df
}

build_marker_table_mode <- function(group_tpm_df, de_df, analysis_cfg, func_anno_df = NULL, contrast_df = NULL) {
  mode <- tolower(as.character(if (!is.null(analysis_cfg$marker_mode)) analysis_cfg$marker_mode else "default"))
  if (mode == "strict") {
    strict_df <- build_marker_table_strict(group_tpm_df, de_df, analysis_cfg, func_anno_df = func_anno_df, contrast_df = contrast_df)
    return(apply_strict_marker_mode(strict_df))
  }
  build_marker_table(group_tpm_df, de_df, analysis_cfg, func_anno_df = func_anno_df, contrast_df = contrast_df)
}

marker_sets_from_table <- function(marker_df) {
  marker_df <- marker_df[isTRUE(marker_df$marker_flag) | marker_df$marker_flag == TRUE, , drop = FALSE]
  split(marker_df$gene_id, marker_df$marker_group)
}
