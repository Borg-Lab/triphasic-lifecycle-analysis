#!/usr/bin/env Rscript
source("R/io.R")

paths <- read_paths_config(Sys.getenv("PATHS_CONFIG", unset = "config/paths.yml"))
analysis <- read_analysis_config(Sys.getenv("ANALYSIS_CONFIG", unset = "config/analysis.yml"))
design <- read_design_table(Sys.getenv("DESIGN_CONFIG", unset = "config/design.tsv"))

out_dir <- file.path(paths$output_root, "module_conservation")
ensure_dir(out_dir)

boot_n <- as.integer(Sys.getenv("MODULE_CONSERVATION_BOOT_N", unset = "100"))
perm_n <- as.integer(Sys.getenv("MODULE_CONSERVATION_PERM_N", unset = "200"))
max_resample_genes <- as.integer(Sys.getenv("MODULE_CONSERVATION_MAX_RESAMPLE_GENES", unset = "500"))
seed <- as.integer(Sys.getenv("MODULE_CONSERVATION_SEED", unset = "20260720"))
set.seed(seed)

power <- if (!is.null(analysis$wgcna_power)) as.integer(analysis$wgcna_power) else 12L
`%||%` <- function(x, y) if (is.null(x)) y else x

membership_path <- file.path(paths$output_root, "modules", "module_membership_labeled.tsv")
if (!file.exists(membership_path)) {
  membership_path <- file.path(paths$output_root, "modules", "module_membership.tsv")
}
modules <- safe_read_table(membership_path)
if (nrow(modules) == 0) stop("Missing module membership table.", call. = FALSE)
if (!"gene_id" %in% colnames(modules) || !"module_id" %in% colnames(modules)) {
  stop("Module membership table must contain gene_id and module_id.", call. = FALSE)
}
modules$gene_id <- normalize_gene_id(modules$gene_id)
if (!"module_label" %in% colnames(modules)) modules$module_label <- modules$module_id
if (!"dominant_stage" %in% colnames(modules)) modules$dominant_stage <- NA_character_
modules <- modules[!is.na(modules$gene_id) & nzchar(modules$gene_id), , drop = FALSE]
modules <- modules[!duplicated(modules$gene_id), , drop = FALSE]

expr_inputs <- prepare_expression_inputs(paths, design, analysis)
counts <- expr_inputs$counts
if (nrow(counts) == 0) stop("No count matrix available after filtering.", call. = FALSE)
mat <- as.matrix(counts[, setdiff(colnames(counts), "gene_id"), drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- normalize_gene_id(counts$gene_id)
mat <- log2(mat + 1)

common_genes <- intersect(rownames(mat), modules$gene_id)
mat <- mat[common_genes, , drop = FALSE]
modules <- modules[match(common_genes, modules$gene_id), , drop = FALSE]

stage_groups <- list(
  CSP = c("Bm3235c_carposporophyte"),
  TETRA = c("Bm3235t_vegtissue", "Bm3235t_stichidia"),
  SPORES = c("Bm3235t_carpospores", "Bm3235t_tetraspores"),
  PERICARP = c("Bm3235f_pericarp"),
  CYSTO_MATERNAL = c("Bm3235f_cystocarps_imm", "Bm3235f_cystocarps_mat", "Bm3235f_pericarp"),
  GAM = c(
    "Bm3235f_vegtissue", "Bm3235f_procarps", "Bm3235m_vegtissue",
    "Bm3235m_spermatangia", "Bm3235m_spermatia"
  )
)

stage_order <- c(
  "Bm3235f_vegtissue", "Bm3235f_procarps", "Bm3235m_vegtissue",
  "Bm3235m_spermatangia", "Bm3235m_spermatia", "Bm3235f_pericarp",
  "Bm3235f_cystocarps_imm", "Bm3235f_cystocarps_mat", "Bm3235c_carposporophyte",
  "Bm3235t_carpospores", "Bm3235t_vegtissue", "Bm3235t_stichidia", "Bm3235t_tetraspores"
)
stage_alone_groups <- stats::setNames(as.list(stage_order), paste0("STAGE__", stage_order))
stage_groups <- c(stage_groups, stage_alone_groups)

samples_for_group <- function(stages) {
  x <- design$sample_id[design$stage %in% stages & as.logical(design$include)]
  intersect(x, colnames(mat))
}
group_samples <- lapply(stage_groups, samples_for_group)
group_counts <- data.frame(
  group = names(group_samples),
  stages = vapply(stage_groups, paste, collapse = ";", FUN.VALUE = character(1)),
  n_samples = vapply(group_samples, length, FUN.VALUE = integer(1)),
  samples = vapply(group_samples, paste, collapse = ";", FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)
write_tsv(group_counts, file.path(out_dir, "module_conservation_group_samples.tsv"))
base_group_names <- c("CSP", "TETRA", "SPORES", "PERICARP", "CYSTO_MATERNAL", "GAM")
stage_group_names <- paste0("STAGE__", stage_order)

module_sort_key <- if ("module_rank" %in% colnames(modules)) modules$module_rank else modules$module_label
module_order <- unique(modules$module_label[order(module_sort_key, modules$module_label)])
module_genes <- split(modules$gene_id, modules$module_label)
module_genes <- module_genes[module_order[module_order %in% names(module_genes)]]
module_meta <- modules[!duplicated(modules$module_label), , drop = FALSE]
rownames(module_meta) <- module_meta$module_label

is_csp_module <- function(label) {
  x <- module_meta[label, , drop = FALSE]
  any(grepl("carposporophyte", c(x$module_label, x$dominant_stage), ignore.case = TRUE), na.rm = TRUE)
}
target_modules <- names(module_genes)[vapply(names(module_genes), is_csp_module, FUN.VALUE = logical(1))]

cor_to_connectivity <- function(expr_gene_sample, sample_ids, power, scale = TRUE) {
  sample_ids <- intersect(sample_ids, colnames(expr_gene_sample))
  genes <- rownames(expr_gene_sample)
  if (length(sample_ids) < 3 || length(genes) < 3) {
    return(list(
      cor_summary = data.frame(
        n_genes = length(genes), n_samples = length(sample_ids),
        mean_signed_cor = NA_real_, median_signed_cor = NA_real_, iqr_signed_cor = NA_real_,
        median_abs_cor = NA_real_, stringsAsFactors = FALSE
      ),
      k = setNames(rep(NA_real_, length(genes)), genes)
    ))
  }

  cm <- suppressWarnings(stats::cor(t(expr_gene_sample[, sample_ids, drop = FALSE]), use = "pairwise.complete.obs"))
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  upper <- cm[upper.tri(cm)]
  adj <- ((1 + cm) / 2) ^ power
  diag(adj) <- 0
  k <- rowSums(adj, na.rm = TRUE)
  if (isTRUE(scale) && is.finite(max(k, na.rm = TRUE)) && max(k, na.rm = TRUE) > 0) {
    k <- k / max(k, na.rm = TRUE)
  }

  list(
    cor_summary = data.frame(
      n_genes = length(genes),
      n_samples = length(sample_ids),
      mean_signed_cor = mean(upper, na.rm = TRUE),
      median_signed_cor = stats::median(upper, na.rm = TRUE),
      iqr_signed_cor = stats::IQR(upper, na.rm = TRUE),
      median_abs_cor = stats::median(abs(upper), na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    k = k
  )
}

rank_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(c(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(stats::cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

observed_stats <- list()
gene_conn <- list()
conn_by_module <- list()

for (module_label in names(module_genes)) {
  genes <- intersect(module_genes[[module_label]], rownames(mat))
  expr_mod <- mat[genes, , drop = FALSE]
  conn_by_module[[module_label]] <- list()

  observed_group_names <- c(base_group_names, stage_group_names)

  for (group_name in observed_group_names) {
    res <- cor_to_connectivity(expr_mod, group_samples[[group_name]], power)
    conn_by_module[[module_label]][[group_name]] <- res$k

    observed_stats[[length(observed_stats) + 1]] <- data.frame(
      module_label = module_label,
      module_id = module_meta[module_label, "module_id"],
      module_rank = if ("module_rank" %in% colnames(module_meta)) module_meta[module_label, "module_rank"] else NA,
      dominant_stage = if ("dominant_stage" %in% colnames(module_meta)) module_meta[module_label, "dominant_stage"] else NA,
      csp_associated = module_label %in% target_modules,
      group = group_name,
      res$cor_summary,
      stringsAsFactors = FALSE
    )

    gene_conn[[length(gene_conn) + 1]] <- data.frame(
      module_label = module_label,
      group = group_name,
      gene_id = names(res$k),
      k_within_scaled = as.numeric(res$k),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, observed_stats)
gene_conn_df <- do.call(rbind, gene_conn)

sim_rows <- list()
for (module_label in names(conn_by_module)) {
  csp_k <- conn_by_module[[module_label]][["CSP"]]
  observed_group_names <- c(base_group_names, stage_group_names)

  for (group_name in observed_group_names) {
    gk <- conn_by_module[[module_label]][[group_name]]
    rc <- if (group_name == "CSP") c(rho = 1, p = NA_real_) else rank_cor(csp_k[names(gk)], gk)
    sim_rows[[length(sim_rows) + 1]] <- data.frame(
      module_label = module_label,
      group = group_name,
      connectivity_rank_cor_vs_CSP = as.numeric(rc["rho"]),
      connectivity_rank_cor_p_vs_CSP = as.numeric(rc["p"]),
      stringsAsFactors = FALSE
    )
  }
}
sim_df <- do.call(rbind, sim_rows)
summary_df <- merge(summary_df, sim_df, by = c("module_label", "group"), all.x = TRUE, sort = FALSE)

get_sim <- function(module_label, group_name) {
  x <- sim_df$connectivity_rank_cor_vs_CSP[sim_df$module_label == module_label & sim_df$group == group_name]
  if (length(x) == 0) NA_real_ else x[1]
}
contrast_rows <- list()
contrast_controls <- c("PERICARP", "CYSTO_MATERNAL", "GAM")
tetra_groups <- c("TETRA", "SPORES")
for (module_label in names(module_genes)) {
  for (tetra_group in tetra_groups) {
    for (control_group in contrast_controls) {
      contrast_rows[[length(contrast_rows) + 1]] <- data.frame(
        module_label = module_label,
        csp_associated = module_label %in% target_modules,
        tetra_group = tetra_group,
        control_group = control_group,
        similarity_tetra = get_sim(module_label, tetra_group),
        similarity_control = get_sim(module_label, control_group),
        delta_connectivity_similarity = get_sim(module_label, tetra_group) - get_sim(module_label, control_group),
        stringsAsFactors = FALSE
      )
    }
  }
}
contrast_df <- do.call(rbind, contrast_rows)

stage_group_names <- paste0("STAGE__", stage_order)
stage_similarity_rows <- list()
for (module_label in names(module_genes)) {
  for (stage_name in stage_order) {
    group_name <- paste0("STAGE__", stage_name)
    stage_similarity_rows[[length(stage_similarity_rows) + 1]] <- data.frame(
      module_label = module_label,
      csp_associated = module_label %in% target_modules,
      stage = stage_name,
      group = group_name,
      n_samples = length(group_samples[[group_name]]),
      similarity_stage = get_sim(module_label, group_name),
      similarity_pericarp = get_sim(module_label, "PERICARP"),
      delta_stage_vs_pericarp = get_sim(module_label, group_name) - get_sim(module_label, "PERICARP"),
      stringsAsFactors = FALSE
    )
  }
}
stage_similarity_df <- do.call(rbind, stage_similarity_rows)

write_tsv(summary_df, file.path(out_dir, "module_conservation_summary.tsv"))
write_tsv(contrast_df, file.path(out_dir, "module_conservation_contrasts.tsv"))
write_tsv(stage_similarity_df, file.path(out_dir, "module_conservation_stage_similarity.tsv"))
write_tsv(gene_conn_df, file.path(out_dir, "module_conservation_gene_connectivity.tsv"))

sample_genes <- function(genes, max_genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) > max_genes) sample(genes, max_genes) else genes
}

bootstrap_rows <- list()
if (boot_n > 0 && length(target_modules) > 0) {
  boot_groups <- c("CSP", "TETRA", "SPORES", "PERICARP", "CYSTO_MATERNAL", "GAM")
  n_balanced <- min(vapply(group_samples[boot_groups], length, FUN.VALUE = integer(1)))
  for (module_label in target_modules) {
    genes <- sample_genes(module_genes[[module_label]], max_resample_genes)
    genes_used <- length(genes)
    expr_mod <- mat[genes, , drop = FALSE]
    for (iter in seq_len(boot_n)) {
      boot_samples <- lapply(group_samples[boot_groups], function(x) sample(x, n_balanced, replace = length(x) < n_balanced))
      k_list <- lapply(boot_samples, function(x) cor_to_connectivity(expr_mod, x, power)$k)
      for (tetra_group in tetra_groups) {
        for (control_group in contrast_controls) {
          sim_tetra <- rank_cor(k_list[["CSP"]], k_list[[tetra_group]])["rho"]
          sim_control <- rank_cor(k_list[["CSP"]], k_list[[control_group]])["rho"]
          bootstrap_rows[[length(bootstrap_rows) + 1]] <- data.frame(
            module_label = module_label,
            iteration = iter,
            n_balanced_samples = n_balanced,
            module_gene_count = length(module_genes[[module_label]]),
            genes_used = genes_used,
            tetra_group = tetra_group,
            control_group = control_group,
            similarity_tetra = as.numeric(sim_tetra),
            similarity_control = as.numeric(sim_control),
            delta_connectivity_similarity = as.numeric(sim_tetra - sim_control),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
bootstrap_df <- if (length(bootstrap_rows) > 0) do.call(rbind, bootstrap_rows) else data.frame()
write_tsv(bootstrap_df, file.path(out_dir, "module_conservation_bootstrap.tsv"))

permutation_rows <- list()
universe <- rownames(mat)
if (perm_n > 0 && length(target_modules) > 0) {
  for (module_label in target_modules) {
    module_gene_count <- length(intersect(module_genes[[module_label]], universe))
    genes_used <- min(module_gene_count, max_resample_genes)
    observed_genes <- sample_genes(module_genes[[module_label]], max_resample_genes)
    expr_obs <- mat[observed_genes, , drop = FALSE]
    perm_group_names <- unique(c("CSP", tetra_groups, contrast_controls))
    k_obs <- lapply(group_samples[perm_group_names], function(x) cor_to_connectivity(expr_obs, x, power)$k)
    observed_delta <- list()
    for (tetra_group in tetra_groups) {
      for (control_group in contrast_controls) {
        key <- paste(tetra_group, control_group, sep = "__")
        sim_tetra <- rank_cor(k_obs[["CSP"]], k_obs[[tetra_group]])["rho"]
        sim_control <- rank_cor(k_obs[["CSP"]], k_obs[[control_group]])["rho"]
        observed_delta[[key]] <- as.numeric(sim_tetra - sim_control)
      }
    }

    nulls <- lapply(observed_delta, function(x) numeric(perm_n))
    for (iter in seq_len(perm_n)) {
      rand_genes <- sample(universe, genes_used, replace = FALSE)
      expr_rand <- mat[rand_genes, , drop = FALSE]
      k_rand <- lapply(group_samples[perm_group_names], function(x) cor_to_connectivity(expr_rand, x, power)$k)
      for (tetra_group in tetra_groups) {
        for (control_group in contrast_controls) {
          key <- paste(tetra_group, control_group, sep = "__")
          sim_tetra <- rank_cor(k_rand[["CSP"]], k_rand[[tetra_group]])["rho"]
          sim_control <- rank_cor(k_rand[["CSP"]], k_rand[[control_group]])["rho"]
          nulls[[key]][iter] <- as.numeric(sim_tetra - sim_control)
        }
      }
    }

    for (key in names(observed_delta)) {
      bits <- strsplit(key, "__", fixed = TRUE)[[1]]
      obs <- observed_delta[[key]]
      null_delta <- nulls[[key]]
      p_greater <- (1 + sum(null_delta >= obs, na.rm = TRUE)) / (1 + sum(is.finite(null_delta)))
      p_two_sided <- (1 + sum(abs(null_delta) >= abs(obs), na.rm = TRUE)) / (1 + sum(is.finite(null_delta)))
      permutation_rows[[length(permutation_rows) + 1]] <- data.frame(
        module_label = module_label,
        module_gene_count = module_gene_count,
        genes_used = genes_used,
        tetra_group = bits[1],
        control_group = bits[2],
        observed_delta = obs,
        null_median = stats::median(null_delta, na.rm = TRUE),
        null_q025 = stats::quantile(null_delta, 0.025, na.rm = TRUE),
        null_q975 = stats::quantile(null_delta, 0.975, na.rm = TRUE),
        p_greater = p_greater,
        p_two_sided = p_two_sided,
        stringsAsFactors = FALSE
      )
    }
  }
}
perm_df <- if (length(permutation_rows) > 0) do.call(rbind, permutation_rows) else data.frame()
if (nrow(perm_df) > 0) {
  perm_df$fdr_greater <- stats::p.adjust(perm_df$p_greater, method = "BH")
  perm_df$fdr_two_sided <- stats::p.adjust(perm_df$p_two_sided, method = "BH")
}
write_tsv(perm_df, file.path(out_dir, "module_conservation_permutation.tsv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 9))

  # Only the manuscript heatmap is exported in this publication repository.
  stage_heat_rows_df <- stage_similarity_df[stage_similarity_df$stage != "Bm3235c_carposporophyte", , drop = FALSE]
  stage_heat_plot_order <- stage_order[stage_order != "Bm3235c_carposporophyte"]
  stage_heat_rows_df$module_label <- factor(stage_heat_rows_df$module_label, levels = rev(names(module_genes)))
  stage_heat_rows_df$stage <- factor(stage_heat_rows_df$stage, levels = stage_heat_plot_order)
  stage_heat_rows_df$rho_label <- ifelse(is.finite(stage_heat_rows_df$similarity_stage), sprintf("%.2f", stage_heat_rows_df$similarity_stage), "")
  stage_heat_rows_df$rho_text_color <- ifelse(is.finite(stage_heat_rows_df$similarity_stage) & abs(stage_heat_rows_df$similarity_stage) >= 0.55, "white", "black")
  p_stage_rows_heat <- ggplot2::ggplot(stage_heat_rows_df, ggplot2::aes(stage, module_label, fill = similarity_stage)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = rho_label, color = rho_text_color), size = 1.55) +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::scale_fill_gradientn(
      colours = c("#4d4d4d", "#9e9e9e", "#d0d0d0", "#ffffff", "#ead8ba", "#c99f70", "#a95555", "#742323"),
      values = scales::rescale(c(-1, -0.66, -0.33, 0, 0.25, 0.5, 0.7, 1)),
      limits = c(-1, 1), name = "Spearman rho"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "Connectivity similarity to carposporophyte by individual stage") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(file.path(out_dir, "module_conservation_similarity_heatmap_12_stages_no_csp_self.pdf"), p_stage_rows_heat, width = 8.2, height = 7.2)

}

cat(sprintf(
  "Module conservation analysis completed: %d modules, %d CSP-associated modules, boot_n=%d, perm_n=%d, max_resample_genes=%d\n",
  length(module_genes), length(target_modules), boot_n, perm_n, max_resample_genes
))
