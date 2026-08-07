open_plot_device <- function(out_file, width_px = 1000, height_px = 700, res = 120) {
  ext <- tolower(tools::file_ext(out_file))
  if (identical(ext, "pdf")) {
    grDevices::pdf(out_file, width = width_px / res, height = height_px / res)
  } else {
    grDevices::png(out_file, width = width_px, height = height_px, res = res)
  }
}

plot_perm_distribution <- function(set_id, observed, perm_vals, out_file) {
  if (length(perm_vals) == 0) return(invisible(NULL))
  open_plot_device(out_file, width_px = 900, height_px = 600, res = 120)
  hist(perm_vals, breaks = 40, col = "grey80", border = "white", main = paste("Permutation NN distances:", set_id), xlab = "Median nearest-neighbor distance")
  abline(v = observed, col = "red", lwd = 2)
  dev.off()
}

plot_ggbio_proximity_blocks <- function(set_id, blocks_df, coords_df, out_file, bin_width = 100000) {
  if (!requireNamespace("ggbio", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    return(invisible(NULL))
  }
  sub_blocks <- blocks_df[blocks_df$set_id == set_id, , drop = FALSE]
  if (nrow(sub_blocks) == 0) return(invisible(NULL))

  chr_lengths_raw <- tapply(coords_df$end, coords_df$chr, max, na.rm = TRUE)
  chr_lengths <- stats::setNames(as.numeric(chr_lengths_raw), names(chr_lengths_raw))
  chr_lengths <- chr_lengths[!is.na(chr_lengths)]
  if (length(chr_lengths) == 0) return(invisible(NULL))

  block_gr <- GenomicRanges::GRanges(
    seqnames = sub_blocks$chr,
    ranges = IRanges::IRanges(start = sub_blocks$block_start, end = sub_blocks$block_end),
    set_id = sub_blocks$set_id,
    block_id = sub_blocks$block_id,
    n_genes = sub_blocks$n_genes
  )

  common_chr <- intersect(names(chr_lengths), as.character(GenomeInfoDb::seqlevels(block_gr)))
  if (length(common_chr) == 0) return(invisible(NULL))
  tmp_seqlens <- stats::setNames(
    rep(NA_real_, length(GenomeInfoDb::seqlevels(block_gr))),
    as.character(GenomeInfoDb::seqlevels(block_gr))
  )
  tmp_seqlens[common_chr] <- as.numeric(chr_lengths[common_chr])
  GenomeInfoDb::seqlengths(block_gr) <- tmp_seqlens
  bins <- GenomicRanges::tileGenome(chr_lengths, tilewidth = as.integer(bin_width), cut.last.tile.in.chrom = TRUE)
  counts <- GenomicRanges::countOverlaps(bins, block_gr)
  S4Vectors::mcols(bins)$counts <- as.numeric(counts)

  p <- ggbio::autoplot(
    bins,
    layout = "karyogram",
    ggplot2::aes(fill = counts),
    geom = "rect"
  ) +
    ggplot2::scale_fill_gradient(low = "white", high = "#b2182b") +
    ggplot2::ggtitle(paste0(set_id, " proximity blocks (", nrow(sub_blocks), " blocks)"))

  open_plot_device(out_file, width_px = 1200, height_px = 800, res = 140)
  print(p)
  grDevices::dev.off()
}

plot_cluster_heatmap <- function(summary_df, out_file) {
  if (nrow(summary_df) == 0) return(invisible(NULL))
  ord <- order(summary_df$z_score)
  z <- matrix(summary_df$z_score[ord], ncol = 1)
  rownames(z) <- summary_df$set_id[ord]
  colnames(z) <- "z_score"

  open_plot_device(out_file, width_px = 700, height_px = max(500, 30 * nrow(z)), res = 120)
  if (nrow(z) >= 2 && ncol(z) >= 2) {
    heatmap(
      z,
      Rowv = NA,
      Colv = NA,
      scale = "none",
      col = colorRampPalette(c("navy", "white", "firebrick"))(50),
      margins = c(8, 10)
    )
  } else {
    barplot(
      z[, 1],
      names.arg = rownames(z),
      las = 2,
      col = "steelblue",
      main = "Cluster z-scores",
      ylab = "z_score"
    )
    abline(h = 0, lty = 2, col = "grey40")
  }
  dev.off()
}

plot_enrichment_dotplot <- function(enrich_df, out_file) {
  if (nrow(enrich_df) == 0) return(invisible(NULL))
  open_plot_device(out_file, width_px = 1000, height_px = 700, res = 120)
  x <- seq_len(nrow(enrich_df))
  y <- -log10(pmax(enrich_df$fdr, 1e-300))
  cex <- pmax(0.8, pmin(3, enrich_df$odds_ratio))
  plot(x, y, pch = 19, cex = cex, xaxt = "n", xlab = "cluster:marker pairs", ylab = "-log10(FDR)", col = "steelblue")
  axis(1, at = x, labels = paste(enrich_df$cluster_id, enrich_df$marker_group, sep = ":"), las = 2, cex.axis = 0.65)
  abline(h = -log10(0.05), lty = 2, col = "red")
  dev.off()
}

plot_chrom_cluster_outputs <- function(chrom_cluster_obj, output_dir) {
  if (is.null(chrom_cluster_obj) || !is.list(chrom_cluster_obj)) return(character())
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out_files <- character()
  summary_df <- chrom_cluster_obj$summary
  blocks_df <- chrom_cluster_obj$blocks
  coords_df <- chrom_cluster_obj$coords
  nulls <- chrom_cluster_obj$nulls

  if (!is.null(summary_df) && nrow(summary_df) > 0) {
    heatmap_file <- file.path(output_dir, "chrom_cluster_zscore_heatmap.pdf")
    plot_cluster_heatmap(summary_df, heatmap_file)
    out_files <- c(out_files, heatmap_file)

    for (i in seq_len(nrow(summary_df))) {
      sid <- summary_df$set_id[i]
      obs <- summary_df$obs_median_nn_dist[i]
      perm_vals <- nulls[[sid]]
      if (is.null(perm_vals) || length(perm_vals) == 0) next
      perm_file <- file.path(output_dir, paste0("chrom_cluster_perm_", sid, ".pdf"))
      plot_perm_distribution(sid, obs, perm_vals, perm_file)
      out_files <- c(out_files, perm_file)

      if (!is.null(blocks_df) && !is.null(coords_df) && nrow(blocks_df) > 0 && nrow(coords_df) > 0) {
        block_file <- file.path(output_dir, paste0("chrom_cluster_blocks_", sid, ".pdf"))
        plot_ggbio_proximity_blocks(sid, blocks_df, coords_df, block_file)
        if (file.exists(block_file)) out_files <- c(out_files, block_file)
      }
    }
  }

  out_files[file.exists(out_files)]
}
