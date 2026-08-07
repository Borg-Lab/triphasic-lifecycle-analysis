# Chromosomal Clustering Method

## Objective
Test whether coexpressed gene sets are genomically clustered.

## Statistic
- For each set, compute nearest-neighbor distance using gene midpoints.
- Primary statistic: median nearest-neighbor distance (smaller = tighter clustering).

## Null model
- Chromosome-stratified permutation with matched per-chromosome gene counts.
- Empirical p-value: `(count(perm <= observed) + 1) / (N + 1)`.

## Multiple testing and calls
- BH-adjusted FDR across tested sets.
- `cluster_flag = TRUE` when `fdr <= alpha` and `z_score < 0`.

## Downstream enrichment
- Fisher's exact test between clustered sets and stage/generation markers.
- BH correction within enrichment family.

## Current Result Tables

Regenerate the combined chromosomal clustering exports with `scripts/47_export_chromosomal_clustering_tables.R`.

- Combined chromosomal clustering statistics are in `runs/triphasic_manuscript/results/chrom_cluster/combined_chrom_cluster_enrichment_stats.tsv`.
  - This combines modules, triphasy generation-biased/repressed sets, and sex-biased gene sets into one normalized table.
  - Key columns are `set_type`, `module_rank`, `set_label`, `n_genes`, `obs_median_nn_dist`, `perm_median`, `z_score`, `empirical_p`, `fdr`, `cluster_flag`, `n_proximity_blocks`, `n_genes_in_proximity_blocks`, `pct_genes_in_proximity_blocks`, and `clustering_strength_log2ratio`.
- Module proximity-block membership with current renamed module labels is in `runs/triphasic_manuscript/results/chrom_cluster/chrom_proximity_blocks_labeled.tsv`.
  - Each row is one numbered genomic cluster/block for one module, with `cluster_id`, coordinates, `n_genes`, and semicolon-separated `gene_ids`.
- A combined proximity-block lookup across modules, triphasy sets, and sex-biased gene sets is in `runs/triphasic_manuscript/results/chrom_cluster/all_chrom_proximity_blocks_labeled.tsv`.
  - Use this when a figure or manuscript table needs the actual genes contained in each cluster.
- Original set-specific block tables remain in `generation/triphasy_up_down_chrom_proximity_blocks.tsv` and `sex_bias/sex_biased_chrom_proximity_blocks.tsv`.

## Current Plotting Notes

- Module chromosome contribution is reported in `runs/triphasic_manuscript/results/chrom_cluster/module_chromosome_contribution.tsv`.
- The heatmap `module_chromosome_contribution_heatmap.pdf` uses the current dendrogram/order-aware module labels from `module_label_map.tsv`.
- Its fill scale is capped at `>=12%` and uses continuous breaks at 0, 1, 3, 6, and 12%, because the current distribution is concentrated below 10% with rare high outliers and no values in the old 20-40% bin.
- Triphasy and sex-biased chromosome contribution heatmaps are separate outputs and should be regenerated from their own current set definitions when those sets change.
