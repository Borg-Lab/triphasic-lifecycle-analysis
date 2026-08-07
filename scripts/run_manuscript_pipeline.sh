#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

echo '[RUN ] 01_validate_inputs.R'
scripts/run_rscript.sh scripts/01_validate_inputs.R
echo '[DONE] 01_validate_inputs.R'

echo '[RUN ] 02_prepare_expression_matrices.R'
scripts/run_rscript.sh scripts/02_prepare_expression_matrices.R
echo '[DONE] 02_prepare_expression_matrices.R'

echo '[RUN ] 03_run_differential_expression.R'
scripts/run_rscript.sh scripts/03_run_differential_expression.R
echo '[DONE] 03_run_differential_expression.R'

echo '[RUN ] 04_define_markers_and_specificity.R'
scripts/run_rscript.sh scripts/04_define_markers_and_specificity.R
echo '[DONE] 04_define_markers_and_specificity.R'

echo '[RUN ] 05_build_wgcna_modules.R'
scripts/run_rscript.sh scripts/05_build_wgcna_modules.R
echo '[DONE] 05_build_wgcna_modules.R'

echo '[RUN ] 06_build_module_label_map.R'
scripts/run_rscript.sh scripts/06_build_module_label_map.R
echo '[DONE] 06_build_module_label_map.R'

echo '[RUN ] 07_build_module_gene_counts.R'
scripts/run_rscript.sh scripts/07_build_module_gene_counts.R
echo '[DONE] 07_build_module_gene_counts.R'

echo '[RUN ] 08_annotate_module_membership.R'
scripts/run_rscript.sh scripts/08_annotate_module_membership.R
echo '[DONE] 08_annotate_module_membership.R'

echo '[RUN ] 09_plot_module_stage_enrichment.R'
scripts/run_rscript.sh scripts/09_plot_module_stage_enrichment.R
echo '[DONE] 09_plot_module_stage_enrichment.R'

echo '[RUN ] 10_plot_module_stage_linepanels.R'
scripts/run_rscript.sh scripts/10_plot_module_stage_linepanels.R
echo '[DONE] 10_plot_module_stage_linepanels.R'

echo '[RUN ] 11_module_conservation_across_stages.R'
scripts/run_rscript.sh scripts/11_module_conservation_across_stages.R
echo '[DONE] 11_module_conservation_across_stages.R'

echo '[RUN ] 12_run_chromosomal_clustering.R'
scripts/run_rscript.sh scripts/12_run_chromosomal_clustering.R
echo '[DONE] 12_run_chromosomal_clustering.R'

echo '[RUN ] 13_plot_module_chromosome_contribution.R'
scripts/run_rscript.sh scripts/13_plot_module_chromosome_contribution.R
echo '[DONE] 13_plot_module_chromosome_contribution.R'

echo '[RUN ] 14_plot_module_chromosome_composite.R'
scripts/run_rscript.sh scripts/14_plot_module_chromosome_composite.R
echo '[DONE] 14_plot_module_chromosome_composite.R'

echo '[RUN ] 15_marker_functional_enrichment.R'
scripts/run_rscript.sh scripts/15_marker_functional_enrichment.R
echo '[DONE] 15_marker_functional_enrichment.R'

echo '[RUN ] 16_module_functional_enrichment.R'
scripts/run_rscript.sh scripts/16_module_functional_enrichment.R
echo '[DONE] 16_module_functional_enrichment.R'

echo '[RUN ] 17_analyze_module_tf_associations.R'
scripts/run_rscript.sh scripts/17_analyze_module_tf_associations.R
echo '[DONE] 17_analyze_module_tf_associations.R'

echo '[RUN ] 18_plot_module_tf_lollipop.R'
scripts/run_rscript.sh scripts/18_plot_module_tf_lollipop.R
echo '[DONE] 18_plot_module_tf_lollipop.R'

echo '[RUN ] 19_plot_de_volcano_panels.R'
scripts/run_rscript.sh scripts/19_plot_de_volcano_panels.R
echo '[DONE] 19_plot_de_volcano_panels.R'

echo '[RUN ] 20_plot_de_left_right_summary.R'
scripts/run_rscript.sh scripts/20_plot_de_left_right_summary.R
echo '[DONE] 20_plot_de_left_right_summary.R'

echo '[RUN ] 21_plot_de_any_zscore_heatmap.R'
scripts/run_rscript.sh scripts/21_plot_de_any_zscore_heatmap.R
echo '[DONE] 21_plot_de_any_zscore_heatmap.R'

echo '[RUN ] 22_plot_marker_expression_heatmap.R'
scripts/run_rscript.sh scripts/22_plot_marker_expression_heatmap.R
echo '[DONE] 22_plot_marker_expression_heatmap.R'

echo '[RUN ] 23_plot_detectability_barplot.R'
scripts/run_rscript.sh scripts/23_plot_detectability_barplot.R
echo '[DONE] 23_plot_detectability_barplot.R'

echo '[RUN ] 24_plot_stage_unique_detectability.R'
scripts/run_rscript.sh scripts/24_plot_stage_unique_detectability.R
echo '[DONE] 24_plot_stage_unique_detectability.R'

echo '[RUN ] 25_stage_unique_functional_enrichment.R'
scripts/run_rscript.sh scripts/25_stage_unique_functional_enrichment.R
echo '[DONE] 25_stage_unique_functional_enrichment.R'

echo '[RUN ] 26_plot_qc_pca_and_correlation.R'
scripts/run_rscript.sh scripts/26_plot_qc_pca_and_correlation.R
echo '[DONE] 26_plot_qc_pca_and_correlation.R'

echo '[RUN ] 27_plot_gene_body_coverage.R'
scripts/run_rscript.sh scripts/27_plot_gene_body_coverage.R
echo '[DONE] 27_plot_gene_body_coverage.R'

echo '[RUN ] 28_compare_smartseq3_nebnext.R'
scripts/run_rscript.sh scripts/28_compare_smartseq3_nebnext.R
echo '[DONE] 28_compare_smartseq3_nebnext.R'

echo '[RUN ] 29_compute_tai_signature.R'
scripts/run_rscript.sh scripts/29_compute_tai_signature.R
echo '[DONE] 29_compute_tai_signature.R'

echo '[RUN ] 30_plot_tai_signature_replots.R'
scripts/run_rscript.sh scripts/30_plot_tai_signature_replots.R
echo '[DONE] 30_plot_tai_signature_replots.R'

echo '[RUN ] 31_plot_phylorank_lifecycle_expression.R'
scripts/run_rscript.sh scripts/31_plot_phylorank_lifecycle_expression.R
echo '[DONE] 31_plot_phylorank_lifecycle_expression.R'

echo '[RUN ] 32_summarize_de_stability.R'
scripts/run_rscript.sh scripts/32_summarize_de_stability.R
echo '[DONE] 32_summarize_de_stability.R'

echo '[RUN ] 33_plot_genome_expression_de_summary.R'
scripts/run_rscript.sh scripts/33_plot_genome_expression_de_summary.R
echo '[DONE] 33_plot_genome_expression_de_summary.R'

echo '[RUN ] 34_plot_constitutive_dynamic_classes.R'
scripts/run_rscript.sh scripts/34_plot_constitutive_dynamic_classes.R
echo '[DONE] 34_plot_constitutive_dynamic_classes.R'

echo '[RUN ] 35_plot_marker_module_phylorank_stacks.R'
scripts/run_rscript.sh scripts/35_plot_marker_module_phylorank_stacks.R
echo '[DONE] 35_plot_marker_module_phylorank_stacks.R'

echo '[RUN ] 36_plot_module_phylorank_enrichment.R'
scripts/run_rscript.sh scripts/36_plot_module_phylorank_enrichment.R
echo '[DONE] 36_plot_module_phylorank_enrichment.R'

echo '[RUN ] 37_define_triphasy_generation_sets.R'
scripts/run_rscript.sh scripts/37_define_triphasy_generation_sets.R
echo '[DONE] 37_define_triphasy_generation_sets.R'

echo '[RUN ] 38_plot_triphasy_expression_profiles.R'
scripts/run_rscript.sh scripts/38_plot_triphasy_expression_profiles.R
echo '[DONE] 38_plot_triphasy_expression_profiles.R'

echo '[RUN ] 39_plot_triphasy_phylorank_enrichment.R'
scripts/run_rscript.sh scripts/39_plot_triphasy_phylorank_enrichment.R
echo '[DONE] 39_plot_triphasy_phylorank_enrichment.R'

echo '[RUN ] 40_triphasy_tap_enrichment.R'
scripts/run_rscript.sh scripts/40_triphasy_tap_enrichment.R
echo '[DONE] 40_triphasy_tap_enrichment.R'

echo '[RUN ] 41_plot_triphasy_chromosomal_clustering.R'
scripts/run_rscript.sh scripts/41_plot_triphasy_chromosomal_clustering.R
echo '[DONE] 41_plot_triphasy_chromosomal_clustering.R'

echo '[RUN ] 42_plot_triphasy_overlap.R'
scripts/run_rscript.sh scripts/42_plot_triphasy_overlap.R
echo '[DONE] 42_plot_triphasy_overlap.R'

echo '[RUN ] 43_analyze_sex_bias_and_uv_expression.R'
scripts/run_rscript.sh scripts/43_analyze_sex_bias_and_uv_expression.R
echo '[DONE] 43_analyze_sex_bias_and_uv_expression.R'

echo '[RUN ] 44_sex_biased_chromosomal_clustering.R'
scripts/run_rscript.sh scripts/44_sex_biased_chromosomal_clustering.R
echo '[DONE] 44_sex_biased_chromosomal_clustering.R'

echo '[RUN ] 45_de_normalization_robustness.R'
scripts/run_rscript.sh scripts/45_de_normalization_robustness.R
echo '[DONE] 45_de_normalization_robustness.R'

echo '[RUN ] 46_export_annotated_group_tpm_means.R'
scripts/run_rscript.sh scripts/46_export_annotated_group_tpm_means.R
echo '[DONE] 46_export_annotated_group_tpm_means.R'

echo '[RUN ] 47_export_chromosomal_clustering_tables.R'
scripts/run_rscript.sh scripts/47_export_chromosomal_clustering_tables.R
echo '[DONE] 47_export_chromosomal_clustering_tables.R'

echo '[RUN ] 48_ase_snp_count_pipeline.sh'
bash scripts/48_ase_snp_count_pipeline.sh
echo '[DONE] 48_ase_snp_count_pipeline.sh'

echo '[RUN ] 49_ase_expression_analysis.R'
scripts/run_rscript.sh scripts/49_ase_expression_analysis.R
echo '[DONE] 49_ase_expression_analysis.R'

