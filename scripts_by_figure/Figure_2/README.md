# Figure_2

Manuscript outputs and generating scripts. Run from this directory with `./run.sh`, or run individual commands from the repository root.

| Panel/Table | Output | Generating script(s) |
| --- | --- | --- |
| H | `../../manuscript_outputs/Figure_2/H_qc_pca_stage_generation.pdf` | `scripts/26_plot_qc_pca_and_correlation.R` |
| I | `../../manuscript_outputs/Figure_2/I_02b_ridge_log2aFC_by_tissue.pdf` | `scripts/48_ase_snp_count_pipeline.sh; scripts/49_ase_expression_analysis.R` |
| J | `../../manuscript_outputs/Figure_2/J_genome_expression_de_stacked_barplot.pdf` | `scripts/32_summarize_de_stability.R; scripts/33_plot_genome_expression_de_summary.R` |
| K | `../../manuscript_outputs/Figure_2/K_expr_detectability_stage_unique_barplot.pdf` | `scripts/24_plot_stage_unique_detectability.R` |
| L | `../../manuscript_outputs/Figure_2/L_de_any_zscore_heatmap.pdf` | `scripts/21_plot_de_any_zscore_heatmap.R` |
| M | `../../manuscript_outputs/Figure_2/M_de_left_right_bar_summary.pdf` | `scripts/20_plot_de_left_right_summary.R` |

## Dependency Note

Panel I runs the ASE shell workflow before the R plotting step. That workflow expects the cluster tools and upstream `nf-core/sarek` joint VCF described in `../../docs/dependencies_and_provenance.md`.
