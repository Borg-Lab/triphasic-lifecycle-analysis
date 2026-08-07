#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

scripts/run_rscript.sh scripts/26_plot_qc_pca_and_correlation.R
bash scripts/48_ase_snp_count_pipeline.sh
scripts/run_rscript.sh scripts/49_ase_expression_analysis.R
scripts/run_rscript.sh scripts/32_summarize_de_stability.R
scripts/run_rscript.sh scripts/33_plot_genome_expression_de_summary.R
scripts/run_rscript.sh scripts/24_plot_stage_unique_detectability.R
scripts/run_rscript.sh scripts/21_plot_de_any_zscore_heatmap.R
scripts/run_rscript.sh scripts/20_plot_de_left_right_summary.R
