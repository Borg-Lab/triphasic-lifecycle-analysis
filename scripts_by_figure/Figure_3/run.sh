#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

scripts/run_rscript.sh scripts/38_plot_triphasy_expression_profiles.R
scripts/run_rscript.sh scripts/37_define_triphasy_generation_sets.R
scripts/run_rscript.sh scripts/42_plot_triphasy_overlap.R
scripts/run_rscript.sh scripts/41_plot_triphasy_chromosomal_clustering.R
scripts/run_rscript.sh scripts/40_triphasy_tap_enrichment.R
scripts/run_rscript.sh scripts/43_analyze_sex_bias_and_uv_expression.R
scripts/run_rscript.sh scripts/39_plot_triphasy_phylorank_enrichment.R
