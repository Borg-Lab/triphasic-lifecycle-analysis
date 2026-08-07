#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

scripts/run_rscript.sh scripts/26_plot_qc_pca_and_correlation.R
scripts/run_rscript.sh scripts/27_plot_gene_body_coverage.R
scripts/run_rscript.sh scripts/28_compare_smartseq3_nebnext.R
