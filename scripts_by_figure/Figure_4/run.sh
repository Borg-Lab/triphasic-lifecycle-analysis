#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

scripts/run_rscript.sh scripts/22_plot_marker_expression_heatmap.R
scripts/run_rscript.sh scripts/09_plot_module_stage_enrichment.R
scripts/run_rscript.sh scripts/10_plot_module_stage_linepanels.R
scripts/run_rscript.sh scripts/14_plot_module_chromosome_composite.R
scripts/run_rscript.sh scripts/13_plot_module_chromosome_contribution.R
scripts/run_rscript.sh scripts/35_plot_marker_module_phylorank_stacks.R
scripts/run_rscript.sh scripts/36_plot_module_phylorank_enrichment.R
