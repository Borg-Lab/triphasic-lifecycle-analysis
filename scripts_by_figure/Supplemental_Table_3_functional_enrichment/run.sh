#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

scripts/run_rscript.sh scripts/16_module_functional_enrichment.R
scripts/run_rscript.sh scripts/15_marker_functional_enrichment.R
scripts/run_rscript.sh scripts/25_stage_unique_functional_enrichment.R
scripts/run_rscript.sh scripts/37_define_triphasy_generation_sets.R
scripts/run_rscript.sh scripts/43_analyze_sex_bias_and_uv_expression.R
