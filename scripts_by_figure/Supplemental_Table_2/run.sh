#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
export PATHS_CONFIG=${PATHS_CONFIG:-config/paths.yml}
export ANALYSIS_CONFIG=${ANALYSIS_CONFIG:-config/analysis.yml}
export DESIGN_CONFIG=${DESIGN_CONFIG:-config/design.tsv}
export CONTRASTS_CONFIG=${CONTRASTS_CONFIG:-config/contrasts.tsv}

scripts/run_rscript.sh scripts/46_export_annotated_group_tpm_means.R
