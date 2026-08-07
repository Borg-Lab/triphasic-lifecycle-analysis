#!/usr/bin/env bash
set -euo pipefail

# Optional Singularity runtime. If TRIPHASIC_R_IMAGE is unset, use Rscript from PATH.
TRIPHASIC_R_IMAGE="${TRIPHASIC_R_IMAGE:-}"
TRIPHASIC_R_LIBS="${TRIPHASIC_R_LIBS:-}"
SINGULARITY_BIND="${SINGULARITY_BIND:-$(pwd)}"

if [[ $# -lt 1 ]]; then
  echo "Usage: scripts/run_rscript.sh <script.R| -e 'expr'> [args...]" >&2
  exit 1
fi

run_rscript() {
  if [[ -n "$TRIPHASIC_R_IMAGE" ]]; then
    local env_args=()
    local prefix=""
    if [[ -n "$TRIPHASIC_R_LIBS" ]]; then
      env_args=(--env "R_LIBS_USER=${TRIPHASIC_R_LIBS}")
      prefix=".libPaths(c('${TRIPHASIC_R_LIBS}', .libPaths())); "
    fi
    singularity exec --bind "$SINGULARITY_BIND" "${env_args[@]}" "$TRIPHASIC_R_IMAGE" \
      Rscript -e "${prefix}${1}" "${@:2}"
  else
    Rscript -e "$1" "${@:2}"
  fi
}

if [[ $# -ge 2 && "$1" == "-e" ]]; then
  EXPR="$2"
  shift 2
  run_rscript "$EXPR" "$@"
elif [[ $# -ge 1 && "$1" == *.R && -f "$1" ]]; then
  SCRIPT_PATH="$1"
  shift
  run_rscript "source('${SCRIPT_PATH}')" "$@"
else
  if [[ -n "$TRIPHASIC_R_IMAGE" ]]; then
    env_args=()
    [[ -n "$TRIPHASIC_R_LIBS" ]] && env_args=(--env "R_LIBS_USER=${TRIPHASIC_R_LIBS}")
    singularity exec --bind "$SINGULARITY_BIND" "${env_args[@]}" "$TRIPHASIC_R_IMAGE" Rscript "$@"
  else
    Rscript "$@"
  fi
fi
