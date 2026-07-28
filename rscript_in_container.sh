#!/usr/bin/env bash
set -euo pipefail
PROJECT="${PROJECT:-$PWD}"
: "${SIF:?Set SIF to the MaAsLin2 analysis container image}"

apptainer exec --cleanenv \
  --bind "$PROJECT:$PROJECT" \
  "$SIF" \
  Rscript "$@"
