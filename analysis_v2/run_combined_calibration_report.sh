#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${MODELS_ROOT:?Set MODELS_ROOT to the directory containing detection_* and continuous_* models}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty report directory}"
: "${REPORT_STATUS:?Set REPORT_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
[[ "$REPORT_STATUS" == "DEVELOPMENT_ONLY" || "$REPORT_STATUS" == "DEFINITIVE" ]] || { echo "[ERROR] Invalid REPORT_STATUS" >&2; exit 1; }
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty" >&2
  exit 1
fi
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_combined_calibration_report.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/make_combined_calibration_report.R \
  --models-root "$MODELS_ROOT" --outdir "$OUTDIR" --report-status "$REPORT_STATUS"
echo "[PASS] Sealed combined calibration reporting package: $OUTDIR"
