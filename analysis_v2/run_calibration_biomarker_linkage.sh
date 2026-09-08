#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${ENDPOINTS_FILE:?Set ENDPOINTS_FILE to a sealed paired_endpoints.tsv}"
: "${PAIRED_RUN:?Set PAIRED_RUN to a sealed paired-biomarker analysis}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty linkage directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
ENDPOINTS_SUCCESS="${ENDPOINTS_SUCCESS:-$(dirname "$ENDPOINTS_FILE")/SUCCESS}"
[[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" || "$ANALYSIS_STATUS" == "DEFINITIVE" ]] || { echo "[ERROR] Invalid ANALYSIS_STATUS" >&2; exit 1; }
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then echo "[ERROR] OUTDIR must be new or empty" >&2; exit 1; fi
for path in "$ENDPOINTS_FILE" "$ENDPOINTS_SUCCESS" "$PAIRED_RUN/SUCCESS" "$PAIRED_RUN/evaluation/SUCCESS" \
  "$PAIRED_RUN/evaluation/biomarker_propagation_metrics.tsv" "$ANALYSIS_SIF"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
secondary_args=()
if [[ -n "${SECONDARY_PAIRED_RUN:-}" ]]; then
  secondary_metrics="$SECONDARY_PAIRED_RUN/evaluation/biomarker_propagation_metrics.tsv"
  for path in "$SECONDARY_PAIRED_RUN/SUCCESS" "$SECONDARY_PAIRED_RUN/evaluation/SUCCESS" "$secondary_metrics"; do
    [[ -s "$path" ]] || { echo "[ERROR] Missing secondary input: $path" >&2; exit 1; }
  done
  secondary_args=(--secondary-biomarker-metrics "$secondary_metrics")
fi
if [[ "$ANALYSIS_STATUS" == "DEFINITIVE" ]] && \
   { [[ -e "$PAIRED_RUN/DEVELOPMENT_ONLY.txt" ]] || [[ -e "$(dirname "$ENDPOINTS_FILE")/DEVELOPMENT_ONLY.txt" ]]; }; then
  echo "[ERROR] Development inputs cannot generate a definitive linkage package" >&2; exit 1
fi
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_calibration_biomarker_linkage.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/link_calibration_to_biomarkers.R \
  --endpoints "$ENDPOINTS_FILE" \
  --biomarker-metrics "$PAIRED_RUN/evaluation/biomarker_propagation_metrics.tsv" \
  "${secondary_args[@]}" \
  --outdir "$OUTDIR" --analysis-status "$ANALYSIS_STATUS"
echo "[PASS] Sealed calibration-to-biomarker linkage: $OUTDIR"
