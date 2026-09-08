#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${PAIRED_RUN:?Set PAIRED_RUN to a sealed paired-biomarker analysis}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty report directory}"
: "${REPORT_STATUS:?Set REPORT_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
[[ "$REPORT_STATUS" == "DEVELOPMENT_ONLY" || "$REPORT_STATUS" == "DEFINITIVE" ]] || { echo "[ERROR] Invalid REPORT_STATUS" >&2; exit 1; }
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then echo "[ERROR] OUTDIR must be new or empty" >&2; exit 1; fi
for path in "$ANALYSIS_SIF" "$PAIRED_RUN/SUCCESS" "$PAIRED_RUN/models/SUCCESS" "$PAIRED_RUN/evaluation/SUCCESS"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_artificial_biomarker_report.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/make_artificial_biomarker_report.R \
  --paired-run "$PAIRED_RUN" --outdir "$OUTDIR" --report-status "$REPORT_STATUS"
echo "[PASS] Sealed artificial-biomarker reporting package: $OUTDIR"
