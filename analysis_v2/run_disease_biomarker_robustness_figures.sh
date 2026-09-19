#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${DISEASE_RUN:?Set DISEASE_RUN to a sealed disease-biomarker analysis}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty robustness-figure directory}"
: "${REPORT_STATUS:?Set REPORT_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"

[[ "$REPORT_STATUS" == "DEVELOPMENT_ONLY" || "$REPORT_STATUS" == "DEFINITIVE" ]] || {
  echo "[ERROR] REPORT_STATUS must be DEVELOPMENT_ONLY or DEFINITIVE" >&2
  exit 1
}
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty: $OUTDIR" >&2
  exit 1
fi

METRICS="$DISEASE_RUN/evaluation/disease_biomarker_propagation_metrics.tsv"
LEDGER="$DISEASE_RUN/evaluation/disease_biomarker_transition_ledger.tsv"
for path in "$ANALYSIS_SIF" "$DISEASE_RUN/SUCCESS" \
            "$DISEASE_RUN/evaluation/SUCCESS" "$METRICS" "$LEDGER"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
if [[ "$REPORT_STATUS" == "DEFINITIVE" && -e "$DISEASE_RUN/DEVELOPMENT_ONLY.txt" ]]; then
  echo "[ERROR] A development disease run cannot produce definitive figures" >&2
  exit 1
fi

apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_disease_biomarker_robustness_figures.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_disease_biomarker_robustness.R \
  --metrics "$METRICS" --ledger "$LEDGER" --outdir "$OUTDIR" \
  --analysis-status "$REPORT_STATUS"

for path in "$OUTDIR/SUCCESS" \
            "$OUTDIR/figures/bystander_biomarker_retention.pdf" \
            "$OUTDIR/figures/bystander_induced_call_rate.pdf" \
            "$OUTDIR/figures/bystander_gained_call_burden.pdf" \
            "$OUTDIR/figures/bystander_retention_atlas.pdf" \
            "$OUTDIR/figures/bystander_transition_effect_change.pdf" \
            "$OUTDIR/provenance/robustness_figures.sha256"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing output: $path" >&2; exit 1; }
done
echo "[PASS] Sealed disease-biomarker robustness figures: $OUTDIR"
