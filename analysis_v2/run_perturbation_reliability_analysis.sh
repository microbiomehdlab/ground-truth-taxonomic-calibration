#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${TRANSITION_LEDGER:?Set TRANSITION_LEDGER to a sealed target-excluded transition ledger}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new reliability-analysis directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
[[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" || "$ANALYSIS_STATUS" == "DEFINITIVE" ]] || {
  echo "[ERROR] invalid ANALYSIS_STATUS" >&2; exit 1;
}
[[ -s "$TRANSITION_LEDGER" ]] || { echo "[ERROR] missing transition ledger" >&2; exit 1; }
[[ -s "$ANALYSIS_SIF" ]] || { echo "[ERROR] missing analysis image" >&2; exit 1; }
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty" >&2; exit 1
fi
if [[ "$ANALYSIS_STATUS" == "DEFINITIVE" ]]; then
  echo "[ERROR] thresholds and held-out validation are not frozen; definitive status is disabled" >&2
  exit 1
fi
score_args=(--transition-ledger "$TRANSITION_LEDGER" --outdir "$OUTDIR/scores")
if [[ -n "${ARTIFACT_SCORES:-}" ]]; then
  [[ -s "$ARTIFACT_SCORES" ]] || { echo "[ERROR] missing ARTIFACT_SCORES" >&2; exit 1; }
  score_args+=(--artifact-scores "$ARTIFACT_SCORES")
fi
if [[ -n "${DISEASE_RESULTS:-}" ]]; then
  [[ -s "$DISEASE_RESULTS" ]] || { echo "[ERROR] missing DISEASE_RESULTS" >&2; exit 1; }
  score_args+=(--disease-results "$DISEASE_RESULTS")
fi
python3 analysis_v2/tests/test_perturbation_reliability_scores.py
python3 analysis_v2/scripts/build_perturbation_reliability_scores.py "${score_args[@]}"
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_perturbation_reliability_figures.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_perturbation_reliability_scores.R \
  --scores "$OUTDIR/scores/taxon_reliability_scores.tsv" --outdir "$OUTDIR/report"
for path in "$OUTDIR/scores/SUCCESS" "$OUTDIR/report/SUCCESS" \
            "$OUTDIR/report/figures/reliability_landscape.pdf" \
            "$OUTDIR/report/figures/reliability_by_replication.pdf"; do
  [[ -s "$path" ]] || { echo "[ERROR] missing output: $path" >&2; exit 1; }
done
printf 'status=DEVELOPMENT_ONLY\nautomatic_feature_removal=NO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
printf 'status=PASS\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Perturbation reliability analysis: $OUTDIR"
