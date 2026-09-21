#!/usr/bin/env bash
# Re-evaluate the frozen partial-snapshot disease calls with the corrected
# target-exclusion policy. Exploratory figure previews only; no model refit.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${SOURCE_RUN:?Set SOURCE_RUN to the sealed partial three-cohort map-reduce run}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the validated analysis image}"
: "${OUTDIR:?Set OUTDIR to a new development-only directory}"
CALLS="$SOURCE_RUN/models/disease/models/primary_disease_da_results.tsv"
for path in "$SOURCE_RUN/SUCCESS" "$SOURCE_RUN/DEVELOPMENT_ONLY.txt" \
            "$SOURCE_RUN/models/disease/SUCCESS" "$CALLS" "$ANALYSIS_SIF" \
            examples/spike_taxon_aliases.csv spikes/spike_panel.tsv; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing source: $path" >&2; exit 1; }
done
[[ ! -e "$OUTDIR" ]] || { echo "[ERROR] OUTDIR exists: $OUTDIR" >&2; exit 1; }
mkdir -p "$OUTDIR"
printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
apptainer exec --cleanenv --bind "$ROOT:$ROOT" --pwd "$ROOT" "$ANALYSIS_SIF" \
  python3 analysis_v2/scripts/evaluate_disease_biomarker_propagation.py \
  --calls "$CALLS" --aliases examples/spike_taxon_aliases.csv \
  --spike-panel spikes/spike_panel.tsv --q-thresholds 0.05,0.10 \
  --outdir "$OUTDIR/evaluation"
test -s "$OUTDIR/evaluation/SUCCESS"
apptainer exec --cleanenv --bind "$ROOT:$ROOT" --pwd "$ROOT" "$ANALYSIS_SIF" \
  python3 analysis_v2/scripts/build_partial_replication_matrix.py \
  --calls "$CALLS" --outdir "$OUTDIR/replication"
test -s "$OUTDIR/replication/SUCCESS"
apptainer exec --cleanenv --bind "$ROOT:$ROOT" --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_disease_biomarker_robustness.R \
  --metrics "$OUTDIR/evaluation/disease_biomarker_propagation_metrics.tsv" \
  --ledger "$OUTDIR/evaluation/disease_biomarker_transition_ledger.tsv" \
  --analysis-status DEVELOPMENT_ONLY --outdir "$OUTDIR/robustness"
for path in "$OUTDIR/robustness/SUCCESS" \
            "$OUTDIR/replication/directional_replication_matrix.svg" \
            "$OUTDIR/robustness/figures/bystander_biomarker_retention.png" \
            "$OUTDIR/robustness/figures/bystander_transition_effect_change.png"; do
  test -s "$path" || { echo "[ERROR] Missing output: $path" >&2; exit 1; }
done
sha256sum "$SOURCE_RUN/SUCCESS" "$SOURCE_RUN/DEVELOPMENT_ONLY.txt" \
  "$CALLS" "$ANALYSIS_SIF" > "$OUTDIR/run_inputs.sha256"
printf 'status\tPASS\nanalysis_status\tDEVELOPMENT_ONLY\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Partial-snapshot target-excluded biomarker figures: $OUTDIR"
