#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${CANONICAL_INPUT:?Set CANONICAL_INPUT to a validated canonical_input.tsv}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty biomarker-analysis directory}"
ALIASES="${ALIASES:-examples/spike_taxon_aliases.csv}"
SPIKE_PANEL="${SPIKE_PANEL:-spikes/spike_panel.tsv}"

if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty: $OUTDIR" >&2
  exit 1
fi
for path in "$CANONICAL_INPUT" "$ANALYSIS_SIF" "$ALIASES" "$SPIKE_PANEL"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
mkdir -p "$OUTDIR"/{input,models,evaluation,provenance}
{
  printf 'field\tvalue\n'
  printf 'status\tIN_PROGRESS\n'
  printf 'created_at\t%s\n' "$(date -Iseconds)"
  printf 'repository_commit\t%s\n' "$(git rev-parse HEAD)"
  printf 'canonical_input\t%s\n' "$(realpath "$CANONICAL_INPUT")"
  printf 'analysis_sif\t%s\n' "$(realpath "$ANALYSIS_SIF")"
  printf 'primary_q_threshold\t0.05\n'
  printf 'sensitivity_q_threshold\t0.10\n'
} > "$OUTDIR/provenance/run_manifest.tsv"

python3 analysis_v2/tests/test_biomarker_abundance_input.py
python3 analysis_v2/scripts/build_biomarker_abundance_input.py \
  --canonical "$CANONICAL_INPUT" --aliases "$ALIASES" --outdir "$OUTDIR/input"

apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_paired_biomarker_models.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/fit_paired_biomarker_models.R \
  --profile-manifest "$OUTDIR/input/biomarker_profile_manifest.tsv" \
  --abundance-long "$OUTDIR/input/biomarker_abundance_long.tsv" \
  --outdir "$OUTDIR/models"

python3 analysis_v2/tests/test_biomarker_propagation.py
python3 analysis_v2/scripts/evaluate_biomarker_propagation.py \
  --calls "$OUTDIR/models/paired_da_results.tsv" --aliases "$ALIASES" \
  --spike-panel "$SPIKE_PANEL" --outdir "$OUTDIR/evaluation" \
  --q-thresholds 0.05,0.10

sha256sum "$CANONICAL_INPUT" "$ANALYSIS_SIF" "$ALIASES" "$SPIKE_PANEL" \
  "$OUTDIR/input/biomarker_profile_manifest.tsv" \
  "$OUTDIR/input/biomarker_abundance_long.tsv" \
  "$OUTDIR/models/paired_da_results.tsv" \
  "$OUTDIR/evaluation/biomarker_propagation_metrics.tsv" \
  > "$OUTDIR/provenance/run_inputs_and_primary_outputs.sha256"
sed -i 's/^status\tIN_PROGRESS$/status\tPASS/' "$OUTDIR/provenance/run_manifest.tsv"
printf 'analysis\tpaired_biomarker_propagation\nstatus\tPASS\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Sealed paired biomarker-propagation analysis: $OUTDIR"
