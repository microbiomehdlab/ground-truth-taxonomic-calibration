#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${CANONICAL_INPUT:?Set CANONICAL_INPUT to a validated canonical_input.tsv}"
: "${SAMPLE_METADATA:?Set SAMPLE_METADATA to the frozen sample manifest with age, sex, and bmi}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty disease-biomarker directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
ALIASES="${ALIASES:-examples/spike_taxon_aliases.csv}"
SPIKE_PANEL="${SPIKE_PANEL:-spikes/spike_panel.tsv}"
CANONICAL_VALIDATION_SUCCESS="${CANONICAL_VALIDATION_SUCCESS:-$(dirname "$CANONICAL_INPUT")/validation/SUCCESS}"

[[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" || "$ANALYSIS_STATUS" == "DEFINITIVE" ]] || {
  echo "[ERROR] ANALYSIS_STATUS must be DEVELOPMENT_ONLY or DEFINITIVE" >&2; exit 1;
}
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty: $OUTDIR" >&2; exit 1
fi
for path in "$CANONICAL_INPUT" "$CANONICAL_VALIDATION_SUCCESS" "$SAMPLE_METADATA" \
            "$ANALYSIS_SIF" "$ALIASES" "$SPIKE_PANEL"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
mkdir -p "$OUTDIR"/{input,models,evaluation,provenance}
{
  printf 'field\tvalue\n'
  printf 'status\tIN_PROGRESS\n'
  printf 'analysis_status\t%s\n' "$ANALYSIS_STATUS"
  printf 'created_at\t%s\n' "$(date -Iseconds)"
  printf 'repository_commit\t%s\n' "$(git rev-parse HEAD)"
  printf 'canonical_input\t%s\n' "$(realpath "$CANONICAL_INPUT")"
  printf 'sample_metadata\t%s\n' "$(realpath "$SAMPLE_METADATA")"
  printf 'analysis_sif\t%s\n' "$(realpath "$ANALYSIS_SIF")"
} > "$OUTDIR/provenance/run_manifest.tsv"
if [[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" ]]; then
  printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
fi

python3 analysis_v2/tests/test_biomarker_abundance_input.py
python3 analysis_v2/scripts/build_biomarker_abundance_input.py \
  --canonical "$CANONICAL_INPUT" --aliases "$ALIASES" \
  --sample-metadata "$SAMPLE_METADATA" --outdir "$OUTDIR/input"

apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_disease_biomarker_models.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/fit_disease_biomarker_models.R \
  --profile-manifest "$OUTDIR/input/biomarker_profile_manifest.tsv" \
  --abundance-long "$OUTDIR/input/biomarker_abundance_long.tsv" \
  --outdir "$OUTDIR/models"

python3 analysis_v2/tests/test_biomarker_propagation.py
python3 analysis_v2/scripts/evaluate_biomarker_propagation.py \
  --calls "$OUTDIR/models/primary_disease_da_results.tsv" --aliases "$ALIASES" \
  --spike-panel "$SPIKE_PANEL" --outdir "$OUTDIR/evaluation" \
  --q-thresholds 0.05,0.10

sha256sum "$CANONICAL_INPUT" "$SAMPLE_METADATA" "$ANALYSIS_SIF" "$ALIASES" "$SPIKE_PANEL" \
  "$OUTDIR/input/biomarker_profile_manifest.tsv" \
  "$OUTDIR/input/biomarker_abundance_long.tsv" \
  "$OUTDIR/models/primary_disease_da_results.tsv" \
  "$OUTDIR/evaluation/biomarker_propagation_metrics.tsv" \
  > "$OUTDIR/provenance/run_inputs_and_primary_outputs.sha256"
sed -i 's/^status\tIN_PROGRESS$/status\tPASS/' "$OUTDIR/provenance/run_manifest.tsv"
printf 'analysis\tnative_disease_biomarker_propagation\nanalysis_status\t%s\nstatus\tPASS\n' \
  "$ANALYSIS_STATUS" > "$OUTDIR/SUCCESS"
echo "[PASS] Sealed disease-biomarker propagation analysis: $OUTDIR"
