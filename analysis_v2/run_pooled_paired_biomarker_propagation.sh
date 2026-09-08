#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${CANONICAL_INPUT:?Set CANONICAL_INPUT}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF}"
: "${OUTDIR:?Set OUTDIR to a new directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS}"
ALIASES="${ALIASES:-examples/spike_taxon_aliases.csv}"; SPIKE_PANEL="${SPIKE_PANEL:-spikes/spike_panel.tsv}"
VALIDATION="${CANONICAL_VALIDATION_SUCCESS:-$(dirname "$CANONICAL_INPUT")/validation/SUCCESS}"
[[ "$ANALYSIS_STATUS" == DEVELOPMENT_ONLY || "$ANALYSIS_STATUS" == DEFINITIVE ]] || { echo '[ERROR] Invalid ANALYSIS_STATUS' >&2; exit 1; }
[[ ! -e "$OUTDIR" || -z "$(find "$OUTDIR" -mindepth 1 -print -quit)" ]] || { echo '[ERROR] OUTDIR must be empty' >&2; exit 1; }
for f in "$CANONICAL_INPUT" "$VALIDATION" "$ANALYSIS_SIF" "$ALIASES" "$SPIKE_PANEL"; do [[ -s "$f" ]] || { echo "[ERROR] Missing $f" >&2; exit 1; }; done
mkdir -p "$OUTDIR"/{input,models,evaluation,provenance}
[[ "$ANALYSIS_STATUS" != DEVELOPMENT_ONLY ]] || printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
python3 analysis_v2/tests/test_biomarker_abundance_input.py
python3 analysis_v2/scripts/build_biomarker_abundance_input.py --canonical "$CANONICAL_INPUT" --aliases "$ALIASES" --outdir "$OUTDIR/input"
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_pooled_paired_biomarker_models.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_paired_biomarker_models.R --profile-manifest "$OUTDIR/input/biomarker_profile_manifest.tsv" --abundance-long "$OUTDIR/input/biomarker_abundance_long.tsv" --condition-mode pooled --outdir "$OUTDIR/models"
python3 analysis_v2/tests/test_biomarker_propagation.py
python3 analysis_v2/scripts/evaluate_biomarker_propagation.py --calls "$OUTDIR/models/paired_da_results.tsv" --aliases "$ALIASES" --spike-panel "$SPIKE_PANEL" --outdir "$OUTDIR/evaluation" --q-thresholds 0.05,0.10
sha256sum "$CANONICAL_INPUT" "$ANALYSIS_SIF" "$OUTDIR/models/paired_da_results.tsv" "$OUTDIR/evaluation/biomarker_propagation_metrics.tsv" > "$OUTDIR/provenance/run.sha256"
printf 'analysis\tpooled_paired_artificial_biomarker\nanalysis_status\t%s\nstatus\tPASS\n' "$ANALYSIS_STATUS" > "$OUTDIR/SUCCESS"
echo "[PASS] Sealed pooled paired biomarker analysis: $OUTDIR"
