#!/usr/bin/env bash
# Slide 3-5 previews from already completed native profiles and corrected endpoints.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${FIGURE_ROOT:?Set FIGURE_ROOT to the validated three-cohort recoverability run}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the validated analysis image}"
: "${OUTDIR:?Set OUTDIR to a new presentation-only output directory}"
NATIVE="$FIGURE_ROOT/partial_snapshot_response/native_abundance"
ENDPOINTS="$FIGURE_ROOT/endpoints/paired_endpoints.tsv"
for path in "$NATIVE/SUCCESS" "$NATIVE/biomarker_profile_manifest.tsv" \
            "$NATIVE/biomarker_abundance_long.tsv" \
            "$FIGURE_ROOT/endpoints/SUCCESS" "$ENDPOINTS" "$ANALYSIS_SIF" \
            examples/spike_taxon_aliases.csv spikes/spike_panel.tsv; do
  test -s "$path" || { echo "[ERROR] Missing prerequisite: $path" >&2; exit 1; }
done
test ! -e "$OUTDIR" || { echo "[ERROR] OUTDIR exists: $OUTDIR" >&2; exit 1; }
mkdir -p "$OUTDIR"
printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
python3 analysis_v2/scripts/build_three_cohort_baseline_figure_input.py \
  --manifest "$NATIVE/biomarker_profile_manifest.tsv" \
  --abundance "$NATIVE/biomarker_abundance_long.tsv" \
  --aliases examples/spike_taxon_aliases.csv \
  --panel spikes/spike_panel.tsv \
  --outdir "$OUTDIR/baseline_input"
test -s "$OUTDIR/baseline_input/SUCCESS"
apptainer exec --cleanenv --bind "$ROOT:$ROOT" --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/plot_institute_presentation_panels.R \
  --baseline "$OUTDIR/baseline_input/baseline_four_taxa_samples.tsv" \
  --endpoints "$ENDPOINTS" \
  --outdir "$OUTDIR/slides"
for stem in slide3_baseline_visibility slide4_focused_recovery slide5_fnuc_detection; do
  test -s "$OUTDIR/slides/$stem.png" || { echo "[ERROR] Missing $stem.png" >&2; exit 1; }
  test -s "$OUTDIR/slides/$stem.pdf" || { echo "[ERROR] Missing $stem.pdf" >&2; exit 1; }
done
test -s "$OUTDIR/slides/SUCCESS"
printf 'status\tPASS\nanalysis_status\tDEVELOPMENT_ONLY\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Institute presentation panels: $OUTDIR/slides"
