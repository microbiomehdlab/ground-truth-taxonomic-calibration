#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${DEV_INPUT_ROOT:?Set DEV_INPUT_ROOT to the sealed legacy native development input}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty output directory}"
CANONICAL_INPUT="$DEV_INPUT_ROOT/canonical_input.tsv"; CANONICAL_VALIDATION_SUCCESS="$DEV_INPUT_ROOT/validation/SUCCESS"
ENDPOINTS_FILE="$DEV_INPUT_ROOT/endpoints/paired_endpoints.tsv"; ENDPOINTS_SUCCESS="$DEV_INPUT_ROOT/endpoints/SUCCESS"
FENG_MANIFEST="${FENG_MANIFEST:-datasets/fengq/manifests/production_manifest.tsv}"; ZELLER_MANIFEST="${ZELLER_MANIFEST:-datasets/zellerg/manifests/production_manifest.tsv}"
for path in "$CANONICAL_INPUT" "$CANONICAL_VALIDATION_SUCCESS" "$ENDPOINTS_FILE" "$ENDPOINTS_SUCCESS" "$DEV_INPUT_ROOT/DEVELOPMENT_ONLY.txt" "$FENG_MANIFEST" "$ZELLER_MANIFEST" "$ANALYSIS_SIF"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing development input: $path" >&2; exit 1; }
done
[[ ! -e "$OUTDIR" || -z "$(find "$OUTDIR" -mindepth 1 -print -quit)" ]] || { echo "[ERROR] OUTDIR must be new or empty" >&2; exit 1; }
mkdir -p "$OUTDIR"/{metadata,models,reports,provenance}; printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
python3 analysis_v2/tests/test_development_sample_metadata.py
python3 analysis_v2/scripts/build_development_sample_metadata.py --canonical "$CANONICAL_INPUT" --feng-manifest "$FENG_MANIFEST" --zeller-manifest "$ZELLER_MANIFEST" --output "$OUTDIR/metadata/sample_metadata.tsv"
common=(CANONICAL_INPUT="$CANONICAL_INPUT" CANONICAL_VALIDATION_SUCCESS="$CANONICAL_VALIDATION_SUCCESS" ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEVELOPMENT_ONLY)
env "${common[@]}" OUTDIR="$OUTDIR/models/artificial_pooled" bash analysis_v2/run_pooled_paired_biomarker_propagation.sh
env "${common[@]}" SAMPLE_METADATA="$OUTDIR/metadata/sample_metadata.tsv" OUTDIR="$OUTDIR/models/artificial_stratified" bash analysis_v2/run_paired_biomarker_propagation.sh
env "${common[@]}" SAMPLE_METADATA="$OUTDIR/metadata/sample_metadata.tsv" OUTDIR="$OUTDIR/models/disease" bash analysis_v2/run_disease_biomarker_propagation.sh
env PAIRED_RUN="$OUTDIR/models/artificial_pooled" SECONDARY_PAIRED_RUN="$OUTDIR/models/artificial_stratified" ANALYSIS_SIF="$ANALYSIS_SIF" REPORT_STATUS=DEVELOPMENT_ONLY OUTDIR="$OUTDIR/reports/artificial" bash analysis_v2/run_artificial_biomarker_report.sh
env DISEASE_RUN="$OUTDIR/models/disease" ANALYSIS_SIF="$ANALYSIS_SIF" REPORT_STATUS=DEVELOPMENT_ONLY OUTDIR="$OUTDIR/reports/disease" bash analysis_v2/run_disease_biomarker_report.sh
env ENDPOINTS_FILE="$ENDPOINTS_FILE" ENDPOINTS_SUCCESS="$ENDPOINTS_SUCCESS" PAIRED_RUN="$OUTDIR/models/artificial_pooled" SECONDARY_PAIRED_RUN="$OUTDIR/models/artificial_stratified" ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEVELOPMENT_ONLY OUTDIR="$OUTDIR/reports/calibration_linkage" bash analysis_v2/run_calibration_biomarker_linkage.sh
sha256sum "$CANONICAL_INPUT" "$ENDPOINTS_FILE" "$OUTDIR/metadata/sample_metadata.tsv" "$OUTDIR/models/artificial_pooled/SUCCESS" "$OUTDIR/models/artificial_stratified/SUCCESS" "$OUTDIR/models/disease/SUCCESS" "$OUTDIR/reports/artificial/SUCCESS" "$OUTDIR/reports/disease/SUCCESS" "$OUTDIR/reports/calibration_linkage/SUCCESS" > "$OUTDIR/provenance/development_run.sha256"
printf 'analysis\tlegacy_feng_zeller_biomarker_development\nstatus\tDEVELOPMENT_ONLY\ncompletion\tPASS\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Legacy Feng/Zeller biomarker development completed: $OUTDIR"
