#!/usr/bin/env bash
# Definitive downstream orchestrator. It never creates or repairs upstream evidence.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${YACHIDA_ENV:?Set YACHIDA_ENV to the frozen strict-production environment}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${ASSEMBLY_SENSITIVITY_ROOT:?Set ASSEMBLY_SENSITIVITY_ROOT to the sealed 30-sample experiment}"
: "${RUN_ROOT:?Set RUN_ROOT to a new definitive output directory}"
source "$YACHIDA_ENV"
MANIFEST="${YACHIDA_MANIFEST:-$ROOT/work/yachida_67x3/metadata/pilot_batched.tsv}"
INDEPENDENT_MANIFEST="${INDEPENDENT_MANIFEST:-$ROOT/work/yachida_67x3/metadata/independent_10_per_condition.tsv}"
SPIKE_ENV="${SPIKE_ENV:-$ROOT/work/yachida_67x3/spikein.env}"
[[ -s "$SPIKE_ENV" ]] && source "$SPIKE_ENV"
SPIKE_PANEL="${SPIKE_PANEL:-$ROOT/spikes/spike_panel.tsv}"
ALIASES="${ALIASES:-$ROOT/examples/spike_taxon_aliases.csv}"
CANONICAL_INPUT="${CANONICAL_INPUT:-$RUN_ROOT/canonical/canonical_input.tsv}"
CANONICAL_SUCCESS="${CANONICAL_VALIDATION_SUCCESS:-$(dirname "$CANONICAL_INPUT")/validation/SUCCESS}"
SENSITIVITY_SUCCESS="${ASSEMBLY_SENSITIVITY_SUCCESS:-$ASSEMBLY_SENSITIVITY_ROOT/experiment_seal/SUCCESS}"
[[ ! -e "$RUN_ROOT" || -z "$(find "$RUN_ROOT" -mindepth 1 -print -quit)" ]] || { echo "[ERROR] RUN_ROOT must be new or empty" >&2; exit 1; }
mkdir -p "$RUN_ROOT"/{canonical,readiness,profiler_semantics,endpoints,models,reports,provenance}

if [[ ! -s "$CANONICAL_INPUT" ]]; then
  [[ "$CANONICAL_INPUT" == "$RUN_ROOT/canonical/canonical_input.tsv" ]] || {
    echo "[ERROR] Externally supplied CANONICAL_INPUT is missing: $CANONICAL_INPUT" >&2
    exit 1
  }
  python3 analysis_v2/scripts/build_yachida_canonical_input.py \
    --manifest "$MANIFEST" --independent-manifest "$INDEPENDENT_MANIFEST" \
    --results-root "$PERSISTENT_RESULTS_ROOT" --spike-panel "$SPIKE_PANEL" \
    --aliases "$ALIASES" --outdir "$RUN_ROOT/canonical"
fi

python3 analysis_v2/scripts/check_yachida_definitive_readiness.py \
  --manifest "$MANIFEST" --state-dir "$YACHIDA_STATE_DIR" \
  --canonical "$CANONICAL_INPUT" --canonical-success "$CANONICAL_SUCCESS" \
  --analysis-sif "$ANALYSIS_SIF" --assembly-sensitivity-success "$SENSITIVITY_SUCCESS" \
  --outdir "$RUN_ROOT/readiness"
[[ "${PREFLIGHT_ONLY:-0}" != 1 ]] || { echo "[PASS] Definitive preflight only"; exit 0; }

audit_args=(--outdir "$RUN_ROOT/profiler_semantics")
while IFS= read -r profile; do
  case "$profile" in
    *.bracken.S.tsv) audit_args+=(--bracken "$profile") ;;
    *.metaphlan.tsv) audit_args+=(--metaphlan "$profile") ;;
    *) echo "[ERROR] Unexpected native profile: $profile" >&2; exit 1 ;;
  esac
done < <(awk -F '\t' 'NR>1 && $22==1 {print $20}' "$CANONICAL_INPUT" | sort -u)
python3 analysis_v2/scripts/audit_profiler_semantics.py "${audit_args[@]}"
python3 analysis_v2/scripts/derive_paired_endpoints.py --input "$CANONICAL_INPUT" --outdir "$RUN_ROOT/endpoints"

ORIGINAL_ROOT="${YACHIDA_ORIGINAL_ROOT:-$(dirname "$PERSISTENT_RESULTS_ROOT")}" \
SENSITIVITY_ROOT="$ASSEMBLY_SENSITIVITY_ROOT" BASELINE_ROOT="${YACHIDA_BASELINE_ROOT:-$(dirname "$PERSISTENT_RESULTS_ROOT")}" \
ANALYSIS_SIF="$ANALYSIS_SIF" OUTDIR="$RUN_ROOT/models/assembly_sensitivity" \
ANALYSIS_STATUS=DEFINITIVE \
  bash analysis_v2/run_assembly_sensitivity.sh

for population in independent community; do
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_detection_dose_response.R \
    --input "$RUN_ROOT/endpoints/paired_endpoints.tsv" --outdir "$RUN_ROOT/models/detection_$population" \
    --cohort yachida --population "$population" --assembly-arm original
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_continuous_dose_response.R \
    --input "$RUN_ROOT/endpoints/paired_endpoints.tsv" --outdir "$RUN_ROOT/models/continuous_$population" \
    --cohort yachida --population "$population" --assembly-arm original
done

common=(CANONICAL_INPUT="$CANONICAL_INPUT" CANONICAL_VALIDATION_SUCCESS="$CANONICAL_SUCCESS" ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEFINITIVE)
env "${common[@]}" OUTDIR="$RUN_ROOT/models/artificial_pooled" bash analysis_v2/run_pooled_paired_biomarker_propagation.sh
env "${common[@]}" SAMPLE_METADATA="$MANIFEST" OUTDIR="$RUN_ROOT/models/artificial_stratified" bash analysis_v2/run_paired_biomarker_propagation.sh
env "${common[@]}" SAMPLE_METADATA="$MANIFEST" OUTDIR="$RUN_ROOT/models/disease" bash analysis_v2/run_disease_biomarker_propagation.sh

env PAIRED_RUN="$RUN_ROOT/models/artificial_pooled" SECONDARY_PAIRED_RUN="$RUN_ROOT/models/artificial_stratified" \
  ANALYSIS_SIF="$ANALYSIS_SIF" REPORT_STATUS=DEFINITIVE OUTDIR="$RUN_ROOT/reports/artificial" \
  bash analysis_v2/run_artificial_biomarker_report.sh
env DISEASE_RUN="$RUN_ROOT/models/disease" ANALYSIS_SIF="$ANALYSIS_SIF" REPORT_STATUS=DEFINITIVE \
  OUTDIR="$RUN_ROOT/reports/disease" bash analysis_v2/run_disease_biomarker_report.sh
env ENDPOINTS_FILE="$RUN_ROOT/endpoints/paired_endpoints.tsv" ENDPOINTS_SUCCESS="$RUN_ROOT/endpoints/SUCCESS" \
  PAIRED_RUN="$RUN_ROOT/models/artificial_pooled" SECONDARY_PAIRED_RUN="$RUN_ROOT/models/artificial_stratified" \
  ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEFINITIVE OUTDIR="$RUN_ROOT/reports/calibration_linkage" \
  bash analysis_v2/run_calibration_biomarker_linkage.sh

sha256sum "$MANIFEST" "$INDEPENDENT_MANIFEST" "$SPIKE_PANEL" "$ALIASES" \
  "$CANONICAL_INPUT" "$ANALYSIS_SIF" "$SENSITIVITY_SUCCESS" \
  "$RUN_ROOT/readiness/SUCCESS" "$RUN_ROOT/profiler_semantics/SUCCESS" "$RUN_ROOT/endpoints/SUCCESS" \
  "$RUN_ROOT/models/assembly_sensitivity/SUCCESS" \
  "$RUN_ROOT/reports/artificial/SUCCESS" "$RUN_ROOT/reports/disease/SUCCESS" \
  "$RUN_ROOT/reports/calibration_linkage/SUCCESS" > "$RUN_ROOT/provenance/definitive_run.sha256"
printf 'analysis\tyachida_definitive_v2\nstatus\tPASS\n' > "$RUN_ROOT/SUCCESS"
echo "[PASS] Definitive Yachida analysis sealed: $RUN_ROOT"
