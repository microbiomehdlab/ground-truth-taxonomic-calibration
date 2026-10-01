#!/usr/bin/env bash
# One definitive downstream workflow for Yachida, Feng, and Zeller.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
source "$ROOT/analysis_v2/lib/metaphlan_reference.sh"
: "${COHORT:?Set COHORT to yachida, feng, or zeller}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${RUN_ROOT:?Set RUN_ROOT to a new output directory}"
case "$COHORT" in
  yachida)
    COHORT_ENV="${COHORT_ENV:-${YACHIDA_ENV:-$ROOT/config/yachida.strict-production.env}}"
    EXPECTED_SAMPLES=201
    : "${ASSEMBLY_SENSITIVITY_ROOT:?Yachida requires the sealed assembly-sensitivity experiment}"
    ;;
  feng)
    COHORT_ENV="${COHORT_ENV:-${CRC_ENV:-$ROOT/config/feng.strict-production.env}}"
    EXPECTED_SAMPLES=154
    ;;
  zeller)
    COHORT_ENV="${COHORT_ENV:-${CRC_ENV:-$ROOT/config/zeller.strict-production.env}}"
    EXPECTED_SAMPLES=156
    ;;
  *) echo "[ERROR] Unsupported cohort: $COHORT" >&2; exit 1 ;;
esac
[[ -s "$COHORT_ENV" ]] || { echo "[ERROR] Missing cohort environment: $COHORT_ENV" >&2; exit 1; }
source "$COHORT_ENV"
if [[ "$COHORT" == yachida ]]; then
  STATE_DIR="$YACHIDA_STATE_DIR"
  SENSITIVITY_SUCCESS="${ASSEMBLY_SENSITIVITY_SUCCESS:-$ASSEMBLY_SENSITIVITY_ROOT/experiment_seal/SUCCESS}"
else
  STATE_DIR="$CRC_STATE_DIR"
fi
UNIFIED_SEAL="$STATE_DIR/production_seal_v2"
[[ -s "$UNIFIED_SEAL/SUCCESS" && -s "$UNIFIED_SEAL/production_seal.sha256" ]] || {
  echo "[ERROR] Missing unified upstream seal: $UNIFIED_SEAL" >&2; exit 1;
}
(cd "$UNIFIED_SEAL" && sha256sum -c --quiet production_seal.sha256)
MANIFEST="$UNIFIED_SEAL/production_manifest.tsv"
INDEPENDENT_MANIFEST="$UNIFIED_SEAL/production_manifest.independent.tsv"
SPIKE_ENV="${SPIKE_ENV:-$ROOT/work/yachida_67x3/spikein.env}"
[[ ! -s "$SPIKE_ENV" ]] || source "$SPIKE_ENV"
SPIKE_PANEL="${SPIKE_PANEL:-$ROOT/spikes/spike_panel.tsv}"
ALIASES="${ALIASES:-$ROOT/examples/spike_taxon_aliases.csv}"
(cd "$ROOT" && sha256sum -c --quiet analysis_v2/taxon_identity_freeze.sha256)
python3 analysis_v2/scripts/check_spike_panel_identity.py \
  --runtime "$SPIKE_PANEL" --frozen "$ROOT/spikes/spike_panel.tsv" --root "$ROOT"
cmp -s "$ALIASES" "$ROOT/examples/spike_taxon_aliases.csv" || {
  echo "[ERROR] Taxon aliases differ from frozen identity policy" >&2; exit 1;
}
CANONICAL_INPUT="${CANONICAL_INPUT:-$RUN_ROOT/canonical/canonical_input.tsv}"
CANONICAL_SUCCESS="${CANONICAL_VALIDATION_SUCCESS:-$(dirname "$CANONICAL_INPUT")/validation/SUCCESS}"
[[ ! -e "$RUN_ROOT" || -z "$(find "$RUN_ROOT" -mindepth 1 -print -quit)" ]] || {
  echo "[ERROR] RUN_ROOT must be new or empty" >&2; exit 1;
}
mkdir -p "$RUN_ROOT"/{canonical,readiness,profiler_semantics,models,reports,provenance}
python3 analysis_v2/scripts/validate_analysis_policy.py \
  --policy analysis_v2/ANALYSIS_POLICY.tsv --outdir "$RUN_ROOT/provenance/analysis_policy"
if [[ ! -s "$CANONICAL_INPUT" ]]; then
  [[ "$CANONICAL_INPUT" == "$RUN_ROOT/canonical/canonical_input.tsv" ]] || {
    echo "[ERROR] Supplied canonical input is missing: $CANONICAL_INPUT" >&2; exit 1;
  }
  python3 analysis_v2/scripts/build_crc_cohort_canonical_input.py \
    --cohort "$COHORT" --manifest "$MANIFEST" --independent-manifest "$INDEPENDENT_MANIFEST" \
    --results-root "$PERSISTENT_RESULTS_ROOT" --spike-panel "$SPIKE_PANEL" \
    --aliases "$ALIASES" --expected-samples "$EXPECTED_SAMPLES" --outdir "$RUN_ROOT/canonical"
fi
readiness_args=(--cohort "$COHORT" --manifest "$MANIFEST" --state-dir "$STATE_DIR"
  --canonical "$CANONICAL_INPUT" --canonical-success "$CANONICAL_SUCCESS"
  --analysis-sif "$ANALYSIS_SIF" --expected-samples "$EXPECTED_SAMPLES"
  --expected-independent 30 --outdir "$RUN_ROOT/readiness")
if [[ "$COHORT" == yachida ]]; then
  readiness_args+=(--assembly-sensitivity-success "$SENSITIVITY_SUCCESS")
fi
python3 analysis_v2/scripts/check_cohort_definitive_readiness.py "${readiness_args[@]}"
if require_reference_inputs "$CANONICAL_INPUT"; then
  echo "[PASS] MetaPhlAn reference inputs validated"
fi
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  echo "[PASS] $COHORT definitive preflight only"
  exit 0
fi
audit_args=(--outdir "$RUN_ROOT/profiler_semantics")
while IFS= read -r profile; do
  case "$profile" in
    *.bracken.S.tsv) audit_args+=(--bracken "$profile") ;;
    *.metaphlan.tsv) audit_args+=(--metaphlan "$profile") ;;
    *) echo "[ERROR] Unexpected native profile: $profile" >&2; exit 1 ;;
  esac
done < <(awk -F '\t' 'NR>1 && $22==1 {print $20}' "$CANONICAL_INPUT" | sort -u)
python3 analysis_v2/scripts/audit_profiler_semantics.py "${audit_args[@]}"
derive_endpoints_with_references "$CANONICAL_INPUT" "$RUN_ROOT/endpoints"
if [[ "$COHORT" == yachida ]]; then
  ORIGINAL_ROOT="${YACHIDA_ORIGINAL_ROOT:-$(dirname "$PERSISTENT_RESULTS_ROOT")}" \
  SENSITIVITY_ROOT="$ASSEMBLY_SENSITIVITY_ROOT" \
  BASELINE_ROOT="${YACHIDA_BASELINE_ROOT:-$(dirname "$PERSISTENT_RESULTS_ROOT")}" \
  ANALYSIS_SIF="$ANALYSIS_SIF" OUTDIR="$RUN_ROOT/models/assembly_sensitivity" \
  ANALYSIS_STATUS=DEFINITIVE bash analysis_v2/run_assembly_sensitivity.sh
fi
for population in independent community; do
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_detection_dose_response.R \
    --input "$CANONICAL_INPUT" --outdir "$RUN_ROOT/models/detection_$population" \
    --cohort "$COHORT" --population "$population" --assembly-arm original
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_continuous_dose_response.R \
    --input "$RUN_ROOT/endpoints/paired_endpoints.tsv" --outdir "$RUN_ROOT/models/continuous_$population" \
    --cohort "$COHORT" --population "$population" --assembly-arm original \
    --reference-scale profiler_scale
done
common=(CANONICAL_INPUT="$CANONICAL_INPUT" CANONICAL_VALIDATION_SUCCESS="$CANONICAL_SUCCESS"
  ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEFINITIVE)
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
hash_inputs=("$MANIFEST" "$INDEPENDENT_MANIFEST" "$SPIKE_PANEL" "$ALIASES"
  "$CANONICAL_INPUT" "$ANALYSIS_SIF" "$RUN_ROOT/readiness/SUCCESS"
  "$RUN_ROOT/profiler_semantics/SUCCESS" "$RUN_ROOT/endpoints/SUCCESS"
  "$RUN_ROOT/provenance/analysis_policy/analysis_policy.tsv"
  "$RUN_ROOT/provenance/analysis_policy/analysis_policy.sha256"
  "$RUN_ROOT/provenance/analysis_policy/SUCCESS"
  "$RUN_ROOT/reports/artificial/SUCCESS" "$RUN_ROOT/reports/disease/SUCCESS"
  "$RUN_ROOT/reports/calibration_linkage/SUCCESS")
if [[ "$COHORT" == yachida ]]; then
  hash_inputs+=("$SENSITIVITY_SUCCESS" "$RUN_ROOT/models/assembly_sensitivity/SUCCESS")
fi
sha256sum "${hash_inputs[@]}" > "$RUN_ROOT/provenance/definitive_run.sha256"
printf 'analysis\t%s_definitive_v2\nstatus\tPASS\n' "$COHORT" > "$RUN_ROOT/SUCCESS"
echo "[PASS] Definitive $COHORT analysis sealed: $RUN_ROOT"
