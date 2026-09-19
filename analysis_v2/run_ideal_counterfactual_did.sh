#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"

: "${DISEASE_RUN:?Set DISEASE_RUN to a sealed disease-biomarker analysis}"
: "${PAIRED_ENDPOINTS:?Set PAIRED_ENDPOINTS to paired_endpoints.tsv}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty ideal-counterfactual directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
# No implicit reference-scale default: the caller states the scale.
: "${REFERENCE_SCALE:?Set REFERENCE_SCALE to profiler_scale or read_proportional}"
case "$REFERENCE_SCALE" in
  profiler_scale|read_proportional) ;;
  *) echo "[ERROR] REFERENCE_SCALE must be profiler_scale or read_proportional, got: $REFERENCE_SCALE" >&2; exit 1 ;;
esac
reference_args=(--reference-scale "$REFERENCE_SCALE")
if [[ "$REFERENCE_SCALE" == "profiler_scale" ]]; then
  # The profiler-scale expected abundance is taken from the matching all-feature
  # response table; it is never reconstructed from F and the focal fraction.
  : "${RESPONSE_TABLE:?REFERENCE_SCALE=profiler_scale requires RESPONSE_TABLE (paired_feature_responses.parquet from the same cohort, population, assembly arm, canonical input and endpoint generation)}"
  if [[ ! -s "$RESPONSE_TABLE" ]]; then
    echo "[ERROR] RESPONSE_TABLE missing or empty: $RESPONSE_TABLE" >&2
    exit 1
  fi
  reference_args+=(--response-table "$RESPONSE_TABLE")
fi

PROFILE_MANIFEST="${PROFILE_MANIFEST:-$DISEASE_RUN/input/biomarker_profile_manifest.tsv}"
ABUNDANCE_LONG="${ABUNDANCE_LONG:-$DISEASE_RUN/input/biomarker_abundance_long.tsv}"
PAIRED_ENDPOINTS_SUCCESS="${PAIRED_ENDPOINTS_SUCCESS:-$(dirname "$PAIRED_ENDPOINTS")/SUCCESS}"

[[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" || "$ANALYSIS_STATUS" == "DEFINITIVE" ]] || {
  echo "[ERROR] ANALYSIS_STATUS must be DEVELOPMENT_ONLY or DEFINITIVE" >&2
  exit 1
}
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "[ERROR] OUTDIR must be new or empty: $OUTDIR" >&2
  exit 1
fi
for path in "$DISEASE_RUN/SUCCESS" "$PROFILE_MANIFEST" "$ABUNDANCE_LONG" \
            "$PAIRED_ENDPOINTS" "$PAIRED_ENDPOINTS_SUCCESS" "$ANALYSIS_SIF"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing required input: $path" >&2; exit 1; }
done
if [[ "$ANALYSIS_STATUS" == "DEFINITIVE" ]]; then
  for marker in \
    "$DISEASE_RUN/DEVELOPMENT_ONLY.txt" \
    "$(dirname "$PAIRED_ENDPOINTS")/DEVELOPMENT_ONLY.txt" \
    "$(dirname "$(dirname "$PAIRED_ENDPOINTS")")/DEVELOPMENT_ONLY.txt"; do
    [[ ! -e "$marker" ]] || {
      echo "[ERROR] Development input cannot produce a definitive analysis: $marker" >&2
      exit 1
    }
  done
fi

mkdir -p "$OUTDIR"/{model,provenance}
{
  printf 'field\tvalue\n'
  printf 'status\tIN_PROGRESS\n'
  printf 'analysis_status\t%s\n' "$ANALYSIS_STATUS"
  printf 'created_at\t%s\n' "$(date -Iseconds)"
  printf 'repository_commit\t%s\n' "$(git rev-parse HEAD)"
  printf 'disease_run\t%s\n' "$(realpath "$DISEASE_RUN")"
  printf 'profile_manifest\t%s\n' "$(realpath "$PROFILE_MANIFEST")"
  printf 'abundance_long\t%s\n' "$(realpath "$ABUNDANCE_LONG")"
  printf 'paired_endpoints\t%s\n' "$(realpath "$PAIRED_ENDPOINTS")"
  printf 'analysis_sif\t%s\n' "$(realpath "$ANALYSIS_SIF")"
} > "$OUTDIR/provenance/run_manifest.tsv"
if [[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" ]]; then
  printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
fi

apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_ideal_counterfactual_did.R
complete_args=()
[[ "$ANALYSIS_STATUS" != "DEFINITIVE" ]] || complete_args=(--require-complete-panel)
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/fit_ideal_counterfactual_did.R \
  --profile-manifest "$PROFILE_MANIFEST" \
  --abundance-long "$ABUNDANCE_LONG" \
  --paired-endpoints "$PAIRED_ENDPOINTS" \
  --outdir "$OUTDIR/model" "${reference_args[@]}" "${complete_args[@]}"

for path in \
  "$OUTDIR/model/SUCCESS" \
  "$OUTDIR/model/ideal_counterfactual_did_results.tsv" \
  "$OUTDIR/model/ideal_counterfactual_panel_audit.tsv" \
  "$OUTDIR/model/ideal_counterfactual_sample_panel_audit.tsv" \
  "$OUTDIR/model/ideal_counterfactual.sha256"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing model output: $path" >&2; exit 1; }
done

sha256sum \
  "$DISEASE_RUN/SUCCESS" "$PROFILE_MANIFEST" "$ABUNDANCE_LONG" \
  "$PAIRED_ENDPOINTS" "$PAIRED_ENDPOINTS_SUCCESS" \
  "$OUTDIR/model/ideal_counterfactual_did_results.tsv" \
  "$OUTDIR/model/ideal_counterfactual_panel_audit.tsv" \
  > "$OUTDIR/provenance/run_inputs_and_primary_outputs.sha256"
sed -i 's/^status\tIN_PROGRESS$/status\tPASS/' "$OUTDIR/provenance/run_manifest.tsv"
printf 'analysis\tideal_counterfactual_difference_in_differences\nanalysis_status\t%s\nstatus\tPASS\n' \
  "$ANALYSIS_STATUS" > "$OUTDIR/SUCCESS"
echo "[PASS] Sealed ideal-counterfactual difference-in-differences: $OUTDIR"
