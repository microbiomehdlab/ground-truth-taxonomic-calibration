#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${COHORT_RESULT_FILES:?Set colon-separated Yachida:Feng:Zeller primary result files}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the frozen downstream image}"
: "${OUTDIR:?Set OUTDIR to a new, empty synthesis directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
EXPECTED_COHORTS="${EXPECTED_COHORTS:-yachida,feng,zeller}"
[[ "$ANALYSIS_STATUS" == "DEVELOPMENT_ONLY" || "$ANALYSIS_STATUS" == "DEFINITIVE" ]] || { echo "[ERROR] Invalid ANALYSIS_STATUS" >&2; exit 1; }
if [[ -e "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then echo "[ERROR] OUTDIR must be new or empty" >&2; exit 1; fi
IFS=: read -r -a inputs <<< "$COHORT_RESULT_FILES"
[[ "${#inputs[@]}" -eq 3 ]] || { echo "[ERROR] Exactly three cohort result files are required" >&2; exit 1; }
for path in "$ANALYSIS_SIF" "${inputs[@]}"; do [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }; done
mkdir -p "$OUTDIR/provenance"
printf 'status\t%s\nrepository_commit\t%s\nexpected_cohorts\t%s\n' "$ANALYSIS_STATUS" "$(git rev-parse HEAD)" "$EXPECTED_COHORTS" > "$OUTDIR/provenance/run_manifest.tsv"
[[ "$ANALYSIS_STATUS" != "DEVELOPMENT_ONLY" ]] || printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_cross_cohort_meta_analysis.R
comma_inputs="$(IFS=,; echo "${inputs[*]}")"
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/meta_analyze_disease_biomarkers.R \
  --inputs "$comma_inputs" --expected-cohorts "$EXPECTED_COHORTS" --outdir "$OUTDIR/results"
sha256sum "$ANALYSIS_SIF" "${inputs[@]}" "$OUTDIR/results/random_effects_meta_analysis.tsv" > "$OUTDIR/provenance/run_inputs_and_primary_outputs.sha256"
printf 'analysis\tcross_cohort_disease_biomarker_synthesis\nanalysis_status\t%s\nstatus\tPASS\n' "$ANALYSIS_STATUS" > "$OUTDIR/SUCCESS"
echo "[PASS] Sealed cross-cohort synthesis: $OUTDIR"
