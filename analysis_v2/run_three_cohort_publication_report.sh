#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${COHORT_RUNS:?Set colon-separated cohort analysis packages}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF}"
: "${OUTDIR:?Set OUTDIR to a new report directory}"
: "${REPORT_STATUS:?Set REPORT_STATUS to DEVELOPMENT_ONLY or DEFINITIVE}"
EXPECTED_COHORTS="${EXPECTED_COHORTS:-yachida,feng,zeller}"
[[ "$REPORT_STATUS" == DEVELOPMENT_ONLY || "$REPORT_STATUS" == DEFINITIVE ]] || { echo "[ERROR] Invalid REPORT_STATUS" >&2; exit 1; }
[[ -s "$ANALYSIS_SIF" ]] || { echo "[ERROR] Missing analysis image: $ANALYSIS_SIF" >&2; exit 1; }
IFS=: read -r -a runs <<< "$COHORT_RUNS"
[[ ${#runs[@]} -gt 0 ]] || { echo "[ERROR] No cohort runs" >&2; exit 1; }
artificial=(); disease=(); linkage=()
for run in "${runs[@]}"; do
  for path in "$run/SUCCESS" "$run/reports/artificial/SUCCESS" "$run/reports/disease/SUCCESS" "$run/reports/calibration_linkage/SUCCESS"; do
    [[ -s "$path" ]] || { echo "[ERROR] Incomplete cohort report package: $path" >&2; exit 1; }
  done
  if [[ "$REPORT_STATUS" == DEFINITIVE ]]; then
    [[ ! -e "$run/DEVELOPMENT_ONLY.txt" && -s "$run/readiness/SUCCESS" && -s "$run/provenance/definitive_run.sha256" ]] || {
      echo "[ERROR] Non-definitive cohort package: $run" >&2; exit 1;
    }
  fi
  artificial+=("$run/reports/artificial"); disease+=("$run/reports/disease"); linkage+=("$run/reports/calibration_linkage")
done
[[ "$REPORT_STATUS" != DEFINITIVE || ${#runs[@]} -eq 3 ]] || { echo "[ERROR] Definitive reporting requires three cohort packages" >&2; exit 1; }
[[ ! -e "$OUTDIR" || -z "$(find "$OUTDIR" -mindepth 1 -print -quit)" ]] || { echo "[ERROR] OUTDIR must be new or empty" >&2; exit 1; }
join_by_comma() { local IFS=,; echo "$*"; }
command=(apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/make_three_cohort_publication_report.R
  --artificial-reports "$(join_by_comma "${artificial[@]}")" --disease-reports "$(join_by_comma "${disease[@]}")"
  --linkage-reports "$(join_by_comma "${linkage[@]}")" --expected-cohorts "$EXPECTED_COHORTS"
  --report-status "$REPORT_STATUS" --outdir "$OUTDIR")
[[ -z "${SYNTHESIS_RUN:-}" ]] || command+=(--synthesis "$SYNTHESIS_RUN")
[[ -z "${ASSEMBLY_REPORT:-}" ]] || command+=(--assembly-report "$ASSEMBLY_REPORT")
if [[ "$REPORT_STATUS" == DEFINITIVE ]]; then
  [[ -n "${SYNTHESIS_RUN:-}" && -s "$SYNTHESIS_RUN/SUCCESS" && -n "${ASSEMBLY_REPORT:-}" && -s "$ASSEMBLY_REPORT/SUCCESS" ]] || {
    echo "[ERROR] Definitive reporting requires sealed synthesis and assembly-sensitivity report" >&2; exit 1;
  }
fi
"${command[@]}"
