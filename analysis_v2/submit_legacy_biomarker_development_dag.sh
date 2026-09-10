#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"

: "${DEV_INPUT_ROOT:?Export DEV_INPUT_ROOT before submitting}"
: "${ANALYSIS_SIF:?Export ANALYSIS_SIF before submitting}"
: "${OUTDIR:?Export OUTDIR to a new output directory}"
SBATCH_BIN="${SBATCH_BIN:-sbatch}"

[[ ! -e "$OUTDIR" ]] || {
  echo "[ERROR] OUTDIR already exists; choose a new directory: $OUTDIR" >&2
  exit 1
}
mkdir -p "$OUTDIR/slurm"
export PROJECT="$ROOT" DEV_INPUT_ROOT ANALYSIS_SIF OUTDIR

submit() {
  local stage="$1" dependency="$2" cpus="$3" memory="$4" time_limit="$5"
  local args=(--parsable --job-name="legacy_${stage}" --cpus-per-task="$cpus"
    --mem="$memory" --time="$time_limit"
    --output="$OUTDIR/slurm/${stage}.%j.out"
    --error="$OUTDIR/slurm/${stage}.%j.err"
    --export="ALL,PIPELINE_STAGE=$stage")
  [[ -z "$dependency" ]] || args+=(--dependency="afterok:$dependency")
  "$SBATCH_BIN" "${args[@]}" analysis_v2/run_legacy_biomarker_development_stage.sbatch
}

setup_job="$(submit setup '' 1 8G 02:00:00)"
pooled_job="$(submit artificial_pooled "$setup_job" 4 64G 2-00:00:00)"
stratified_job="$(submit artificial_stratified "$setup_job" 4 64G 2-00:00:00)"
disease_job="$(submit disease "$setup_job" 4 64G 2-00:00:00)"
artificial_report_job="$(submit artificial_report "$pooled_job:$stratified_job" 2 24G 12:00:00)"
disease_report_job="$(submit disease_report "$disease_job" 2 24G 12:00:00)"
linkage_job="$(submit calibration_linkage "$pooled_job:$stratified_job" 2 24G 12:00:00)"
seal_job="$(submit seal "$artificial_report_job:$disease_report_job:$linkage_job" 1 4G 01:00:00)"

cat <<EOF
[SUBMITTED] Parallel legacy biomarker-development DAG
OUTDIR=$OUTDIR
setup=$setup_job
artificial_pooled=$pooled_job
artificial_stratified=$stratified_job
disease=$disease_job
artificial_report=$artificial_report_job
disease_report=$disease_report_job
calibration_linkage=$linkage_job
seal=$seal_job

Monitor:
  squeue -j $setup_job,$pooled_job,$stratified_job,$disease_job,$artificial_report_job,$disease_report_job,$linkage_job,$seal_job
Final success:
  test -s '$OUTDIR/SUCCESS' && echo '[PASS] Parallel run complete'
EOF
