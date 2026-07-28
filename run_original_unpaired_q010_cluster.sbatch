#!/usr/bin/env bash
#SBATCH --job-name=orig_unpaired_q010
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=orig_unpaired_q010_%j.out
#SBATCH --error=orig_unpaired_q010_%j.err

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
: "${SIF:?Set SIF to the MaAsLin2 analysis container image}"
RUN_ROOT="${RUN_ROOT:-RUNS_publication_original_unpaired_q010_${SLURM_JOB_ID:-manual}}"

cd "$PROJECT"
export PROJECT SIF

echo "[INFO] Job ID:   ${SLURM_JOB_ID:-not-running-under-slurm}"
echo "[INFO] Host:     $(hostname)"
echo "[INFO] Project:  $PROJECT"
echo "[INFO] SIF:      $SIF"
echo "[INFO] Run root: $RUN_ROOT"
echo "[INFO] Started:  $(date --iso-8601=seconds)"

R_BIN="$PROJECT/rscript_in_container.sh" \
bash cluster_preflight_original_unpaired.sh

PROJECT="$PROJECT" \
SIF="$SIF" \
RUN_ROOT="$RUN_ROOT" \
R_BIN="$PROJECT/rscript_in_container.sh" \
ALLOW_EXISTING="${ALLOW_EXISTING:-false}" \
bash run_publication_original_unpaired_q010.sh
