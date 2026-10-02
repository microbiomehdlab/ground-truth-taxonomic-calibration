#!/bin/bash
#SBATCH --job-name=da2_signflip
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=da2_signflip.%A_%a.out
#SBATCH --error=da2_signflip.%A_%a.err
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA2_PAIRED_RESULTS:?}" "${DA2_NULL_ROOT:?}" "${SLURM_ARRAY_TASK_ID:?}"
cd "$PROJECT"
printf -v context 'pilot_%04d' "$SLURM_ARRAY_TASK_ID"
source_dir="$DA2_PAIRED_RESULTS/$context"
out="$DA2_NULL_ROOT/$context"
test -s "$source_dir/SUCCESS"
test ! -e "$out"
(cd "$source_dir"; sha256sum -c --quiet SHA256SUMS)
mkdir -p "$DA2_NULL_ROOT"
attempt="$(mktemp -d "$DA2_NULL_ROOT/${context}_attempt_XXXXXX")"
sha256sum "$source_dir/SHA256SUMS" > "$attempt/input.sha256"
sha256sum analysis_v2/lib/paired_null.R analysis_v2/scripts/run_da2_signflip.R \
  analysis_v2/run_da2_signflip.sbatch > "$attempt/code.sha256"
sha256sum "$ANALYSIS_SIF" > "$attempt/image.sha256"
apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1 \
  --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/scripts/run_da2_signflip.R "$source_dir" "$attempt/results" "$PROJECT"
(cd "$source_dir"; sha256sum -c --quiet SHA256SUMS)
sha256sum -c --quiet "$attempt/input.sha256" "$attempt/code.sha256" "$attempt/image.sha256"
git rev-parse HEAD > "$attempt/source_commit.txt"
printf 'status\tPASS_CONDITIONAL_SIGNFLIP_DIAGNOSTIC\ndefinitive\t0\n' > "$attempt/SUCCESS"
(cd "$attempt"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
mv -T -n "$attempt" "$out"
test ! -e "$attempt"
