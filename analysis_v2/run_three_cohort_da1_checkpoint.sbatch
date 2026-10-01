#!/bin/bash
#SBATCH --job-name=da1_checkpoint
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=da1_checkpoint.%j.out
#SBATCH --error=da1_checkpoint.%j.err
set -euo pipefail
: "${PROJECT:?Export project}"
: "${ANALYSIS_SIF:?Export analysis image}"
: "${DA_INVENTORY_ROOT:?Export verified inventory root}"
: "${DA1_CHECKPOINT_ROOT:?Export fresh checkpoint root}"
cd "$PROJECT"
test ! -e "$DA1_CHECKPOINT_ROOT"
(cd "$DA_INVENTORY_ROOT"; sha256sum -c --quiet SHA256SUMS)
for cohort in yachida feng zeller; do
    apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
        python3 -B analysis_v2/scripts/build_da1_checkpoint.py \
        --inventory "$DA_INVENTORY_ROOT/$cohort" --cohort "$cohort" \
        --outdir "$DA1_CHECKPOINT_ROOT/$cohort"
    (cd "$DA1_CHECKPOINT_ROOT/$cohort"; sha256sum -c --quiet SHA256SUMS)
    cat "$DA1_CHECKPOINT_ROOT/$cohort/model_design.tsv"
done
git rev-parse HEAD > "$DA1_CHECKPOINT_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$DA1_CHECKPOINT_ROOT/analysis_image.sha256"
printf 'status\tPASS_DA1_CHECKPOINT\nmodels_run\t0\n' > "$DA1_CHECKPOINT_ROOT/SUCCESS"
(cd "$DA1_CHECKPOINT_ROOT"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Three-cohort DA1 checkpoint; no models: $DA1_CHECKPOINT_ROOT"
