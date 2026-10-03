#!/usr/bin/env bash
#SBATCH --job-name=maaslin_timing
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=maaslin_timing.%A_%a.out
#SBATCH --error=maaslin_timing.%A_%a.err
set -euo pipefail
: "${MAASLIN_TIMING_ROOT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}"
cd "$MAASLIN_TIMING_ROOT/source"
mode="${1:?prepare, run or collect}"
if [[ "$mode" == prepare ]]; then
    sha256sum "$ANALYSIS_SIF" > "$MAASLIN_TIMING_ROOT/image.sha256"
    sha256sum analysis_v2/scripts/maaslin_timing.py \
        analysis_v2/scripts/time_maaslin_context.R analysis_v2/lib/maaslin_context.R \
        analysis_v2/lib/maaslin_contract.R analysis_v2/lib/unpaired_null.R \
        analysis_v2/scripts/audit_bracken_denominators.py \
        analysis_v2/scripts/build_biomarker_abundance_input.py \
        analysis_v2/scripts/prepare_da_pilot.py analysis_v2/scripts/plan_da3_canary.py \
        > "$MAASLIN_TIMING_ROOT/source_files.sha256"
fi
sha256sum -c --quiet "$MAASLIN_TIMING_ROOT/image.sha256"
sha256sum -c --quiet "$MAASLIN_TIMING_ROOT/source_files.sha256"
runtime=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PWD"
    --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,REQUIRE_MAASLIN_TIMING_BACKEND=1
    "$ANALYSIS_SIF")
case "$mode" in
    prepare)
        "${runtime[@]}" Rscript analysis_v2/tests/test_maaslin_timing.R
        "${runtime[@]}" python3 -B analysis_v2/scripts/maaslin_timing.py prepare \
            --source "$DA3_BATCH_PLAN" --out "$MAASLIN_TIMING_ROOT/plan"
        ;;
    run)
        "${runtime[@]}" python3 -B analysis_v2/scripts/maaslin_timing.py run \
            --plan "$MAASLIN_TIMING_ROOT/plan" --results "$MAASLIN_TIMING_ROOT/results" \
            --repo "$PWD" --index "${SLURM_ARRAY_TASK_ID:?}"
        ;;
    collect)
        "${runtime[@]}" python3 -B analysis_v2/scripts/maaslin_timing.py collect \
            --plan "$MAASLIN_TIMING_ROOT/plan" --results "$MAASLIN_TIMING_ROOT/results" \
            --out "$MAASLIN_TIMING_ROOT/REPORT"
        ;;
    *) echo "Unknown mode: $mode" >&2; exit 1 ;;
esac
