#!/usr/bin/env bash
# Generate and profile one independent and one community spike without deletion.
set -euo pipefail
IFS=$'\n\t'
umask 002

: "${PROJECT:?set PROJECT to the repository root}"
: "${YACHIDA_ENV:?set YACHIDA_ENV to config/yachida.env}"
: "${SAMPLE_ID:?set the retained preprocessing-smoke sample ID}"

source "$YACHIDA_ENV"
SPIKE_ENV="${SPIKE_ENV:-$PROJECT/work/yachida_67x3/spikein.env}"
SAMPLE_WORK="${SAMPLE_WORK:-$PROJECT/work/yachida_67x3/smoke_scratch/$SAMPLE_ID}"
SMOKE_ROOT="${YACHIDA_SPIKE_SMOKE_ROOT:-$PROJECT/work/yachida_67x3/spike_smoke/$SAMPLE_ID}"
[[ -s "$SPIKE_ENV" ]] || { echo "[ERROR] Missing spike environment: $SPIKE_ENV" >&2; exit 1; }
requested_sampling_mode="${SAMPLING_MODE:-}"
source "$SPIKE_ENV"
if [[ -n "$requested_sampling_mode" ]]; then
  SAMPLING_MODE="$requested_sampling_mode"
fi
SAMPLING_MODE="${SAMPLING_MODE:-single_pass}"
[[ "$SAMPLING_MODE" == legacy || "$SAMPLING_MODE" == single_pass ]] || {
  echo "[ERROR] SAMPLING_MODE must be legacy or single_pass" >&2
  exit 1
}

: "${IMG:?set IMG in the spike environment}"
: "${SPIKE_PANEL:?set SPIKE_PANEL in the spike environment}"
: "${POOLS_DIR:?set POOLS_DIR in the spike environment}"
: "${SEED_BASE:?set SEED_BASE in the spike environment}"

handoff="$SAMPLE_WORK/metashotgunprep_outputs.env"
[[ -s "$handoff" ]] || { echo "[ERROR] Missing preprocessing handoff: $handoff" >&2; exit 1; }
source "$handoff"
for file in "$CLEAN_R1" "$CLEAN_R2"; do
  [[ -s "$file" ]] || { echo "[ERROR] Missing/empty cleaned mate: $file" >&2; exit 1; }
  gzip -t -- "$file"
done

[[ -s "$POOLS_DIR/pool_files.sha256" ]] || { echo "[ERROR] Pools have not been finalized" >&2; exit 1; }
[[ "$(wc -l < "$POOLS_DIR/pool_files.sha256")" -eq 20 ]] || {
  echo "[ERROR] Expected 20 entries in pool checksum manifest" >&2
  exit 1
}

mkdir -p "$SMOKE_ROOT"
samples="$SMOKE_ROOT/cleaned_sample.tsv"
printf 'sample_id\tR1\tR2\n%s\t%s\t%s\n' "$SAMPLE_ID" "$CLEAN_R1" "$CLEAN_R2" > "$samples"

community="$SMOKE_ROOT/community.tsv"
printf 'label\tweight\n' > "$community"
awk -F'\t' 'NR > 1 && NF {print $1 "\t" ($5 == "" ? 1 : $5)}' "$SPIKE_PANEL" >> "$community"

seed_helper="$PROJECT/spikes/scripts/spikein/stable_seed.py"
independent_work="$SMOKE_ROOT/independent/Fnuc"
community_work="$SMOKE_ROOT/community/CRCpanel"
mkdir -p "$independent_work" "$community_work"
rm -f -- "$SMOKE_ROOT/SUCCESS"

echo "[STEP 1/4] Generate Fnuc independent spike at 0.01%"
env \
  IMG="$IMG" WORK="$independent_work" SAMPLES_TSV="$samples" \
  POOL1="$POOLS_DIR/Fnuc.pool_1.fq" POOL2="$POOLS_DIR/Fnuc.pool_2.fq" \
  LABEL=Fnuc FRACTIONS=0.0001 SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" \
  BIND="${BIND:-}" SAMPLING_MODE="$SAMPLING_MODE" \
  SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
  bash "$PROJECT/spikes/scripts/spikein/spike_one_taxon_array.sbatch"

echo "[STEP 2/4] Generate ten-member community spike at total 0.01%"
env \
  IMG="$IMG" WORK="$community_work" SAMPLES_TSV="$samples" \
  COMMUNITY_TSV="$community" POOLS_DIR="$POOLS_DIR" COMMUNITY_LABEL=CRCpanel \
  FULL_FRACTIONS="${COMMUNITY_FRACTIONS:-0.0001,0.0005,0.001,0.005,0.01,0.05,0.10}" \
  RUN_FRACTIONS=0.0001 SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" \
  BIND="${BIND:-}" SAMPLING_MODE="$SAMPLING_MODE" SLURM_ARRAY_TASK_ID=1 \
  bash "$PROJECT/spikes/scripts/spikein/spike_community_array.sbatch"

ind_id="${SAMPLE_ID}_Fnuc_f0p0001"
com_id="${SAMPLE_ID}_CRCpanel_f0p0001"
ind_r1="$independent_work/spiked_fastqs/${ind_id}_1.fq.gz"
ind_r2="$independent_work/spiked_fastqs/${ind_id}_2.fq.gz"
com_r1="$community_work/spiked_fastqs/${com_id}_1.fq.gz"
com_r2="$community_work/spiked_fastqs/${com_id}_2.fq.gz"
for file in "$ind_r1" "$ind_r2" "$com_r1" "$com_r2"; do
  [[ -s "$file" ]] || { echo "[ERROR] Missing/empty spiked mate: $file" >&2; exit 1; }
  gzip -t -- "$file"
done

echo "[STEP 3/4] Profile independent spike with both profilers"
env \
  YACHIDA_ENV="$YACHIDA_ENV" SAMPLE_ID="$ind_id" SAMPLE_WORK="$SAMPLE_WORK" \
  PROFILE_R1="$ind_r1" PROFILE_R2="$ind_r2" \
  YACHIDA_BASELINE_SMOKE_ROOT="$SMOKE_ROOT/profiles" \
  bash "$PROJECT/datasets/yachida/run_baseline_profiling_smoke.sh"

echo "[STEP 4/4] Profile community spike with both profilers"
env \
  YACHIDA_ENV="$YACHIDA_ENV" SAMPLE_ID="$com_id" SAMPLE_WORK="$SAMPLE_WORK" \
  PROFILE_R1="$com_r1" PROFILE_R2="$com_r2" \
  YACHIDA_BASELINE_SMOKE_ROOT="$SMOKE_ROOT/profiles" \
  bash "$PROJECT/datasets/yachida/run_baseline_profiling_smoke.sh"

provenance="$SMOKE_ROOT/smoke_provenance.tsv"
{
  printf 'field\tvalue\n'
  printf 'sample_id\t%s\n' "$SAMPLE_ID"
  printf 'independent_target\tFnuc\n'
  printf 'independent_fraction\t0.0001\n'
  printf 'community_label\tCRCpanel\n'
  printf 'community_total_fraction\t0.0001\n'
  printf 'pool_checksum_manifest_sha256\t%s\n' "$(sha256sum "$POOLS_DIR/pool_files.sha256" | awk '{print $1}')"
  printf 'seed_scheme\tstable-seed-v1\n'
  printf 'sampling_mode\t%s\n' "$SAMPLING_MODE"
} > "$provenance"

touch "$SMOKE_ROOT/SUCCESS"
echo "[PASS] Yachida spike-generation and dual-profiler smoke test completed"
echo "[OK] Sampling mode: $SAMPLING_MODE"
echo "[KEEP] Synthetic FASTQs retained under: $SMOKE_ROOT"
echo "[OK] Profiles: $SMOKE_ROOT/profiles"
