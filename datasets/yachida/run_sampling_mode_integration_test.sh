#!/usr/bin/env bash
# Compare legacy and single-pass modes through the production spike scripts.
set -euo pipefail
IFS=$'\n\t'
umask 002

: "${PROJECT:?set PROJECT to repository root}"
: "${YACHIDA_ENV:?set YACHIDA_ENV}"
: "${SAMPLE_ID:?set SAMPLE_ID}"
source "$YACHIDA_ENV"
SPIKE_ENV="${SPIKE_ENV:-$PROJECT/work/yachida_67x3/spikein.env}"
source "$SPIKE_ENV"

: "${IMG:?set IMG in spike environment}"
: "${POOLS_DIR:?set POOLS_DIR in spike environment}"
: "${SPIKE_PANEL:?set SPIKE_PANEL in spike environment}"
: "${SEED_BASE:?set SEED_BASE in spike environment}"

SAMPLE_WORK="${SAMPLE_WORK:-$PROJECT/work/yachida_67x3/smoke_scratch/$SAMPLE_ID}"
source "$SAMPLE_WORK/metashotgunprep_outputs.env"
TEST_ROOT="${YACHIDA_SAMPLING_TEST_ROOT:-$PROJECT/work/yachida_67x3/sampling_mode_integration/$SAMPLE_ID}"
BACKGROUND_PAIRS="${BACKGROUND_PAIRS:-100000}"
POOL_PAIRS="${POOL_PAIRS:-1000}"
FRACTIONS="${TEST_FRACTIONS:-0.0001,0.0005,0.001}"
seed_helper="$PROJECT/spikes/scripts/spikein/stable_seed.py"

rm -rf -- "$TEST_ROOT.tmp"
mkdir -p "$TEST_ROOT.tmp/background" "$TEST_ROOT.tmp/pools"
root="$TEST_ROOT.tmp"

python3 - "$CLEAN_R1" "$CLEAN_R2" "$BACKGROUND_PAIRS" "$root/background" <<'PY'
import gzip, pathlib, sys
r1, r2, pairs, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), pathlib.Path(sys.argv[4])
for mate, source in enumerate((r1, r2), 1):
    with gzip.open(source, "rt") as inp, gzip.open(out / f"background_{mate}.fq.gz", "wt") as dst:
        for index, line in enumerate(inp):
            if index >= pairs * 4:
                break
            dst.write(line)
PY

while IFS=$'\t' read -r label _taxon _assembly _fasta _weight _url; do
  [[ "$label" == label || -z "$label" ]] && continue
  for mate in 1 2; do
    head -n "$((POOL_PAIRS * 4))" "$POOLS_DIR/${label}.pool_${mate}.fq" \
      > "$root/pools/${label}.pool_${mate}.fq"
  done
done < "$SPIKE_PANEL"

samples="$root/samples.tsv"
printf 'sample_id\tR1\tR2\n%s\t%s\t%s\n' \
  "${SAMPLE_ID}_integration" "$root/background/background_1.fq.gz" \
  "$root/background/background_2.fq.gz" > "$samples"
community="$root/community.tsv"
printf 'label\tweight\n' > "$community"
awk -F'\t' 'NR > 1 && NF {print $1 "\t" ($5 == "" ? 1 : $5)}' "$SPIKE_PANEL" >> "$community"

for mode in legacy single_pass; do
  echo "[RUN] Independent production script: $mode"
  env IMG="$IMG" WORK="$root/$mode/independent/Fnuc" SAMPLES_TSV="$samples" \
    POOL1="$root/pools/Fnuc.pool_1.fq" POOL2="$root/pools/Fnuc.pool_2.fq" \
    LABEL=Fnuc FRACTIONS="$FRACTIONS" SEED_BASE="$SEED_BASE" \
    SEED_HELPER="$seed_helper" BIND="${BIND:-}" SAMPLING_MODE="$mode" \
    FASTQ_ASSEMBLY_MODE=recompress SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
    bash "$PROJECT/spikes/scripts/spikein/spike_one_taxon_array.sbatch"

  echo "[RUN] Community production script: $mode"
  env IMG="$IMG" WORK="$root/$mode/community/CRCpanel" SAMPLES_TSV="$samples" \
    COMMUNITY_TSV="$community" POOLS_DIR="$root/pools" COMMUNITY_LABEL=CRCpanel \
    FULL_FRACTIONS="$FRACTIONS" RUN_FRACTIONS="$FRACTIONS" \
    SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" BIND="${BIND:-}" \
    SAMPLING_MODE="$mode" FASTQ_ASSEMBLY_MODE=recompress SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
    bash "$PROJECT/spikes/scripts/spikein/spike_community_array.sbatch"
done

echo "[RUN] Independent production script: single_pass + gzip_members"
env IMG="$IMG" WORK="$root/gzip_members/independent/Fnuc" SAMPLES_TSV="$samples" \
  POOL1="$root/pools/Fnuc.pool_1.fq" POOL2="$root/pools/Fnuc.pool_2.fq" \
  LABEL=Fnuc FRACTIONS="$FRACTIONS" SEED_BASE="$SEED_BASE" \
  SEED_HELPER="$seed_helper" BIND="${BIND:-}" SAMPLING_MODE=single_pass \
  FASTQ_ASSEMBLY_MODE=gzip_members SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
  bash "$PROJECT/spikes/scripts/spikein/spike_one_taxon_array.sbatch"

echo "[RUN] Community production script: single_pass + gzip_members"
env IMG="$IMG" WORK="$root/gzip_members/community/CRCpanel" SAMPLES_TSV="$samples" \
  COMMUNITY_TSV="$community" POOLS_DIR="$root/pools" COMMUNITY_LABEL=CRCpanel \
  FULL_FRACTIONS="$FRACTIONS" RUN_FRACTIONS="$FRACTIONS" \
  SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" BIND="${BIND:-}" \
  SAMPLING_MODE=single_pass FASTQ_ASSEMBLY_MODE=gzip_members \
  SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
  bash "$PROJECT/spikes/scripts/spikein/spike_community_array.sbatch"

compare_fastqs() {
  local legacy_root="$1" single_root="$2" count=0 file relative other left right
  while IFS= read -r file; do
    relative="${file#"$legacy_root/"}"
    other="$single_root/$relative"
    [[ -s "$other" ]] || { echo "[ERROR] Missing single-pass output: $other" >&2; return 1; }
    left="$(gzip -cd "$file" | sha256sum | cut -d' ' -f1)"
    right="$(gzip -cd "$other" | sha256sum | cut -d' ' -f1)"
    [[ "$left" == "$right" ]] || { echo "[ERROR] Sequence mismatch: $relative" >&2; return 1; }
    count=$((count + 1))
  done < <(find "$legacy_root" -type f -name '*.fq.gz' | sort)
  [[ "$count" -gt 0 ]] || { echo "[ERROR] No legacy FASTQs found" >&2; return 1; }
  echo "[OK] $count compressed FASTQ products have identical decompressed content"
}

compare_fastqs "$root/legacy/independent/Fnuc/spiked_fastqs" "$root/single_pass/independent/Fnuc/spiked_fastqs"
compare_fastqs "$root/legacy/community/CRCpanel/spiked_fastqs" "$root/single_pass/community/CRCpanel/spiked_fastqs"
compare_fastqs "$root/single_pass/independent/Fnuc/spiked_fastqs" "$root/gzip_members/independent/Fnuc/spiked_fastqs"
compare_fastqs "$root/single_pass/community/CRCpanel/spiked_fastqs" "$root/gzip_members/community/CRCpanel/spiked_fastqs"

cut -f1-10 "$root/legacy/independent/Fnuc/logs/${SAMPLE_ID}_integration.spike_design.tsv" > "$root/legacy.independent.key.tsv"
cut -f1-10 "$root/single_pass/independent/Fnuc/logs/${SAMPLE_ID}_integration.spike_design.tsv" > "$root/single.independent.key.tsv"
cmp "$root/legacy.independent.key.tsv" "$root/single.independent.key.tsv"
cmp "$root/legacy/community/CRCpanel/logs/${SAMPLE_ID}_integration.CRCpanel.design.tsv" \
    "$root/single_pass/community/CRCpanel/logs/${SAMPLE_ID}_integration.CRCpanel.design.tsv"
cut -f1-10 "$root/gzip_members/independent/Fnuc/logs/${SAMPLE_ID}_integration.spike_design.tsv" > "$root/gzip_members.independent.key.tsv"
cmp "$root/single.independent.key.tsv" "$root/gzip_members.independent.key.tsv"
cmp "$root/single_pass/community/CRCpanel/logs/${SAMPLE_ID}_integration.CRCpanel.design.tsv" \
    "$root/gzip_members/community/CRCpanel/logs/${SAMPLE_ID}_integration.CRCpanel.design.tsv"

if [[ "${RUN_PROFILE_EQUIVALENCE:-0}" == 1 ]]; then
  echo "[RUN] Profiler equivalence: recompress versus gzip_members"
  tag="$(python3 - "${FRACTIONS%%,*}" <<'PY'
import sys
value = float(sys.argv[1])
print("f" + (f"{value:.6f}").rstrip("0").rstrip(".").replace(".", "p"))
PY
)"
  for assembly in single_pass gzip_members; do
    profile_id="assembly_${assembly}"
    input="$root/$assembly/independent/Fnuc/spiked_fastqs/${SAMPLE_ID}_integration_Fnuc_${tag}"
    env YACHIDA_ENV="$YACHIDA_ENV" SAMPLE_ID="$profile_id" SAMPLE_WORK="$root" \
      PROFILE_R1="${input}_1.fq.gz" PROFILE_R2="${input}_2.fq.gz" \
      YACHIDA_BASELINE_SMOKE_ROOT="$root/profiles/$assembly" \
      bash "$PROJECT/datasets/yachida/run_baseline_profiling_smoke.sh"
  done
  left="$root/profiles/single_pass/assembly_single_pass"
  right="$root/profiles/gzip_members/assembly_gzip_members"
  cmp "$left/assembly_single_pass.kraken2.report" "$right/assembly_gzip_members.kraken2.report"
  cmp "$left/assembly_single_pass.bracken.S.tsv" "$right/assembly_gzip_members.bracken.S.tsv"
  grep -v '^#' "$left/assembly_single_pass.metaphlan.tsv" > "$root/metaphlan.recompress.data.tsv"
  grep -v '^#' "$right/assembly_gzip_members.metaphlan.tsv" > "$root/metaphlan.gzip_members.data.tsv"
  cmp "$root/metaphlan.recompress.data.tsv" "$root/metaphlan.gzip_members.data.tsv"
  echo "[OK] Both profilers produced identical biological result tables"
fi

printf 'field\tvalue\nbackground_pairs\t%s\npool_pairs\t%s\nfractions\t%s\ncomparison\tdecompressed_sha256\n' \
  "$BACKGROUND_PAIRS" "$POOL_PAIRS" "$FRACTIONS" > "$root/integration_test.tsv"
touch "$root/SUCCESS"
rm -rf -- "$TEST_ROOT"
mv "$root" "$TEST_ROOT"
echo "[PASS] Sampling and FASTQ-assembly modes produce sequence-identical outputs"
echo "[OK] Results: $TEST_ROOT"
