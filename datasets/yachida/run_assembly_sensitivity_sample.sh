#!/usr/bin/env bash
# Add two clean-assembly independent-spike arms without replacing production.
set -euo pipefail
IFS=$'\n\t'
umask 002

: "${PROJECT:?set PROJECT}"
: "${YACHIDA_ENV:?set YACHIDA_ENV}"
: "${SAMPLE_ID:?provided by stream_sample.py}"
: "${SAMPLE_WORK:?provided by stream_sample.py}"
: "${RECEIPT:?provided by stream_sample.py}"
: "${STUDY:?provided by stream_sample.py}"
: "${TARGET_CONDITION:?provided by stream_sample.py}"
source "$YACHIDA_ENV"

: "${UPSTREAM_SIF:?set in YACHIDA_ENV}"
: "${PERSISTENT_QC_ROOT:?set in YACHIDA_ENV}"
: "${ASSEMBLY_SENSITIVITY_ROOT:?set in YACHIDA_ENV}"
: "${ASSEMBLY_SENSITIVITY_POOLS_DIR:?set in YACHIDA_ENV}"
ASSEMBLY_SENSITIVITY_ARMS="${ASSEMBLY_SENSITIVITY_ARMS:-$PROJECT/datasets/yachida/assembly_sensitivity_arms.tsv}"
SPIKE_ENV="${SPIKE_ENV:-$PROJECT/work/yachida_67x3/spikein.env}"
INDEPENDENT_MANIFEST="${INDEPENDENT_MANIFEST:-$PROJECT/work/yachida_67x3/metadata/independent_10_per_condition.tsv}"
[[ -s "$ASSEMBLY_SENSITIVITY_ARMS" ]] || { echo "[ERROR] Missing arm manifest" >&2; exit 1; }
[[ -s "$SPIKE_ENV" ]] || { echo "[ERROR] Missing spike environment" >&2; exit 1; }
[[ -s "$INDEPENDENT_MANIFEST" ]] || { echo "[ERROR] Missing independent manifest" >&2; exit 1; }
source "$SPIKE_ENV"

FASTQ_ASSEMBLY_MODE="${FASTQ_ASSEMBLY_MODE:-gzip_members}"
PROFILE_CONCURRENCY="${PROFILE_CONCURRENCY:-1}"
BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-10}"
INDEPENDENT_FRACTIONS="${INDEPENDENT_FRACTIONS:-0.0001,0.0005,0.001,0.005,0.01,0.05}"
[[ "$BRACKEN_THRESHOLD" == 10 ]] || { echo "[ERROR] BRACKEN_THRESHOLD must be 10" >&2; exit 1; }
[[ "$PROFILE_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid PROFILE_CONCURRENCY" >&2; exit 1; }

python3 - "$INDEPENDENT_MANIFEST" "$SAMPLE_ID" <<'PY'
import csv, pathlib, sys
path, sample = pathlib.Path(sys.argv[1]), sys.argv[2]
with path.open(newline="", encoding="utf-8") as handle:
    matches = [row for row in csv.DictReader(handle, delimiter="\t") if row["sample_id"] == sample]
if len(matches) != 1:
    raise SystemExit(f"[ERROR] {sample} is not uniquely present in the frozen independent subset")
PY

pair_index="$ASSEMBLY_SENSITIVITY_POOLS_DIR/pool_pair_counts.tsv"
seal="$ASSEMBLY_SENSITIVITY_POOLS_DIR/pool_pair_counts.tsv.sha256"
[[ -s "$pair_index" && -s "$seal" ]] || { echo "[ERROR] Missing sealed replacement-pool index" >&2; exit 1; }
(cd "$ASSEMBLY_SENSITIVITY_POOLS_DIR" && sha256sum -c "$(basename "$seal")")

sample_root="$ASSEMBLY_SENSITIVITY_ROOT/results/$STUDY/$SAMPLE_ID"
profile_root="$sample_root/profiles/independent"
design_root="$sample_root/spike_design/independent"
quarantine_root="$ASSEMBLY_SENSITIVITY_ROOT/quarantine/$STUDY/$SAMPLE_ID"
mkdir -p "$profile_root" "$design_root"

verify_profile() {
  python3 - "$1" <<'PY'
import csv, hashlib, pathlib, sys
root = pathlib.Path(sys.argv[1])
receipt = root / "retained_outputs.tsv"
if not (root / "SUCCESS").is_file() or not receipt.is_file():
    raise SystemExit(1)
with receipt.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if not rows:
    raise SystemExit(1)
for row in rows:
    path = pathlib.Path(row["path"])
    if not path.is_file() or path.stat().st_size != int(row["bytes"]):
        raise SystemExit(1)
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != row["sha256"]:
        raise SystemExit(1)
PY
}

profile_pair() {
  local group="$1" profile_id="$2" r1="$3" r2="$4" expected
  expected="$group/$profile_id"
  if verify_profile "$expected"; then echo "[SKIP] Verified profile: $profile_id"; return; fi
  if [[ -e "$expected" ]]; then
    mkdir -p "$quarantine_root"
    mv "$expected" "$quarantine_root/${profile_id}.$(date -u +%Y%m%dT%H%M%SZ).$$"
  fi
  env YACHIDA_ENV="$YACHIDA_ENV" SAMPLE_ID="$profile_id" SAMPLE_WORK="$SAMPLE_WORK" \
    PROFILE_R1="$r1" PROFILE_R2="$r2" YACHIDA_BASELINE_SMOKE_ROOT="$group" \
    bash "$PROJECT/datasets/yachida/run_baseline_profiling_smoke.sh"
  verify_profile "$expected" || { echo "[ERROR] Profile verification failed: $profile_id" >&2; exit 1; }
}

handoff="$SAMPLE_WORK/metashotgunprep_outputs.env"
if [[ -s "$handoff" ]]; then
  source "$handoff"
fi
if [[ ! -s "${CLEAN_R1:-}" || ! -s "${CLEAN_R2:-}" ]] || ! gzip -t "$CLEAN_R1" "$CLEAN_R2"; then
  bash "$PROJECT/datasets/yachida/run_metashotgunprep.sh"
  source "$handoff"
fi

samples="$SAMPLE_WORK/cleaned_sample.tsv"
printf 'sample_id\tR1\tR2\n%s\t%s\t%s\n' "$SAMPLE_ID" "$CLEAN_R1" "$CLEAN_R2" > "$samples"
seed_helper="$PROJECT/spikes/scripts/spikein/stable_seed.py"

frac_tag() {
  python3 - "$1" <<'PY'
import sys
text = (f"{float(sys.argv[1]):.6f}").rstrip("0").rstrip(".")
print("f" + text.replace(".", "p"))
PY
}

arm_count=0
while IFS=$'\t' read -r arm_label target_label assembly_accession extra; do
  [[ "$arm_label" == arm_label || -z "$arm_label" ]] && continue
  [[ "$arm_label" =~ ^[A-Za-z0-9._-]+$ && "$target_label" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "[ERROR] Unsafe assembly-arm label" >&2; exit 1;
  }
  pool_size="$(awk -F'\t' -v label="$target_label" '$1 == label {print $4}' "$pair_index")"
  [[ "$pool_size" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Missing pool size: $target_label" >&2; exit 1; }
  IFS=',' read -r -a fractions <<< "$INDEPENDENT_FRACTIONS"
  missing=()
  for fraction in "${fractions[@]}"; do
    tag="$(frac_tag "$fraction")"
    profile_id="${SAMPLE_ID}_${arm_label}_${tag}"
    if ! verify_profile "$profile_root/$arm_label/$profile_id"; then missing+=("$fraction"); fi
  done
  if (( ${#missing[@]} > 0 )); then
    work="$SAMPLE_WORK/spikes/assembly_sensitivity/$arm_label"
    mkdir -p "$work/logs/spiked_manifest_parts"
    fractions_csv="$(IFS=,; printf '%s' "${missing[*]}")"
    env IMG="$IMG" WORK="$work" SAMPLES_TSV="$samples" \
      POOL1="$ASSEMBLY_SENSITIVITY_POOLS_DIR/${target_label}.pool_1.fq" \
      POOL2="$ASSEMBLY_SENSITIVITY_POOLS_DIR/${target_label}.pool_2.fq" \
      POOL_SIZE="$pool_size" LABEL="$target_label" OUTPUT_LABEL="$arm_label" SEED_LABEL="$target_label" \
      FRACTIONS="$fractions_csv" SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" \
      BIND="${BIND:-}" SAMPLING_MODE=single_pass FASTQ_ASSEMBLY_MODE="$FASTQ_ASSEMBLY_MODE" \
      SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
      bash "$PROJECT/spikes/scripts/spikein/spike_one_taxon_array.sbatch"
    cp -f "$work/logs/${SAMPLE_ID}.spike_design.tsv" "$design_root/${arm_label}.tsv"
    for fraction in "${missing[@]}"; do
      tag="$(frac_tag "$fraction")"
      profile_id="${SAMPLE_ID}_${arm_label}_${tag}"
      r1="$work/spiked_fastqs/${profile_id}_1.fq.gz"
      r2="$work/spiked_fastqs/${profile_id}_2.fq.gz"
      profile_pair "$profile_root/$arm_label" "$profile_id" "$r1" "$r2"
      rm -f -- "$r1" "$r2"
    done
  fi
  ((arm_count+=1))
done < "$ASSEMBLY_SENSITIVITY_ARMS"

fraction_count="$(awk -F',' '{print NF}' <<< "$INDEPENDENT_FRACTIONS")"
expected_profiles=$((arm_count * fraction_count))
observed_profiles="$(find "$profile_root" -name SUCCESS -type f | wc -l)"
[[ "$arm_count" -eq 2 && "$expected_profiles" -eq 12 && "$observed_profiles" -eq 12 ]] || {
  echo "[ERROR] Expected 12 clean-assembly profiles; observed $observed_profiles" >&2; exit 1;
}
while IFS= read -r directory; do verify_profile "$directory"; done \
  < <(find "$profile_root" -name SUCCESS -type f -printf '%h\n' | sort)

cp -f "$ASSEMBLY_SENSITIVITY_ARMS" "$sample_root/assembly_sensitivity_arms.tsv"
printf 'field\tvalue\nsample_id\t%s\nstudy\t%s\ncondition\t%s\nexpected_profiles\t12\nobserved_profiles\t%s\n' \
  "$SAMPLE_ID" "$STUDY" "$TARGET_CONDITION" "$observed_profiles" > "$sample_root/sample_completion.tsv"
touch "$sample_root/SUCCESS"

python3 - "$sample_root" "$PERSISTENT_QC_ROOT/$STUDY/$SAMPLE_ID" "$RECEIPT" <<'PY'
import hashlib, pathlib, sys
roots = [pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve()]
receipt = pathlib.Path(sys.argv[3]).resolve()
files = sorted({path for root in roots for path in root.rglob("*") if path.is_file() and path.stat().st_size})
if not files:
    raise SystemExit("[ERROR] No persistent outputs found")
receipt.parent.mkdir(parents=True, exist_ok=True)
with receipt.open("w", encoding="utf-8") as handle:
    handle.write("path\tsha256\tbytes\n")
    for path in files:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                digest.update(block)
        handle.write(f"{path}\t{digest.hexdigest()}\t{path.stat().st_size}\n")
print(f"[OK] Recorded {len(files)} persistent sensitivity outputs")
PY
echo "[PASS] Assembly-sensitivity sample completed: $SAMPLE_ID"
