#!/usr/bin/env bash
# Complete one sample, persist compact outputs, and authorize outer cleanup.
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

: "${PERSISTENT_QC_ROOT:?set in YACHIDA_ENV}"
: "${PERSISTENT_RESULTS_ROOT:?set in YACHIDA_ENV}"
: "${UPSTREAM_SIF:?set in YACHIDA_ENV}"
SPIKE_ENV="${SPIKE_ENV:-$PROJECT/work/yachida_67x3/spikein.env}"
INDEPENDENT_MANIFEST="${INDEPENDENT_MANIFEST:-$PROJECT/work/yachida_67x3/metadata/independent_10_per_condition.tsv}"
[[ -s "$SPIKE_ENV" ]] || { echo "[ERROR] Missing spike environment: $SPIKE_ENV" >&2; exit 1; }
[[ -s "$INDEPENDENT_MANIFEST" ]] || { echo "[ERROR] Missing independent manifest: $INDEPENDENT_MANIFEST" >&2; exit 1; }
source "$SPIKE_ENV"
: "${SPIKE_PANEL:?set in spike environment}"
: "${POOLS_DIR:?set in spike environment}"
: "${SEED_BASE:?set in spike environment}"

SAMPLING_MODE=single_pass
INDEPENDENT_FRACTIONS="${INDEPENDENT_FRACTIONS:-0.0001,0.0005,0.001,0.005,0.01,0.05}"
COMMUNITY_FRACTIONS="${COMMUNITY_FRACTIONS:-0.0001,0.0005,0.001,0.005,0.01,0.05,0.10}"
pair_index="$POOLS_DIR/pool_pair_counts.tsv"
[[ -s "$pair_index" && -s "$POOLS_DIR/pool_pair_counts.tsv.sha256" ]] || {
  echo "[ERROR] Run index_finalized_pools.sh before streaming samples" >&2
  exit 1
}
sha256sum -c "$POOLS_DIR/pool_pair_counts.tsv.sha256"

sample_root="$PERSISTENT_RESULTS_ROOT/$STUDY/$SAMPLE_ID"
profile_root="$sample_root/profiles"
design_root="$sample_root/spike_design"
quarantine_root="$PERSISTENT_RESULTS_ROOT/quarantine/$STUDY/$SAMPLE_ID"
mkdir -p "$profile_root" "$design_root/independent" "$design_root/community"
python3 - "$SAMPLE_WORK" "$sample_root" "$PERSISTENT_QC_ROOT" <<'PY'
import pathlib, sys
scratch, results, qc = (pathlib.Path(value).resolve() for value in sys.argv[1:])
for persistent in (results, qc):
    if persistent == scratch or scratch in persistent.parents:
        raise SystemExit("[ERROR] Persistent output cannot be inside disposable sample scratch")
PY

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
  local profile_group="$1" profile_id="$2" r1="$3" r2="$4"
  local expected="$profile_group/$profile_id"
  if verify_profile "$expected"; then
    echo "[SKIP] Verified profile: $profile_id"
    return 0
  fi
  if [[ -e "$expected" ]]; then
    mkdir -p "$quarantine_root"
    moved="$quarantine_root/${profile_id}.$(date -u +%Y%m%dT%H%M%SZ).$$"
    mv "$expected" "$moved"
    echo "[WARN] Moved incomplete profile to quarantine: $moved" >&2
  fi
  env YACHIDA_ENV="$YACHIDA_ENV" SAMPLE_ID="$profile_id" SAMPLE_WORK="$SAMPLE_WORK" \
    PROFILE_R1="$r1" PROFILE_R2="$r2" YACHIDA_BASELINE_SMOKE_ROOT="$profile_group" \
    bash "$PROJECT/datasets/yachida/run_baseline_profiling_smoke.sh"
  verify_profile "$expected" || { echo "[ERROR] Profile verification failed: $profile_id" >&2; exit 1; }
}

handoff="$SAMPLE_WORK/metashotgunprep_outputs.env"
reuse_cleaned=0
if [[ -s "$handoff" ]]; then
  source "$handoff"
  if [[ -s "${CLEAN_R1:-}" && -s "${CLEAN_R2:-}" ]] && gzip -t "$CLEAN_R1" "$CLEAN_R2"; then
    reuse_cleaned=1
    echo "[SKIP] Reusing validated cleaned mates: $SAMPLE_ID"
  fi
fi
if (( reuse_cleaned == 0 )); then
  bash "$PROJECT/datasets/yachida/run_metashotgunprep.sh"
  source "$handoff"
fi

profile_pair "$profile_root/baseline" "$SAMPLE_ID" "$CLEAN_R1" "$CLEAN_R2"

samples="$SAMPLE_WORK/cleaned_sample.tsv"
printf 'sample_id\tR1\tR2\n%s\t%s\t%s\n' "$SAMPLE_ID" "$CLEAN_R1" "$CLEAN_R2" > "$samples"
seed_helper="$PROJECT/spikes/scripts/spikein/stable_seed.py"

is_independent=0
if python3 - "$INDEPENDENT_MANIFEST" "$SAMPLE_ID" <<'PY'
import csv, pathlib, sys
path, sample = pathlib.Path(sys.argv[1]), sys.argv[2]
with path.open(newline="", encoding="utf-8") as handle:
    matches = [row for row in csv.DictReader(handle, delimiter="\t") if row["sample_id"] == sample]
raise SystemExit(0 if len(matches) == 1 else 1)
PY
then
  is_independent=1
fi

frac_tag() {
  python3 - "$1" <<'PY'
import sys
value = float(sys.argv[1])
text = (f"{value:.6f}").rstrip("0").rstrip(".")
print("f" + text.replace(".", "p"))
PY
}

join_by_comma() {
  local IFS=,
  printf '%s' "$*"
}

if (( is_independent == 1 )); then
  echo "[INFO] $SAMPLE_ID belongs to the nested independent-spike subset"
  while IFS=$'\t' read -r label _taxon _assembly _fasta _weight _url; do
    [[ "$label" == label || -z "$label" ]] && continue
    pool_size="$(awk -F'\t' -v label="$label" '$1 == label {print $4}' "$pair_index")"
    [[ "$pool_size" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Missing pool size: $label" >&2; exit 1; }
    IFS=',' read -r -a fractions <<< "$INDEPENDENT_FRACTIONS"
    missing=()
    for fraction in "${fractions[@]}"; do
      tag="$(frac_tag "$fraction")"
      profile_id="${SAMPLE_ID}_${label}_${tag}"
      profile_group="$profile_root/independent/$label"
      if verify_profile "$profile_group/$profile_id"; then
        echo "[SKIP] Verified independent profile: $profile_id"
      else
        missing+=("$fraction")
      fi
    done
    if (( ${#missing[@]} > 0 )); then
      work="$SAMPLE_WORK/spikes/independent/$label"
      mkdir -p "$work/logs/spiked_manifest_parts"
      if [[ -s "$design_root/independent/${label}.tsv" &&
            ! -s "$work/logs/${SAMPLE_ID}.spike_design.tsv" ]]; then
        cp -f "$design_root/independent/${label}.tsv" "$work/logs/${SAMPLE_ID}.spike_design.tsv"
        echo "[RESUME] Restored persistent independent design: $label"
      fi
      env IMG="$IMG" WORK="$work" SAMPLES_TSV="$samples" \
        POOL1="$POOLS_DIR/${label}.pool_1.fq" POOL2="$POOLS_DIR/${label}.pool_2.fq" \
        POOL_SIZE="$pool_size" LABEL="$label" FRACTIONS="$(join_by_comma "${missing[@]}")" \
        SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" BIND="${BIND:-}" \
        SAMPLING_MODE="$SAMPLING_MODE" SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
        bash "$PROJECT/spikes/scripts/spikein/spike_one_taxon_array.sbatch"
      mkdir -p "$design_root/independent"
      cp -f "$work/logs/${SAMPLE_ID}.spike_design.tsv" "$design_root/independent/${label}.tsv"
      for fraction in "${missing[@]}"; do
        tag="$(frac_tag "$fraction")"
        profile_id="${SAMPLE_ID}_${label}_${tag}"
        profile_group="$profile_root/independent/$label"
        r1="$work/spiked_fastqs/${profile_id}_1.fq.gz"
        r2="$work/spiked_fastqs/${profile_id}_2.fq.gz"
        profile_pair "$profile_group" "$profile_id" "$r1" "$r2"
        rm -f -- "$r1" "$r2"
        echo "[CLEANED] Verified independent spike FASTQs: $profile_id"
      done
    fi
  done < "$SPIKE_PANEL"
else
  echo "[INFO] $SAMPLE_ID is not in the nested independent-spike subset"
fi

community="$SAMPLE_WORK/community.tsv"
printf 'label\tweight\n' > "$community"
awk -F'\t' 'NR > 1 && NF {print $1 "\t" ($5 == "" ? 1 : $5)}' "$SPIKE_PANEL" >> "$community"
IFS=',' read -r -a community_fractions <<< "$COMMUNITY_FRACTIONS"
missing=()
for fraction in "${community_fractions[@]}"; do
  tag="$(frac_tag "$fraction")"
  profile_id="${SAMPLE_ID}_CRCpanel_${tag}"
  profile_group="$profile_root/community"
  if verify_profile "$profile_group/$profile_id"; then
    echo "[SKIP] Verified community profile: $profile_id"
  else
    missing+=("$fraction")
  fi
done
if (( ${#missing[@]} > 0 )); then
  work="$SAMPLE_WORK/spikes/community/CRCpanel"
  mkdir -p "$work/logs/spiked_manifest_parts"
  if [[ -s "$design_root/community/CRCpanel.tsv" &&
        ! -s "$work/logs/${SAMPLE_ID}.CRCpanel.design.tsv" ]]; then
    cp -f "$design_root/community/CRCpanel.tsv" "$work/logs/${SAMPLE_ID}.CRCpanel.design.tsv"
    echo "[RESUME] Restored persistent community design"
  fi
  env IMG="$IMG" WORK="$work" SAMPLES_TSV="$samples" COMMUNITY_TSV="$community" \
    POOLS_DIR="$POOLS_DIR" COMMUNITY_LABEL=CRCpanel FULL_FRACTIONS="$COMMUNITY_FRACTIONS" \
    RUN_FRACTIONS="$(join_by_comma "${missing[@]}")" SEED_BASE="$SEED_BASE" SEED_HELPER="$seed_helper" \
    BIND="${BIND:-}" SAMPLING_MODE="$SAMPLING_MODE" SLURM_ARRAY_TASK_ID=1 KEEP_TMP=0 \
    bash "$PROJECT/spikes/scripts/spikein/spike_community_array.sbatch"
  cp -f "$work/logs/${SAMPLE_ID}.CRCpanel.design.tsv" "$design_root/community/CRCpanel.tsv"
  for fraction in "${missing[@]}"; do
    tag="$(frac_tag "$fraction")"
    profile_id="${SAMPLE_ID}_CRCpanel_${tag}"
    profile_group="$profile_root/community"
    r1="$work/spiked_fastqs/${profile_id}_1.fq.gz"
    r2="$work/spiked_fastqs/${profile_id}_2.fq.gz"
    profile_pair "$profile_group" "$profile_id" "$r1" "$r2"
    rm -f -- "$r1" "$r2"
    echo "[CLEANED] Verified community spike FASTQs: $profile_id"
  done
fi

expected_profiles=8
(( is_independent == 1 )) && expected_profiles=68
observed_profiles="$(find "$profile_root" -name SUCCESS -type f | wc -l)"
[[ "$observed_profiles" -eq "$expected_profiles" ]] || {
  echo "[ERROR] Expected $expected_profiles successful profiles; found $observed_profiles" >&2
  exit 1
}
while IFS= read -r directory; do verify_profile "$directory"; done \
  < <(find "$profile_root" -name SUCCESS -type f -printf '%h\n' | sort)
community_design_rows="$(awk 'NR > 1 {n++} END {print n+0}' "$design_root/community/CRCpanel.tsv")"
[[ "$community_design_rows" -eq 7 ]] || {
  echo "[ERROR] Expected 7 community design rows; found $community_design_rows" >&2
  exit 1
}
independent_design_rows=0
if (( is_independent == 1 )); then
  independent_design_rows="$(awk 'FNR > 1 {n++} END {print n+0}' "$design_root"/independent/*.tsv)"
  [[ "$independent_design_rows" -eq 60 ]] || {
    echo "[ERROR] Expected 60 independent design rows; found $independent_design_rows" >&2
    exit 1
  }
fi

printf 'field\tvalue\n' > "$sample_root/sample_completion.tsv"
printf '%s\t%s\n' \
  sample_id "$SAMPLE_ID" \
  study "$STUDY" \
  condition "$TARGET_CONDITION" \
  batch_id "${BATCH_ID:-}" \
  batch_position "${BATCH_POSITION:-}" \
  sampling_mode "$SAMPLING_MODE" \
  independent_subset "$is_independent" \
  expected_profiles "$expected_profiles" \
  observed_profiles "$observed_profiles" \
  community_design_rows "$community_design_rows" \
  independent_design_rows "$independent_design_rows" \
  >> "$sample_root/sample_completion.tsv"
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
print(f"[OK] Recorded {len(files)} persistent outputs")
PY
echo "[PASS] Full streaming sample completed: $SAMPLE_ID"
echo "[OK] Persistent results: $sample_root"
