#!/usr/bin/env bash
# Compare legacy repeated seqtk scans with an exact single-read fan-out.
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  benchmark_single_pass_sampling.sh --env spikein.env [options]

Options:
  --label LABEL       panel label (default: Fnuc)
  --sample-id ID      stable seed identity (default: sampling_benchmark)
  --subset-pairs N    paired records copied for the benchmark (default: 100000)
  --outdir DIR        output directory (default: WORK/sampling_benchmark)

This benchmark never changes production pools. It creates a small paired prefix,
selects six subsets with the existing stable seeds by (a) scanning repeatedly
and (b) broadcasting one pool read to concurrent seqtk processes, and requires
byte-identical R1/R2 selections for every fraction.
EOF
}

ENV=""
LABEL="Fnuc"
SAMPLE_ID="sampling_benchmark"
SUBSET_PAIRS=100000
OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --subset-pairs) SUBSET_PAIRS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$ENV" && -s "$ENV" ]] || { usage; exit 2; }
source "$ENV"
: "${IMG:?set IMG in env}"
: "${WORK:?set WORK in env}"
: "${POOLS_DIR:?set POOLS_DIR in env}"
: "${SEED_BASE:?set SEED_BASE in env}"

[[ "$SUBSET_PAIRS" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid --subset-pairs" >&2; exit 2; }
OUTDIR="${OUTDIR:-$WORK/sampling_benchmark/$LABEL}"
pool1="$POOLS_DIR/${LABEL}.pool_1.fq"
pool2="$POOLS_DIR/${LABEL}.pool_2.fq"
[[ -s "$pool1" && -s "$pool2" ]] || { echo "[ERROR] Missing finalized pool for $LABEL" >&2; exit 1; }
command -v apptainer >/dev/null 2>&1 || { echo "[ERROR] apptainer not found" >&2; exit 1; }

# Apptainer bind destinations must be absolute. Accept a convenient relative
# --outdir from users, but canonicalize it before constructing any bind.
mkdir -p "$(dirname "$OUTDIR")"
OUTDIR="$(cd "$(dirname "$OUTDIR")" && pwd -P)/$(basename "$OUTDIR")"

rm -rf -- "$OUTDIR.tmp"
mkdir -p "$OUTDIR.tmp"/{subset,legacy,single,tmp}
bench="$OUTDIR.tmp"
lines=$((SUBSET_PAIRS * 4))

echo "[STEP 1/5] Copy first $SUBSET_PAIRS synchronized pool pairs"
head -n "$lines" "$pool1" > "$bench/subset/pool_1.fq"
head -n "$lines" "$pool2" > "$bench/subset/pool_2.fq"
[[ "$(awk 'END{print int(NR/4)}' "$bench/subset/pool_1.fq")" -eq "$SUBSET_PAIRS" ]]
[[ "$(awk 'END{print int(NR/4)}' "$bench/subset/pool_2.fq")" -eq "$SUBSET_PAIRS" ]]

fractions=(0.0001 0.0005 0.001 0.005 0.01 0.05)
# Selection sizes span the useful benchmark range and remain below SUBSET_PAIRS.
sizes=(100 500 1000 5000 10000 50000)
manifest="$bench/selection_manifest.tsv"
printf 'tag\tfraction\tseed\tn\n' > "$manifest"
for i in "${!fractions[@]}"; do
  f="${fractions[$i]}"
  n="${sizes[$i]}"
  (( n <= SUBSET_PAIRS )) || { echo "[ERROR] subset-pairs must be at least $n" >&2; exit 1; }
  canonical="$(python3 - "$f" <<'PY'
import decimal, sys
print(format(decimal.Decimal(sys.argv[1]).normalize(), "f"))
PY
)"
  seed="$(python3 "$(dirname "${BASH_SOURCE[0]}")/stable_seed.py" --base "$SEED_BASE" \
    --namespace spike-independent-v1 "$SAMPLE_ID" "$LABEL" "$canonical")"
  tag="f$(printf '%s' "$f" | sed 's/[.]\([0-9]*\)$/p\1/')"
  printf '%s\t%s\t%s\t%s\n' "$tag" "$f" "$seed" "$n" >> "$manifest"
done

extract_mate2() {
  local selected1="$1" pool_r2="$2" output2="$3" stem="$4"
  local ids1="$bench/tmp/${stem}.ids1" ids2="$bench/tmp/${stem}.ids2"
  awk 'NR%4==1 {gsub(/^@/,"",$1); print $1}' "$selected1" > "$ids1"
  awk '{o=$1; print o; id=o; if(id~/\/1$/){sub(/\/1$/,"/2",id); print id} else if(id~/_1$/){sub(/_1$/,"_2",id); print id} else if(id~/[.]1$/){sub(/[.]1$/,".2",id); print id}}' \
    "$ids1" | awk '!seen[$0]++' > "$ids2"
  apptainer exec --cleanenv -B "$bench:$bench" "$IMG" seqtk subseq "$pool_r2" "$ids2" > "$output2"
}

echo "[STEP 2/5] Legacy repeated scans"
legacy_start=$(date +%s)
while IFS=$'\t' read -r tag _fraction seed n; do
  [[ "$tag" == "tag" ]] && continue
  apptainer exec --cleanenv -B "$bench:$bench" "$IMG" seqtk sample -s "$seed" "$bench/subset/pool_1.fq" "$n" > "$bench/legacy/${tag}_1.fq"
  extract_mate2 "$bench/legacy/${tag}_1.fq" "$bench/subset/pool_2.fq" "$bench/legacy/${tag}_2.fq" "legacy.$tag"
done < "$manifest"
legacy_seconds=$(( $(date +%s) - legacy_start ))

echo "[STEP 3/5] One-read fan-out scans"
single_start=$(date +%s)
declare -a fifos=() pids=() tags=()
while IFS=$'\t' read -r tag _fraction seed n; do
  [[ "$tag" == "tag" ]] && continue
  fifo="$bench/tmp/${tag}.r1.fifo"; mkfifo "$fifo"; fifos+=("$fifo"); tags+=("$tag")
  apptainer exec --cleanenv -B "$bench:$bench" "$IMG" seqtk sample -s "$seed" "$fifo" "$n" > "$bench/single/${tag}_1.fq" &
  pids+=("$!")
done < "$manifest"
tee "${fifos[@]}" < "$bench/subset/pool_1.fq" > /dev/null
for pid in "${pids[@]}"; do wait "$pid"; done

fifos=(); pids=()
for tag in "${tags[@]}"; do
  ids1="$bench/tmp/single.${tag}.ids1"; ids2="$bench/tmp/single.${tag}.ids2"
  awk 'NR%4==1 {gsub(/^@/,"",$1); print $1}' "$bench/single/${tag}_1.fq" > "$ids1"
  awk '{o=$1; print o; id=o; if(id~/\/1$/){sub(/\/1$/,"/2",id); print id} else if(id~/_1$/){sub(/_1$/,"_2",id); print id} else if(id~/[.]1$/){sub(/[.]1$/,".2",id); print id}}' \
    "$ids1" | awk '!seen[$0]++' > "$ids2"
  fifo="$bench/tmp/${tag}.r2.fifo"; mkfifo "$fifo"; fifos+=("$fifo")
  apptainer exec --cleanenv -B "$bench:$bench" "$IMG" seqtk subseq "$fifo" "$ids2" > "$bench/single/${tag}_2.fq" &
  pids+=("$!")
done
tee "${fifos[@]}" < "$bench/subset/pool_2.fq" > /dev/null
for pid in "${pids[@]}"; do wait "$pid"; done
single_seconds=$(( $(date +%s) - single_start ))

echo "[STEP 4/5] Require byte-identical selections"
comparison="$bench/comparison.tsv"
printf 'tag\tmate\tlegacy_sha256\tsingle_sha256\tidentical\n' > "$comparison"
for tag in "${tags[@]}"; do
  for mate in 1 2; do
    a="$bench/legacy/${tag}_${mate}.fq"; b="$bench/single/${tag}_${mate}.fq"
    sha_a="$(sha256sum "$a" | awk '{print $1}')"; sha_b="$(sha256sum "$b" | awk '{print $1}')"
    identical=false; [[ "$sha_a" == "$sha_b" ]] && cmp -s "$a" "$b" && identical=true
    printf '%s\t%s\t%s\t%s\t%s\n' "$tag" "$mate" "$sha_a" "$sha_b" "$identical" >> "$comparison"
    [[ "$identical" == true ]] || { echo "[FAIL] Non-identical output for $tag mate $mate" >&2; exit 1; }
  done
done

echo "[STEP 5/5] Save benchmark report"
{
  printf 'metric\tvalue\n'
  printf 'label\t%s\n' "$LABEL"
  printf 'subset_pairs\t%s\n' "$SUBSET_PAIRS"
  printf 'fractions\t%s\n' "${#fractions[@]}"
  printf 'legacy_seconds\t%s\n' "$legacy_seconds"
  printf 'single_pass_seconds\t%s\n' "$single_seconds"
  python3 - "$legacy_seconds" "$single_seconds" <<'PY'
import sys
old, new = map(float, sys.argv[1:])
print(f"speedup\t{old/new:.3f}" if new else "speedup\tinf")
PY
} > "$bench/timing.tsv"

rm -rf -- "$bench/tmp"
rm -rf -- "$OUTDIR"
mv "$bench" "$OUTDIR"
echo "[PASS] Single-pass sampling is byte-identical for all 12 selected mate files"
cat "$OUTDIR/timing.tsv"
echo "[OK] Results: $OUTDIR"
