#!/usr/bin/env bash
#SBATCH -J spike_one_taxon
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 08:00:00

set -euo pipefail
IFS=$'\n\t'
umask 002

: "${IMG:?set IMG}"
: "${WORK:?set WORK (output dir for this taxon)}"
: "${SAMPLES_TSV:?set SAMPLES_TSV}"
: "${POOL1:?set POOL1}"
: "${POOL2:?set POOL2}"
: "${LABEL:?set LABEL}"
: "${FRACTIONS:?set FRACTIONS (comma-separated)}"
: "${SEED_BASE:?set SEED_BASE}"
: "${SEED_HELPER:?set SEED_HELPER to stable_seed.py}"
POOL_SIZE="${POOL_SIZE:-}"
MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"
KEEP_TMP="${KEEP_TMP:-0}"
BIND="${BIND:-}"

if command -v singularity >/dev/null 2>&1; then
  CTR="singularity"
elif command -v apptainer >/dev/null 2>&1; then
  CTR="apptainer"
else
  echo "[ERROR] singularity/apptainer not found on compute node" >&2
  exit 1
fi

run_in_container() {
  local cmd="$1"
  local cache="$WORK/tmp/xdg_cache_${SLURM_JOB_ID:-nojid}_${SLURM_ARRAY_TASK_ID:-notask}"
  mkdir -p "$cache"
  local bind_args=()
  [[ -n "$BIND" ]] && bind_args=(-B "$BIND")

  "$CTR" exec \
    "${bind_args[@]}" \
    -B "$WORK:$WORK" \
    --env XDG_CACHE_HOME="$cache" \
    "$IMG" bash -lc "$cmd"
}

count_reads() {
  local fq="$1"
  if [[ "$fq" == *.gz ]]; then
    gzip -cd -- "$fq" | awk 'END{print int(NR/4)}'
  else
    awk 'END{print int(NR/4)}' "$fq"
  fi
}

extract_first_n_lines_gz() {
  local in="$1"
  local n="$2"
  local out="$3"

  python3 - "$in" "$n" "$out" <<'PY'
import gzip, sys
inp, n_raw, outp = sys.argv[1:4]
n = int(n_raw)

opener_in = gzip.open if inp.endswith(".gz") else open
with opener_in(inp, "rt", encoding="utf-8", errors="ignore") as fin, gzip.open(outp, "wt") as fout:
    for i, line in enumerate(fin):
        if i >= n:
            break
        fout.write(line)
PY
}

task="${SLURM_ARRAY_TASK_ID:?need SLURM_ARRAY_TASK_ID}"
line_num=$((task + 1))
line="$(sed -n "${line_num}p" "$SAMPLES_TSV" || true)"
[[ -n "$line" ]] || { echo "[ERROR] No line $line_num in $SAMPLES_TSV" >&2; exit 1; }

sample_id="$(printf '%s\n' "$line" | cut -f1)"
r1="$(printf '%s\n' "$line" | cut -f2)"
r2="$(printf '%s\n' "$line" | cut -f3)"
rest="$(printf '%s\n' "$line" | cut -f4- || true)"

[[ -f "$r1" ]] || { echo "[ERROR] Missing R1 for $sample_id: $r1" >&2; exit 1; }
[[ -f "$r2" ]] || { echo "[ERROR] Missing R2 for $sample_id: $r2" >&2; exit 1; }
[[ -s "$POOL1" ]] || { echo "[ERROR] Missing/empty POOL1: $POOL1" >&2; exit 1; }
[[ -s "$POOL2" ]] || { echo "[ERROR] Missing/empty POOL2: $POOL2" >&2; exit 1; }

if [[ -z "$POOL_SIZE" ]]; then
  p1n="$(count_reads "$POOL1")"
  p2n="$(count_reads "$POOL2")"
  [[ "$p1n" == "$p2n" ]] || {
    echo "[ERROR] Pool mate count mismatch for $LABEL: POOL1=$p1n POOL2=$p2n" >&2
    exit 1
  }
  POOL_SIZE="$p1n"
fi

mkdir -p \
  "$WORK/spiked_fastqs" \
  "$WORK/logs" \
  "$WORK/logs/spiked_manifest_parts" \
  "$WORK/tmp"

safe_sample="$(printf '%s' "$sample_id" | tr '/ ' '__')"
tmp="$WORK/tmp/${safe_sample}"
mkdir -p "$tmp"

cleanup() {
  if [[ "$KEEP_TMP" == "1" ]]; then
    echo "[INFO] KEEP_TMP=1, preserving tmp dir: $tmp"
  else
    rm -rf "$tmp"
  fi
}
trap cleanup EXIT

R1n="$(count_reads "$r1")"
R2n="$(count_reads "$r2")"
if (( R1n < R2n )); then
  R="$R1n"
else
  R="$R2n"
fi
[[ "$R" -gt 0 ]] || {
  echo "[ERROR] Zero synchronized pairs for $sample_id (R1=$R1n R2=$R2n)" >&2
  exit 1
}

bg1="$tmp/bg_1.fq.gz"
bg2="$tmp/bg_2.fq.gz"
nlines=$((R * 4))

extract_first_n_lines_gz "$r1" "$nlines" "$bg1"
extract_first_n_lines_gz "$r2" "$nlines" "$bg2"

design="$WORK/logs/${safe_sample}.spike_design.tsv"
manifest_part="$WORK/logs/spiked_manifest_parts/${safe_sample}.tsv"

echo -e "sample_id\tlabel\tfraction\tR1\tR2\tR\tN_requested\tN_inserted\tf_hat\tseed\toutput_r1\toutput_r2\tpool1\tpool2" > "$design"
: > "$manifest_part"

IFS=',' read -r -a fracs <<< "$FRACTIONS"
for f in "${fracs[@]}"; do
  # Stable under manifest reordering, array splitting, and fraction subsetting.
  # The canonical decimal representation prevents 0.10 and 0.1 diverging.
  canonical_fraction="$(python3 - "$f" <<'PY'
import decimal, sys
value = decimal.Decimal(sys.argv[1]).normalize()
print(format(value, "f"))
PY
)"
  seed="$(python3 "$SEED_HELPER" --base "$SEED_BASE" \
    --namespace spike-independent-v1 \
    "$sample_id" "$LABEL" "$canonical_fraction")"

  python3 - <<PY >/dev/null 2>&1
f=float("$f")
assert 0.0 < f < 1.0
PY
  [[ $? -eq 0 ]] || { echo "[ERROR] Invalid fraction: $f" >&2; exit 1; }

  N="$(python3 - <<PY
R=int($R)
f=float("$f")
print(int(round((f/(1.0-f))*R)))
PY
)"
  if (( N <= 0 )); then
    echo "[WARN] $sample_id $LABEL f=$f gives N=0; skipping" >&2
    continue
  fi

  if (( N > POOL_SIZE )); then
    echo "[ERROR] Requested N=$N for $LABEL but pool only has $POOL_SIZE pairs. Increase COVERAGE when building pools." >&2
    exit 1
  fi

  f_tag="$(python3 - <<PY
f=float("$f")
s=("{:.6f}".format(f)).rstrip("0").rstrip(".")
print("f" + s.replace(".", "p"))
PY
)"

  sp1="$tmp/${LABEL}.${f_tag}.spike_1.fq"
  sp2="$tmp/${LABEL}.${f_tag}.spike_2.fq"
  ids1="$tmp/${LABEL}.${f_tag}.ids1.txt"
  ids2="$tmp/${LABEL}.${f_tag}.ids2.txt"

  run_in_container "micromamba run -n ${MAMBA_ENV} seqtk sample -s ${seed} '$POOL1' ${N} > '$sp1'"

  awk 'NR%4==1 {gsub(/^@/,"",$1); print $1}' "$sp1" > "$ids1"

  awk '{
    o=$1; print o;
    id=o;
    if(id ~ /\/1$/){sub(/\/1$/,"/2",id); print id;}
    else if(id ~ /_1$/){sub(/_1$/,"_2",id); print id;}
    else if(id ~ /\.1$/){sub(/\.1$/,".2",id); print id;}
  }' "$ids1" | awk '!seen[$0]++' > "$ids2"

  run_in_container "micromamba run -n ${MAMBA_ENV} seqtk subseq '$POOL2' '$ids2' > '$sp2'"

  n1="$(awk 'END{print int(NR/4)}' "$sp1")"
  n2="$(awk 'END{print int(NR/4)}' "$sp2")"

  if [[ "$n2" != "$n1" ]]; then
    echo "[WARN] $sample_id $LABEL $f_tag mate mismatch after transformed IDs (n1=$n1 n2=$n2). Retrying with ids1." >&2
    run_in_container "micromamba run -n ${MAMBA_ENV} seqtk subseq '$POOL2' '$ids1' > '$sp2'"
    n2="$(awk 'END{print int(NR/4)}' "$sp2")"
  fi

  [[ "$n1" -eq "$N" ]] || {
    echo "[ERROR] seqtk sample returned $n1 reads but N_requested=$N for $LABEL" >&2
    exit 1
  }
  [[ "$n2" -eq "$N" ]] || {
    echo "[ERROR] mate extraction returned $n2 reads but N_requested=$N for $LABEL" >&2
    exit 1
  }

  out1="$WORK/spiked_fastqs/${sample_id}_${LABEL}_${f_tag}_1.fq.gz"
  out2="$WORK/spiked_fastqs/${sample_id}_${LABEL}_${f_tag}_2.fq.gz"

  ( gzip -cd -- "$bg1"; cat "$sp1" ) | gzip -c > "$out1"
  ( gzip -cd -- "$bg2"; cat "$sp2" ) | gzip -c > "$out2"

  f_hat="$(python3 - <<PY
R=int($R)
N=int($n1)
print("{:.8f}".format(N/(R+N)))
PY
)"

  echo -e "${sample_id}\t${LABEL}\t${f}\t${R1n}\t${R2n}\t${R}\t${N}\t${n1}\t${f_hat}\t${seed}\t${out1}\t${out2}\t${POOL1}\t${POOL2}" >> "$design"

  spiked_id="${sample_id}_${LABEL}_${f_tag}"
  if [[ -n "$rest" && "$rest" != "$line" ]]; then
    echo -e "${spiked_id}\t${out1}\t${out2}\t${rest}" >> "$manifest_part"
  else
    echo -e "${spiked_id}\t${out1}\t${out2}" >> "$manifest_part"
  fi
done

echo "[OK] Done sample=$sample_id label=$LABEL"
echo "[OK] Design log: $design"
