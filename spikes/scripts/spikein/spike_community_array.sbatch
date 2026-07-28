#!/usr/bin/env bash
#SBATCH -J spike_community
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -t 32:00:00

set -euo pipefail
IFS=$'\n\t'

: "${IMG:?set IMG}"
: "${WORK:?set WORK (e.g. \$WORK/community/CRCpanel)}"
: "${SAMPLES_TSV:?set SAMPLES_TSV}"
: "${COMMUNITY_TSV:?set COMMUNITY_TSV (label\\tweight)}"
: "${POOLS_DIR:?set POOLS_DIR}"
: "${COMMUNITY_LABEL:?set COMMUNITY_LABEL}"
: "${SEED_BASE:?set SEED_BASE}"

# Backward-compatible:
# - FULL_FRACTIONS preserves the original full list for seed indexing
# - RUN_FRACTIONS is the subset to execute in this submission
FULL_FRACTIONS="${FULL_FRACTIONS:-${FRACTIONS:?set FULL_FRACTIONS or FRACTIONS}}"
RUN_FRACTIONS="${RUN_FRACTIONS:-$FULL_FRACTIONS}"

MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"
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
  # Avoid micromamba/libmamba global lock contention across array tasks:
  # give each task its own cache directory.
  local cache="$WORK/tmp/xdg_cache_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  mkdir -p "$cache"
  local bind_args=()
  [[ -n "$BIND" ]] && bind_args=(-B "$BIND")

  "$CTR" exec \
    "${bind_args[@]}" \
    -B "$WORK:$WORK" \
    --env XDG_CACHE_HOME="$cache" \
    "$IMG" bash -lc "$1"
}

count_reads() {
  local fq="$1"
  if [[ "$fq" == *.gz ]]; then
    zcat "$fq" | awk 'END{print int(NR/4)}'
  else
    awk 'END{print int(NR/4)}' "$fq"
  fi
}

frac_tag() {
  python3 - "$1" <<'PY'
import sys
f=float(sys.argv[1])
s=("{:.6f}".format(f)).rstrip("0").rstrip(".")
print("f"+s.replace(".","p"))
PY
}

task="${SLURM_ARRAY_TASK_ID:?need array}"
line_num=$((task + 1))
line="$(sed -n "${line_num}p" "$SAMPLES_TSV" || true)"
[[ -n "$line" ]] || { echo "[ERROR] No line $line_num in $SAMPLES_TSV" >&2; exit 1; }

sample_id="$(echo "$line" | cut -f1)"
r1="$(echo "$line" | cut -f2)"
r2="$(echo "$line" | cut -f3)"
rest="$(echo "$line" | cut -f4- || true)"

[[ -f "$r1" ]] || { echo "[ERROR] Missing R1 for $sample_id: $r1" >&2; exit 1; }
[[ -f "$r2" ]] || { echo "[ERROR] Missing R2 for $sample_id: $r2" >&2; exit 1; }

mkdir -p "$WORK/spiked_fastqs" "$WORK/logs" "$WORK/logs/spiked_manifest_parts" "$WORK/tmp"

safe_sample="$(echo "$sample_id" | tr '/ ' '__')"
tmp="$WORK/tmp/${safe_sample}"
mkdir -p "$tmp"

IFS=',' read -r -a full_fracs <<< "$FULL_FRACTIONS"
IFS=',' read -r -a run_fracs <<< "$RUN_FRACTIONS"
(( ${#run_fracs[@]} > 0 )) || { echo "[ERROR] RUN_FRACTIONS is empty" >&2; exit 1; }

orig_frac_idx() {
  local target="$1"
  local i
  for i in "${!full_fracs[@]}"; do
    [[ "${full_fracs[$i]}" == "$target" ]] && { echo "$i"; return 0; }
  done
  echo "[ERROR] Fraction '$target' not found in FULL_FRACTIONS='$FULL_FRACTIONS'" >&2
  exit 1
}

# read community weights
labels=()
weights=()
sumw=0.0
while IFS=$'\t' read -r lab w; do
  [[ -z "${lab:-}" ]] && continue
  [[ "$lab" =~ ^# ]] && continue
  [[ "$lab" == "label" ]] && continue
  w="${w:-1}"
  labels+=("$lab")
  weights+=("$w")
  sumw="$(python3 - <<PY
s=float("$sumw"); s+=float("$w"); print(s)
PY
)"
done < "$COMMUNITY_TSV"
(( ${#labels[@]} > 0 )) || { echo "[ERROR] COMMUNITY_TSV has no taxa" >&2; exit 1; }

R1n="$(count_reads "$r1")"
R2n="$(count_reads "$r2")"
if (( R1n < R2n )); then R="$R1n"; else R="$R2n"; fi
[[ "$R" -gt 0 ]] || { echo "[ERROR] Zero pairs for $sample_id" >&2; exit 1; }
nlines=$((R*4))

bg1="$tmp/bg_1.fq.gz"
bg2="$tmp/bg_2.fq.gz"
if [[ "$r1" == *.gz ]]; then zcat "$r1" | head -n "$nlines" | gzip -c > "$bg1"; else head -n "$nlines" "$r1" | gzip -c > "$bg1"; fi
if [[ "$r2" == *.gz ]]; then zcat "$r2" | head -n "$nlines" | gzip -c > "$bg2"; else head -n "$nlines" "$r2" | gzip -c > "$bg2"; fi

design="$WORK/logs/${safe_sample}.${COMMUNITY_LABEL}.design.tsv"
manifest_part="$WORK/logs/spiked_manifest_parts/${safe_sample}.tsv"

if [[ ! -f "$design" ]]; then
  echo -e "sample_id\tcommunity\tfraction\tR1\tR2\tR\tN_total\tf_hat\tseed_base" > "$design"
fi
touch "$manifest_part"

for f in "${run_fracs[@]}"; do
  idx="$(orig_frac_idx "$f")"
  seed_base=$((SEED_BASE + task*1000 + idx*10))

  N_total="$(python3 - <<PY
import math
R=int($R)
f=float("$f")
print(int(round((f/(1.0-f))*R)))
PY
)"
  if (( N_total <= 0 )); then
    echo "[WARN] $sample_id $COMMUNITY_LABEL f=$f gives N_total=0; skipping" >&2
    continue
  fi

  f_tag="$(frac_tag "$f")"

  sp1_all="$tmp/${COMMUNITY_LABEL}.${f_tag}.spike_1.fq"
  sp2_all="$tmp/${COMMUNITY_LABEL}.${f_tag}.spike_2.fq"
  : > "$sp1_all"; : > "$sp2_all"

  allocated=0
  for i in "${!labels[@]}"; do
    lab="${labels[$i]}"
    w="${weights[$i]}"

    if (( i == ${#labels[@]} - 1 )); then
      Ni=$((N_total - allocated))
    else
      Ni="$(python3 - <<PY
import math
N=int($N_total)
w=float("$w"); s=float("$sumw")
print(int(round(N*(w/s))))
PY
)"
      (( Ni < 0 )) && Ni=0
      allocated=$((allocated + Ni))
    fi
    (( Ni <= 0 )) && continue

    p1="$POOLS_DIR/${lab}.pool_1.fq"
    p2="$POOLS_DIR/${lab}.pool_2.fq"
    [[ -s "$p1" && -s "$p2" ]] || { echo "[ERROR] Missing pool for $lab: $p1 / $p2" >&2; exit 1; }

    seed=$((seed_base + i))
    sp1="$tmp/${lab}.${f_tag}.sp1.fq"
    sp2="$tmp/${lab}.${f_tag}.sp2.fq"
    ids1="$tmp/${lab}.${f_tag}.ids1.txt"
    ids2="$tmp/${lab}.${f_tag}.ids2.txt"

    run_in_container "micromamba run -n ${MAMBA_ENV} seqtk sample -s ${seed} '$p1' ${Ni} > '$sp1'"
    awk 'NR%4==1 {gsub(/^@/,"",$1); print $1}' "$sp1" > "$ids1"
    awk '{
      o=$1; print o;
      id=o;
      if(id ~ /\/1$/){sub(/\/1$/,"/2",id); print id;}
      else if(id ~ /_1$/){sub(/_1$/,"_2",id); print id;}
      else if(id ~ /\.1$/){sub(/\.1$/,".2",id); print id;}
    }' "$ids1" | awk '!seen[$0]++' > "$ids2"
    run_in_container "micromamba run -n ${MAMBA_ENV} seqtk subseq '$p2' '$ids2' > '$sp2'"

    n1="$(awk 'END{print int(NR/4)}' "$sp1")"
    n2="$(awk 'END{print int(NR/4)}' "$sp2")"
    [[ "$n1" == "$n2" ]] || { echo "[ERROR] $sample_id $lab $f_tag mate mismatch (n1=$n1 n2=$n2)" >&2; exit 1; }

    cat "$sp1" >> "$sp1_all"
    cat "$sp2" >> "$sp2_all"
  done

  n1_all="$(awk 'END{print int(NR/4)}' "$sp1_all")"
  n2_all="$(awk 'END{print int(NR/4)}' "$sp2_all")"
  [[ "$n1_all" == "$n2_all" ]] || { echo "[ERROR] Community mate mismatch total (n1=$n1_all n2=$n2_all)" >&2; exit 1; }

  out1="$WORK/spiked_fastqs/${sample_id}_${COMMUNITY_LABEL}_${f_tag}_1.fq.gz"
  out2="$WORK/spiked_fastqs/${sample_id}_${COMMUNITY_LABEL}_${f_tag}_2.fq.gz"
  ( zcat "$bg1"; cat "$sp1_all" ) | gzip -c > "$out1"
  ( zcat "$bg2"; cat "$sp2_all" ) | gzip -c > "$out2"

  f_hat="$(python3 - <<PY
R=int($R); N=int($n1_all)
print("{:.8f}".format(N/(R+N)))
PY
)"
  spiked_id="${sample_id}_${COMMUNITY_LABEL}_${f_tag}"

  # Remove any existing entry for this exact fraction, then append the regenerated one.
  awk -F'\t' -v c="$COMMUNITY_LABEL" -v f="$f" \
    'NR==1 || !($2==c && $3==f)' \
    "$design" > "${design}.tmp" && mv "${design}.tmp" "$design"

  awk -F'\t' -v id="$spiked_id" \
    '$1!=id' \
    "$manifest_part" > "${manifest_part}.tmp" && mv "${manifest_part}.tmp" "$manifest_part"

  echo -e "${sample_id}\t${COMMUNITY_LABEL}\t${f}\t${R1n}\t${R2n}\t${R}\t${n1_all}\t${f_hat}\t${seed_base}" >> "$design"

  if [[ -n "$rest" && "$rest" != "$line" ]]; then
    echo -e "${spiked_id}\t${out1}\t${out2}\t${rest}" >> "$manifest_part"
  else
    echo -e "${spiked_id}\t${out1}\t${out2}" >> "$manifest_part"
  fi
done

echo "[OK] Done $sample_id community=$COMMUNITY_LABEL run_fractions=$RUN_FRACTIONS"
