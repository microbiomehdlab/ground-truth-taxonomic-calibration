#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
umask 0022
export LC_ALL=C

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

: "${DATASET:?DATASET is required}"
: "${SAMPLE:?SAMPLE is required}"

source "$ROOT/config/global.env"

KEEP_READ_LEVEL="${KEEP_READ_LEVEL:-0}"
KEEP_QC="${KEEP_QC:-0}"
KEEP_OPEN_CLOSED="${KEEP_OPEN_CLOSED:-0}"

SAMPLE="$(printf '%s' "$SAMPLE" | tr -d '\r' | sed 's/[[:space:]]\+$//; s/^[[:space:]]\+//')"
DATASET="$(printf '%s' "$DATASET" | tr -d '\r' | sed 's/[[:space:]]\+$//; s/^[[:space:]]\+//')"

OUT_DIR="$ROOT/results/$DATASET/kraken2_bracken"
mkdir -p "$OUT_DIR"

BRACKEN_S="$OUT_DIR/${SAMPLE}.bracken.S.tsv"
QC="$OUT_DIR/${SAMPLE}.qc.tsv"
OPEN="$OUT_DIR/${SAMPLE}.open_world.csv"
CLOSED="$OUT_DIR/${SAMPLE}.closed_world.csv"
SPEC="$OUT_DIR/${SAMPLE}.species.reads.tsv"

KRAKEN_OUT="$OUT_DIR/${SAMPLE}.kraken2.out"
R2C_ALL="$OUT_DIR/${SAMPLE}.kraken_r2c.tsv"
R2C_SPEC="$OUT_DIR/${SAMPLE}.kraken_r2c.species.tsv"
R2C_GENUS="$OUT_DIR/${SAMPLE}.kraken_r2c.genus.tsv"
PATHS_ALL="$OUT_DIR/${SAMPLE}.taxonomy.paths.all.tsv"
PATHS_SPEC="$OUT_DIR/${SAMPLE}.taxonomy.paths.species.tsv"

# Main final table
python3 - "$BRACKEN_S" "$SPEC" <<'PY'
import sys, csv
inp, outp = sys.argv[1], sys.argv[2]
with open(outp, "w") as fo:
    fo.write("species\treads\n")
    with open(inp, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            if (r.get("taxonomy_lvl","") or "").strip().upper() != "S":
                continue
            name = (r.get("name","") or "").strip()
            try:
                est = float(r.get("new_est_reads","0") or 0)
            except Exception:
                est = 0.0
            if name:
                fo.write(f"{name}\t{int(round(est))}\n")
PY

# Optional QC
if [[ "$KEEP_QC" == "1" ]]; then
  TOTAL=0
  if [[ -s "$OUT_DIR/${SAMPLE}.total_reads.txt" ]]; then
    raw="$(< "$OUT_DIR/${SAMPLE}.total_reads.txt")"
    TOTAL="${raw//[^0-9]/}"
  fi
  [[ "$TOTAL" =~ ^[0-9]+$ ]] || TOTAL=0

  ASSIGNED=$(awk -F'\t' 'NR>1 {s+=$2} END{printf "%d\n", s+0}' "$SPEC")
  UNCLASS=$(( TOTAL - ASSIGNED ))
  (( UNCLASS < 0 )) && UNCLASS=0
  NSPEC=$(awk 'END{print (NR>0?NR-1:0)}' "$SPEC")
  MAPRATE=$(awk -v a="$ASSIGNED" -v t="$TOTAL" 'BEGIN{printf "%.6f",(t>0?a/t:0)}')
  UFRAC=$(awk -v u="$UNCLASS" -v t="$TOTAL" 'BEGIN{printf "%.6f",(t>0?u/t:0)}')

  {
    printf "sample\ttotal_reads\tassigned_reads\tmapping_rate\t#species\tunclassified_reads\tunclassified_frac\n"
    printf "%s\t%d\t%d\t%.6f\t%d\t%d\t%.6f\n" "$SAMPLE" "$TOTAL" "$ASSIGNED" "$MAPRATE" "$NSPEC" "$UNCLASS" "$UFRAC"
  } > "$QC"
fi

# Optional profile exports
if [[ "$KEEP_OPEN_CLOSED" == "1" ]]; then
  awk -F'\t' -v SAMPLE="$SAMPLE" '
    NR==1 { next }
    { name=$1; cnt=$2+0; if (name=="") next; i++; names[i]=name; counts[i]=cnt; total+=cnt }
    END {
      printf "sample"
      for (j=1;j<=i;j++) printf ",%s", names[j]
      printf "\n%s", SAMPLE
      for (j=1;j<=i;j++) printf ",%.6f", (total>0 ? counts[j]*100.0/total : 0)
      printf "\n"
    }
  ' "$SPEC" > "$CLOSED"

  cp -f "$CLOSED" "$OPEN"
fi

# Optional read-level outputs only when explicitly requested
if [[ "$KEEP_READ_LEVEL" == "1" ]]; then
  : "${K2_DB:?K2_DB not set}"
  NODES_DMP="${NODES_DMP:-$K2_DB/taxonomy/nodes.dmp}"
  NAMES_DMP="${NAMES_DMP:-$K2_DB/taxonomy/names.dmp}"
  HELPER="$ROOT/workflows/kraken2_bracken/kraken_r2c.py"
  [[ -s "$HELPER" ]] || { echo "[error] Missing helper: $HELPER" >&2; exit 1; }
  [[ -s "$KRAKEN_OUT" ]] || { echo "[error] Missing $KRAKEN_OUT but KEEP_READ_LEVEL=1" >&2; exit 1; }

  cp -f "$SPEC" "$PATHS_SPEC"
  cp -f "$SPEC" "$PATHS_ALL"

  python3 "$HELPER" \
    --kraken-out "$KRAKEN_OUT" \
    --nodes "$NODES_DMP" \
    --names "$NAMES_DMP" \
    --out-species "$R2C_SPEC" \
    --out-genus "$R2C_GENUS" \
    --out-all "$R2C_ALL"
else
  rm -f "$PATHS_ALL" "$PATHS_SPEC" "$R2C_ALL" "$R2C_SPEC" "$R2C_GENUS" 2>/dev/null || true
fi

# Remove temp QC helper if QC is off
[[ "$KEEP_QC" == "1" ]] || rm -f "$OUT_DIR/${SAMPLE}.total_reads.txt" 2>/dev/null || true

echo "[postprocess done] $SAMPLE"
printf "  %s\n" "$SPEC"