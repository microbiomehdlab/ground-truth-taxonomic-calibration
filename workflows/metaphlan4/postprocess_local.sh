#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

: "${DATASET:?}"; : "${SAMPLE:?}"

# Normalize SAMPLE/DATASET to avoid stray spaces/CRLF
SAMPLE="$(printf '%s' "$SAMPLE" | tr -d '\r' | sed 's/[[:space:]]\+$//; s/^[[:space:]]\+//')"
DATASET="$(printf '%s' "$DATASET" | tr -d '\r' | sed 's/[[:space:]]\+$//; s/^[[:space:]]\+//')"


OUT_DIR="$ROOT/results/$DATASET/metaphlan4"
RAW="$OUT_DIR/${SAMPLE}.metaphlan.tsv"
TOTAL_FILE="$OUT_DIR/${SAMPLE}.total_reads.txt"

# ---- TOTAL reads (plain integer) ----
if [[ -s "$TOTAL_FILE" ]]; then
  TOTAL_RAW="$(< "$TOTAL_FILE")"
  TOTAL="${TOTAL_RAW//[^0-9]/}"
  : "${TOTAL:=0}"
else
  TOTAL=0
fi
[[ "$TOTAL" =~ ^[0-9]+$ ]] || { echo "[error] Could not parse TOTAL reads"; exit 1; }

SPEC="$OUT_DIR/${SAMPLE}.species.reads.tsv"
QC="$OUT_DIR/${SAMPLE}.qc.tsv"
OPEN="$OUT_DIR/${SAMPLE}.open_world.split.csv"
CLOSED="$OUT_DIR/${SAMPLE}.closed_world.csv"

# ---- Extract species-only (terminal s__... WITHOUT any |t__) ----
# MetaPhlAn columns: clade_name \t ncbi_tax_id \t relative_abundance \t additional_species
# Keep only rows where clade_name ends with s__Species (no trailing |t__...).
tmp_species="$OUT_DIR/${SAMPLE}.species.tmp.tsv"
awk -F'\t' '
  BEGIN{OFS="\t"}
  # species-only rows: clade_name ends with s__... and has no further "|"
  $1 ~ /(^|[|])s__[^|]+$/ {
    sp=$1
    sub(/^.*\|s__/, "", sp)     # drop path, keep species token after last "|s__"
    sub(/^s__/, "", sp)         # also handle case where row starts at species level
    gsub(/_/," ", sp)           # Escherichia_coli -> Escherichia coli
    ra=$3+0                     # column 3 = relative abundance (percent)
    if (ra>0) print sp, ra
  }
' "$RAW" | sort -k1,1 > "$tmp_species" || :  # empty ok

# ---- species.reads.tsv (reads = round(TOTAL * pct/100)) ----
echo -e "species\treads" > "$SPEC"
if [[ -s "$tmp_species" ]]; then
  awk -F'\t' -v T="$TOTAL" '{reads=int(T*($2/100.0)+0.5); print $1 "\t" reads}' "$tmp_species" >> "$SPEC"
fi

# ---- QC and fractions ----
SUM_PCT=$(awk -F'\t' '{s+=$2}END{printf "%.6f", (s+0)}' "$tmp_species")
# Assigned reads from TOTAL * SUM_PCT/100 (avoid rounding drift)
ASSIGNED=$(awk -v t="$TOTAL" -v p="$SUM_PCT" 'BEGIN{printf "%d", int(t*(p/100.0)+0.5)}')
NSPEC=$(awk 'END{print NR-1}' "$SPEC")
UNCLASS_PCT=$(awk -v s="$SUM_PCT" 'BEGIN{u=100.0-s; if(u<0)u=0; if(u>100)u=100; printf "%.6f", u}')
UNCLASS=$(( TOTAL - ASSIGNED )); (( UNCLASS < 0 )) && UNCLASS=0
MAPRATE=$(awk -v a="$ASSIGNED" -v t="$TOTAL" 'BEGIN{printf "%.6f", (t>0? a/t:0)}')
UFRAC=$(awk -v u="$UNCLASS" -v t="$TOTAL" 'BEGIN{printf "%.6f", (t>0? u/t:0)}')

printf "sample\ttotal_reads\tassigned_reads\tmapping_rate\t#species\tunclassified_reads\tunclassified_frac\n" > "$QC"
printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\n" "$SAMPLE" "$TOTAL" "$ASSIGNED" "$MAPRATE" "$NSPEC" "$UNCLASS" "$UFRAC" >> "$QC"

# ---- Open-world CSV: species % (from MetaPhlAn) + AboveSpecies(=0) + Unclassified% ----
{
  printf "sample"
  if [[ -s "$tmp_species" ]]; then awk -F'\t' '{printf ",%s",$1}' "$tmp_species"; fi
  printf ",AboveSpecies,Unclassified\n"

  printf "%s" "$SAMPLE"
  if [[ -s "$tmp_species" ]]; then awk -F'\t' '{printf ",%.6f",$2}' "$tmp_species"; fi
  printf ",%.6f,%.6f\n" 0.0 "$UNCLASS_PCT"
} > "$OPEN"

# ---- Closed-world CSV: species renormalized to 100% ----
{
  printf "sample"
  if [[ -s "$tmp_species" ]]; then awk -F'\t' '{printf ",%s",$1}' "$tmp_species"; fi
  printf "\n"

  printf "%s" "$SAMPLE"
  if [[ -s "$tmp_species" ]]; then
    awk -F'\t' -v s="$SUM_PCT" '{val=(s>0? ($2*100.0/s):0); printf ",%.6f", val}' "$tmp_species"
  fi
  printf "\n"
} > "$CLOSED"

rm -f "$tmp_species"
echo "[postprocess done] $SAMPLE"
