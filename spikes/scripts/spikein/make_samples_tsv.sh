#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  make_samples_tsv.sh <root_dir> [<root_dir2> ...] [options] > samples.tsv

Creates a paired-end manifest (TSV) on stdout:
  sample_id<TAB>fastq1<TAB>fastq2

This version is dataset-agnostic:
- It does NOT assume SID* folders.
- It finds sample FASTQs by searching for one (or more) "sample subdirs"
  (default: preprocessed) anywhere under each root.
- It infers R1/R2 pairs using common conventions:
    *_1.fastq.gz    <-> *_2.fastq.gz
    *_1.fq.gz       <-> *_2.fq.gz
    *_R1_*.fastq.gz <-> *_R2_*.fastq.gz
    *_R1*.fastq.gz  <-> *_R2*.fastq.gz
    *read1*         <-> *read2*
    *.R1.*          <-> *.R2.*

If multiple pairs are found for a sample folder (e.g., multi-lane), choose behavior with --multi.

Positional args:
  <root_dir> ...   One or more roots to scan. Typically:
                  - data/processed/metaprep   (scans FengQ, ZellerG, etc.)
                  - or a single cohort: .../FengQ/sequencing

Options:
  --subdir NAME[,NAME2,...]   Subdir(s) that contain FASTQs (default: preprocessed)
  --maxdepth N                Max depth to search for subdir (default: 6)
  --fastq-glob GLOB           FASTQ filename glob (default: *.fastq.gz,*.fq.gz)
  --id-depth N                Build sample_id from the last N path components above SUBDIR
                              relative to the root (default: 1)
                              Example (root=metaprep, path=FengQ/sequencing/SID31004/preprocessed):
                                id-depth 1 -> SID31004
                                id-depth 2 -> sequencing__SID31004
                                id-depth 3 -> FengQ__sequencing__SID31004
  --id-sep STR                Separator for id-depth joining (default: __)
  --paired-only true|false    If true, missing pair counts as an error (default: true)
  --strict true|false         If true, exit non-zero on any errors (default: true)
  --multi error|largest|emit  What to do if multiple R1/R2 pairs are found in one folder:
                              - error   : report error (default)
                              - largest : pick largest pair (by total bytes) and warn
                              - emit    : output one row per pair, with unique sample_id suffix

Examples:
  # Scan everything under metaprep (FengQ + ZellerG + others), sample_id includes cohort:
  make_samples_tsv.sh data/processed/metaprep --subdir preprocessed --maxdepth 8 --id-depth 3 > samples.tsv

  # Scan just one cohort:
  make_samples_tsv.sh data/processed/metaprep/ZellerG/sequencing --subdir preprocessed > samples.tsv

  # Multi-lane: emit one row per lane/pair
  make_samples_tsv.sh data/processed/metaprep --multi emit --id-depth 3 > samples.tsv
EOF
}

# ---------------- parse roots ----------------
ROOTS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --*) break ;;
    *) ROOTS+=("$1"); shift ;;
  esac
done

[[ ${#ROOTS[@]} -gt 0 ]] || { usage; exit 2; }

# ---------------- defaults ----------------
SUBDIRS=("preprocessed")
MAXDEPTH=6
FASTQ_GLOB="*.fastq.gz,*.fq.gz"
ID_DEPTH=1
ID_SEP="__"
PAIRED_ONLY="true"
STRICT="true"
MULTI="error"

# ---------------- parse options ----------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --subdir) IFS=',' read -r -a SUBDIRS <<< "$2"; shift 2 ;;
    --maxdepth) MAXDEPTH="$2"; shift 2 ;;
    --fastq-glob) FASTQ_GLOB="$2"; shift 2 ;;
    --id-depth) ID_DEPTH="$2"; shift 2 ;;
    --id-sep) ID_SEP="$2"; shift 2 ;;
    --paired-only) PAIRED_ONLY="$2"; shift 2 ;;
    --strict) STRICT="$2"; shift 2 ;;
    --multi) MULTI="$2"; shift 2 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

# Validate enums
case "$MULTI" in error|largest|emit) ;; *) echo "[ERROR] --multi must be error|largest|emit" >&2; exit 2 ;; esac
case "$PAIRED_ONLY" in true|false) ;; *) echo "[ERROR] --paired-only must be true|false" >&2; exit 2 ;; esac
case "$STRICT" in true|false) ;; *) echo "[ERROR] --strict must be true|false" >&2; exit 2 ;; esac

# Expand FASTQ glob list
IFS=',' read -r -a FASTQ_PATTERNS <<< "$FASTQ_GLOB"

errs=0

echo -e "sample_id\tfastq1\tfastq2"

# --- helpers ---
abs_path() {
  # realpath is standard on Linux; fallback to python if absent
  if command -v realpath >/dev/null 2>&1; then
    realpath "$1"
  else
    python3 - <<PY
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
  fi
}

join_last_components() {
  local root_abs="$1"
  local parent_abs="$2"
  local depth="$3"
  local sep="$4"

  # compute relative path from root to parent
  local rel="$parent_abs"
  if [[ "$parent_abs" == "$root_abs"* ]]; then
    rel="${parent_abs#"$root_abs"/}"
  fi

  # split rel by '/'
  IFS='/' read -r -a parts <<< "$rel"
  local n="${#parts[@]}"
  if (( n == 0 )); then
    echo "$(basename "$parent_abs")"
    return
  fi

  local start=$(( n - depth ))
  (( start < 0 )) && start=0

  local out=""
  for ((i=start; i<n; i++)); do
    if [[ -z "$out" ]]; then out="${parts[i]}"; else out="${out}${sep}${parts[i]}"; fi
  done
  echo "$out"
}

pair_for_r1() {
  # prints matched r2 path on stdout, or empty if none
  local f="$1"
  local bn
  bn="$(basename "$f")"

  # Try safe, specific patterns first
  local candidates=()

  # Illumina style: _R1_001.fastq.gz
  if [[ "$bn" == *"_R1_"*".fastq.gz" ]]; then
    candidates+=("${f/_R1_/_R2_}")
  fi

  # Generic R1 token
  if [[ "$bn" == *"_R1"*".fastq.gz" || "$bn" == *"_R1"*".fq.gz" ]]; then
    candidates+=("${f/_R1/_R2}")
  fi

  # Dot token: .R1.
  if [[ "$bn" == *".R1."* ]]; then
    candidates+=("${f/.R1./.R2.}")
  fi

  # read1/read2
  if [[ "$bn" == *"read1"* ]]; then
    candidates+=("${f/read1/read2}")
  fi

  # Suffix _1_001.fastq.gz
  if [[ "$bn" == *"_1_001.fastq.gz" ]]; then
    candidates+=("${f/_1_001.fastq.gz/_2_001.fastq.gz}")
  fi

  # Suffix _1.fastq.gz / _1.fq.gz
  if [[ "$bn" == *"_1.fastq.gz" ]]; then
    candidates+=("${f/_1.fastq.gz/_2.fastq.gz}")
  fi
  if [[ "$bn" == *"_1.fq.gz" ]]; then
    candidates+=("${f/_1.fq.gz/_2.fq.gz}")
  fi

  # Less strict fallback: replace _1_ with _2_
  if [[ "$bn" == *"_1_"* ]]; then
    candidates+=("${f/_1_/_2_}")
  fi

  # Return the first existing candidate that differs from f
  for c in "${candidates[@]}"; do
    [[ "$c" != "$f" && -f "$c" ]] && { echo "$c"; return 0; }
  done
  echo ""
  return 0
}

# --- main scan ---
for root in "${ROOTS[@]}"; do
  [[ -d "$root" ]] || { echo "[ERROR] Not a directory: $root" >&2; errs=$((errs+1)); continue; }

  root_abs="$(abs_path "$root")"

  # Find candidate FASTQ-holding directories
  scan_dirs=()
  for sub in "${SUBDIRS[@]}"; do
    while read -r d; do
      scan_dirs+=("$d")
    done < <(find "$root" -maxdepth "$MAXDEPTH" -type d -name "$sub" 2>/dev/null | sort)
  done

  # If none found, fall back to using the root itself
  if [[ ${#scan_dirs[@]} -eq 0 ]]; then
    scan_dirs+=("$root")
  fi

  for fqdir in "${scan_dirs[@]}"; do
    # Determine sample folder (= parent of fqdir if fqdir is a SUBDIR match; else fqdir itself)
    is_subdir_match=false
    for sub in "${SUBDIRS[@]}"; do
      [[ "$(basename "$fqdir")" == "$sub" ]] && is_subdir_match=true
    done

    if [[ "$is_subdir_match" == true ]]; then
      sample_parent="$(dirname "$fqdir")"
    else
      sample_parent="$fqdir"
    fi

    parent_abs="$(abs_path "$sample_parent")"
    sid="$(join_last_components "$root_abs" "$parent_abs" "$ID_DEPTH" "$ID_SEP")"

    # Collect FASTQs in this directory (not recursive at this level)
    files=()
    for pat in "${FASTQ_PATTERNS[@]}"; do
      while read -r f; do files+=("$f"); done < <(find "$fqdir" -maxdepth 1 -type f -name "$pat" 2>/dev/null | sort)
    done

    # Deduplicate
    if [[ ${#files[@]} -gt 1 ]]; then
      mapfile -t files < <(printf "%s\n" "${files[@]}" | awk '!seen[$0]++')
    fi

    # Try to infer pairs
    declare -A seen_r1=()
    pairs_r1=()
    pairs_r2=()

    for f in "${files[@]}"; do
      [[ -f "$f" ]] || continue
      bn="$(basename "$f")"
      # Heuristic: only attempt pairing from likely-R1 names
      if [[ "$bn" == *"_R1"* || "$bn" == *"_1"* || "$bn" == *"read1"* || "$bn" == *".R1."* ]]; then
        r2="$(pair_for_r1 "$f")"
        if [[ -n "$r2" ]]; then
          # avoid duplicates
          if [[ -z "${seen_r1[$f]+x}" ]]; then
            seen_r1["$f"]=1
            pairs_r1+=("$f")
            pairs_r2+=("$r2")
          fi
        fi
      fi
    done

    if [[ ${#pairs_r1[@]} -eq 0 ]]; then
      # If no pairs found, but directory had fastqs, warn once.
      if [[ ${#files[@]} -gt 0 ]]; then
        echo "[WARN] No paired-end FASTQs inferred for sample '$sid' in: $fqdir" >&2
      fi
      [[ "$PAIRED_ONLY" == "true" ]] && errs=$((errs+1))
      continue
    fi

    if [[ ${#pairs_r1[@]} -eq 1 ]]; then
      echo -e "${sid}\t${pairs_r1[0]}\t${pairs_r2[0]}"
      continue
    fi

    # Multiple pairs found
    case "$MULTI" in
      error)
        echo "[WARN] Multiple R1/R2 pairs found for sample '$sid' in: $fqdir" >&2
        for i in "${!pairs_r1[@]}"; do
          echo "[WARN]   ${pairs_r1[$i]}  ||  ${pairs_r2[$i]}" >&2
        done
        errs=$((errs+1))
        ;;
      largest)
        best_i=0
        best_sz=-1
        for i in "${!pairs_r1[@]}"; do
          sz1=$(stat -c%s "${pairs_r1[$i]}" 2>/dev/null || echo 0)
          sz2=$(stat -c%s "${pairs_r2[$i]}" 2>/dev/null || echo 0)
          tot=$((sz1+sz2))
          if (( tot > best_sz )); then best_sz=$tot; best_i=$i; fi
        done
        echo "[WARN] Multiple pairs for '$sid' in $fqdir; picking largest (bytes=$best_sz)." >&2
        echo -e "${sid}\t${pairs_r1[$best_i]}\t${pairs_r2[$best_i]}"
        ;;
      emit)
        echo "[WARN] Multiple pairs for '$sid' in $fqdir; emitting one row per pair." >&2
        for i in "${!pairs_r1[@]}"; do
          stem="$(basename "${pairs_r1[$i]}")"
          stem="${stem%.fastq.gz}"; stem="${stem%.fq.gz}"
          echo -e "${sid}${ID_SEP}${stem}\t${pairs_r1[$i]}\t${pairs_r2[$i]}"
        done
        ;;
    esac
  done
done

if [[ "$STRICT" == "true" && "$errs" -gt 0 ]]; then
  echo "[ERROR] Encountered $errs errors while building samples.tsv" >&2
  exit 1
fi