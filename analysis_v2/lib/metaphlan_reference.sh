#!/usr/bin/env bash
# Shared resolution of the profiler-specific MetaPhlAn reference inputs.
#
# Bracken keeps the read-proportional reference e = (1-F)o + f unchanged.
# MetaPhlAn uses the decided genome-equivalent reference, which needs the exact
# implanted assembly lengths (G_t) and the audited sample-wide effective
# community genome sizes (G_eff,i).
#
# Sourced by the definitive runners. There is no default: a canonical input
# containing MetaPhlAn rows without both reference tables is a hard failure,
# because silently defaulting is how the read-fraction / genome-equivalent
# estimand mismatch entered the analysis in the first place.
#
# Required environment when MetaPhlAn rows are present:
#   TARGET_GENOME_SIZES     TSV: target_label, genome_size_bp
#   EFFECTIVE_GENOME_SIZES  TSV: cohort, sample_id, effective_genome_size_bp
# Optional:
#   METAPHLAN_READ_REFERENCE_SENSITIVITY=0  disables the labelled sensitivity arm

# A canonical table with no profiler column is malformed. Treating it as
# Bracken-only would silently skip the MetaPhlAn reference requirement.
canonical_has_metaphlan() {
  local canonical="$1"
  if ! awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "profiler") ok = 1 }
                    NR == 1 { exit !ok }' "$canonical"; then
    echo "[ERROR] canonical input has no 'profiler' column: $canonical" >&2
    exit 1
  fi
  awk -F '\t' '
    NR == 1 { for (i = 1; i <= NF; i++) if ($i == "profiler") col = i; next }
    col && $col == "metaphlan4" { found = 1; exit }
    END { exit !found }
  ' "$canonical"
}

# Header check by column name; positional indexing breaks silently if the
# schema is reordered.
require_columns() {
  local path="$1"; shift
  local missing
  missing="$(awk -F '\t' -v want="$*" '
    NR == 1 { for (i = 1; i <= NF; i++) seen[$i] = 1
              n = split(want, cols, " ")
              for (j = 1; j <= n; j++) if (!(cols[j] in seen)) printf "%s ", cols[j]
              exit }' "$path")"
  if [[ -n "$missing" ]]; then
    echo "[ERROR] $path is missing column(s): $missing" >&2
    exit 1
  fi
}

# Every canonical cohort/sample carrying MetaPhlAn rows must have a G_eff entry.
require_geff_coverage() {
  local canonical="$1" geff="$2"
  local missing
  missing="$(awk -F '\t' '
    NR == FNR {
      if (FNR == 1) { for (i = 1; i <= NF; i++) g[$i] = i; next }
      have[$g["cohort"] SUBSEP $g["sample_id"]] = 1; next
    }
    FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $c["profiler"] == "metaphlan4" {
      key = $c["cohort"] SUBSEP $c["sample_id"]
      if (!(key in have)) need[key] = 1
    }
    END { n = 0; for (k in need) n++; print n }
  ' "$geff" "$canonical")"
  if [[ "${missing:-0}" != "0" ]]; then
    echo "[ERROR] $missing MetaPhlAn cohort/sample(s) absent from the G_eff table: $geff" >&2
    exit 1
  fi
}

require_reference_inputs() {
  local canonical="$1"
  if ! canonical_has_metaphlan "$canonical"; then
    echo "[INFO] no MetaPhlAn rows in $canonical; Bracken-only run" >&2
    return 1
  fi
  : "${TARGET_GENOME_SIZES:?MetaPhlAn rows present: set TARGET_GENOME_SIZES to the exact implanted assembly length table}"
  : "${EFFECTIVE_GENOME_SIZES:?MetaPhlAn rows present: set EFFECTIVE_GENOME_SIZES to the audited effective_genome_size.tsv}"
  local path
  for path in "$TARGET_GENOME_SIZES" "$EFFECTIVE_GENOME_SIZES"; do
    if [[ ! -s "$path" ]]; then
      echo "[ERROR] reference table missing or empty: $path" >&2
      exit 1
    fi
  done
  # Schema, not merely non-emptiness.
  require_columns "$TARGET_GENOME_SIZES" target_label genome_size_bp
  require_columns "$EFFECTIVE_GENOME_SIZES" cohort sample_id effective_genome_size_bp
  require_geff_coverage "$canonical" "$EFFECTIVE_GENOME_SIZES"
  return 0
}

# derive_endpoints_with_references <canonical> <endpoint_root>
#
# Writes the primary genome-equivalent result to <endpoint_root> and, unless
# disabled, the labelled read-proportional sensitivity to
# <endpoint_root>_read_reference_sensitivity. The two never share a directory,
# so a sensitivity result cannot be mistaken for, or overwrite, the primary.
derive_endpoints_with_references() {
  local canonical="$1" endpoint_root="$2"
  local sensitivity_root="${endpoint_root}_read_reference_sensitivity"

  if ! require_reference_inputs "$canonical"; then
    python3 analysis_v2/scripts/derive_paired_endpoints.py \
      --input "$canonical" --outdir "$endpoint_root"
    return
  fi

  if [[ -e "$endpoint_root" || -e "$sensitivity_root" ]]; then
    echo "[ERROR] endpoint output already exists; use a new run root: $endpoint_root" >&2
    exit 1
  fi

  local target_resolved effective_resolved
  target_resolved="$(readlink -f "$TARGET_GENOME_SIZES")"
  effective_resolved="$(readlink -f "$EFFECTIVE_GENOME_SIZES")"

  python3 analysis_v2/scripts/derive_paired_endpoints.py \
    --input "$canonical" --outdir "$endpoint_root" \
    --metaphlan-reference genome_equivalent \
    --target-genome-sizes "$target_resolved" \
    --effective-genome-sizes "$effective_resolved"

  if [[ ! -s "$endpoint_root/SUCCESS" ]]; then
    echo "[ERROR] primary endpoint stage produced no SUCCESS: $endpoint_root" >&2
    exit 1
  fi

  {
    printf 'field\tvalue\n'
    printf 'metaphlan_reference\tgenome_equivalent\n'
    printf 'role\tPRIMARY\n'
    printf 'target_genome_sizes\t%s\n' "$target_resolved"
    printf 'effective_genome_sizes\t%s\n' "$effective_resolved"
    printf 'target_genome_sizes_sha256\t%s\n' "$(sha256sum "$target_resolved" | cut -d' ' -f1)"
    printf 'effective_genome_sizes_sha256\t%s\n' "$(sha256sum "$effective_resolved" | cut -d' ' -f1)"
    printf 'canonical_input\t%s\n' "$(readlink -f "$canonical")"
    printf 'canonical_input_sha256\t%s\n' "$(sha256sum "$canonical" | cut -d' ' -f1)"
    printf 'command\t%s\n' "derive_paired_endpoints.py --metaphlan-reference genome_equivalent"
    printf 'generated_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "$endpoint_root/metaphlan_reference_provenance.tsv"

  if [[ "${METAPHLAN_READ_REFERENCE_SENSITIVITY:-1}" == "1" ]]; then
    python3 analysis_v2/scripts/derive_paired_endpoints.py \
      --input "$canonical" --outdir "$sensitivity_root" \
      --metaphlan-reference read_proportional
    {
      printf 'field\tvalue\n'
      printf 'metaphlan_reference\tread_proportional\n'
      printf 'role\tSENSITIVITY\n'
      printf 'primary_endpoints\t%s\n' "$(readlink -f "$endpoint_root")"
      printf 'canonical_input\t%s\n' "$(readlink -f "$canonical")"
      printf 'generated_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "$sensitivity_root/metaphlan_reference_provenance.tsv"
    printf 'SENSITIVITY_ONLY\tmetaphlan_read_proportional\n' \
      > "$sensitivity_root/SENSITIVITY_ONLY.txt"
    if [[ ! -s "$sensitivity_root/SUCCESS" ]]; then
      echo "[ERROR] sensitivity stage produced no SUCCESS: $sensitivity_root" >&2
      exit 1
    fi
    echo "[INFO] read-proportional MetaPhlAn sensitivity: $sensitivity_root" >&2
  else
    printf 'sensitivity_disabled\tMETAPHLAN_READ_REFERENCE_SENSITIVITY=0\n' \
      >> "$endpoint_root/metaphlan_reference_provenance.tsv"
    echo "[INFO] read-proportional sensitivity disabled by request" >&2
  fi
}
