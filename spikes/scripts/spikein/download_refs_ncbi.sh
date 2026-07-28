#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage(){
  cat <<'USAGE'
Usage:
  download_refs_ncbi.sh --panel spike_panel.tsv --outdir refs/

Downloads missing genome FASTAs for spike-in taxa.

Panel TSV header expected (extra columns OK):
  label  taxon_name  assembly  fasta  weight  url

Rules:
  - If fasta exists -> skip
  - Else if url provided -> download url (supports .gz)
  - Else if assembly provided -> use NCBI Datasets CLI:
      datasets download genome accession <assembly> --include genome --filename <zip>
    then extract the genomic FASTA (*_genomic.fna or *_genomic.fna.gz).

Requires:
  - datasets CLI installed and in PATH (for assembly downloads)
  - unzip
  - wget or curl (for url downloads)

Notes:
  - If your HPC blocks outbound internet on compute nodes, run this on a login/transfer node.
USAGE
}

PANEL=""; OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --panel) PANEL="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$PANEL" && -n "$OUTDIR" ]] || { usage; exit 2; }
[[ -f "$PANEL" ]] || { echo "[ERROR] Missing panel: $PANEL" >&2; exit 1; }
mkdir -p "$OUTDIR"

command -v unzip >/dev/null 2>&1 || { echo "[ERROR] unzip not found in PATH." >&2; exit 1; }

have_wget=false; command -v wget >/dev/null 2>&1 && have_wget=true
have_curl=false; command -v curl >/dev/null 2>&1 && have_curl=true
have_datasets=false; command -v datasets >/dev/null 2>&1 && have_datasets=true

download_url(){
  local url="$1"
  local out="$2"
  local tmp="${out}.download"

  if $have_wget; then
    wget -O "$tmp" "$url"
  elif $have_curl; then
    curl -L -o "$tmp" "$url"
  else
    echo "[ERROR] Need wget or curl to download URLs." >&2
    return 1
  fi

  if [[ "$tmp" == *.gz || "$url" == *.gz ]]; then
    gunzip -c "$tmp" > "$out"
    rm -f "$tmp"
  else
    mv "$tmp" "$out"
  fi

  [[ -s "$out" ]] || { echo "[ERROR] URL download produced empty fasta: $out" >&2; return 1; }
}

datasets_download(){
  local asm="$1"
  local zip="$2"

  $have_datasets || { echo "[ERROR] 'datasets' CLI not found in PATH." >&2; return 1; }

  # Try with --no-progressbar (newer CLI), then retry without it (older CLI)
  set +e
  datasets download genome accession "$asm" --include genome --filename "$zip" --no-progressbar >/dev/null 2>&1
  rc=$?
  set -e
  if [[ $rc -eq 0 ]]; then return 0; fi

  datasets download genome accession "$asm" --include genome --filename "$zip"
}

extract_genomic_fna(){
  local unzipped_dir="$1"
  local out="$2"

  local fna
  fna="$(find "$unzipped_dir" -type f \( -name "*_genomic.fna" -o -name "*_genomic.fna.gz" \) | head -n 1 || true)"
  [[ -n "$fna" ]] || { echo "[ERROR] Could not find *_genomic.fna in $unzipped_dir" >&2; return 1; }

  if [[ "$fna" == *.gz ]]; then
    gunzip -c "$fna" > "$out"
  else
    cp "$fna" "$out"
  fi

  [[ -s "$out" ]] || { echo "[ERROR] Extracted fasta is empty: $out" >&2; return 1; }
}

download_assembly(){
  local asm="$1"
  local out="$2"
  local label="$3"

  local tmpdir zip unz
  tmpdir="$(mktemp -d)"
  zip="$tmpdir/${label}_${asm}.zip"
  unz="$tmpdir/unzipped"
  mkdir -p "$unz"

  echo "[DL] $label assembly=$asm"
  datasets_download "$asm" "$zip"

  # Validate zip
  unzip -t "$zip" >/dev/null 2>&1 || { echo "[ERROR] Downloaded zip looks invalid for $asm (maybe blocked/truncated)." >&2; rm -rf "$tmpdir"; return 1; }

  unzip -q "$zip" -d "$unz"
  extract_genomic_fna "$unz" "$out"

  rm -rf "$tmpdir"
}

# columns: label taxon_name assembly fasta weight url ...
tail -n +2 "$PANEL" | while IFS=$'\t' read -r label taxon assembly fasta weight url rest; do
  [[ -z "${label:-}" ]] && continue
  [[ "$label" =~ ^# ]] && continue

  if [[ -z "${fasta:-}" || "$fasta" == "-" ]]; then
    fasta="$OUTDIR/${label}.fa"
  fi
  mkdir -p "$(dirname "$fasta")"

  if [[ -s "$fasta" ]]; then
    echo "[SKIP] $label fasta exists: $fasta"
    continue
  fi

  if [[ -n "${url:-}" ]]; then
    echo "[GET] $label via URL -> $fasta"
    download_url "$url" "$fasta"
    continue
  fi

  if [[ -n "${assembly:-}" ]]; then
    echo "[GET] $label via assembly -> $fasta"
    download_assembly "$assembly" "$fasta" "$label"
    continue
  fi

  echo "[WARN] $label has no fasta, no url, and no assembly accession; cannot fetch automatically." >&2
done
