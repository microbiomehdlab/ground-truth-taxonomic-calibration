#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage(){
  cat <<'USAGE'
Usage:
  download_refs_ncbi.sh --panel spike_panel.tsv --outdir refs/ --image taxonomic-tools.sif

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
  - Apptainer image containing the pinned NCBI Datasets CLI (recommended), or
    datasets installed on the host
  - unzip
  - wget or curl (for url downloads)

Notes:
  - If your HPC blocks outbound internet on compute nodes, run this on a login/transfer node.
USAGE
}

PANEL=""; OUTDIR=""; IMAGE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --panel) PANEL="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
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
DATASETS_CMD=()
if [[ -n "$IMAGE" ]]; then
  [[ -s "$IMAGE" ]] || { echo "[ERROR] Missing image: $IMAGE" >&2; exit 1; }
  command -v apptainer >/dev/null 2>&1 || { echo "[ERROR] apptainer not found" >&2; exit 1; }
  DATASETS_CMD=(apptainer exec --cleanenv "$IMAGE" datasets)
elif command -v datasets >/dev/null 2>&1; then
  DATASETS_CMD=(datasets)
fi
(( ${#DATASETS_CMD[@]} > 0 )) || {
  echo "[ERROR] Provide --image containing datasets, or install datasets on the host." >&2
  exit 1
}
mkdir -p "$OUTDIR/provenance"
"${DATASETS_CMD[@]}" version > "$OUTDIR/provenance/ncbi_datasets_version.txt"

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

  # Try with --no-progressbar (newer CLI), then retry without it (older CLI)
  set +e
  "${DATASETS_CMD[@]}" download genome accession "$asm" --include genome --filename "$zip" --no-progressbar >/dev/null 2>&1
  rc=$?
  set -e
  if [[ $rc -eq 0 ]]; then return 0; fi

  "${DATASETS_CMD[@]}" download genome accession "$asm" --include genome --filename "$zip"
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

  report="$(find "$unz" -type f -name 'assembly_data_report.jsonl' | head -n 1 || true)"
  catalog="$(find "$unz" -type f -name 'dataset_catalog.json' | head -n 1 || true)"
  [[ -n "$report" ]] && cp "$report" "$OUTDIR/provenance/${label}.${asm}.assembly_data_report.jsonl"
  [[ -n "$catalog" ]] && cp "$catalog" "$OUTDIR/provenance/${label}.${asm}.dataset_catalog.json"

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

python3 - "$PANEL" "$OUTDIR/reference_genome_checksums.tsv" <<'PY'
import csv
import hashlib
import pathlib
import sys

panel = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])

def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()

def normalized_content(path):
    records = []
    current = []
    total = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequence = "".join(current).upper()
                    records.append(hashlib.sha256(sequence.encode()).hexdigest())
                    total += len(sequence)
                    current = []
            else:
                current.append("".join(line.split()))
    if current:
        sequence = "".join(current).upper()
        records.append(hashlib.sha256(sequence.encode()).hexdigest())
        total += len(sequence)
    digest = hashlib.sha256("\n".join(sorted(records)).encode()).hexdigest()
    return digest, len(records), total

with panel.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
with output.open("w", encoding="utf-8", newline="") as handle:
    fields = ["label", "assembly", "fasta", "bytes", "sha256", "normalized_sequence_sha256", "contigs", "bases"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        path = pathlib.Path(row["fasta"])
        if not path.is_file() or not path.stat().st_size:
            raise SystemExit(f"[ERROR] Missing downloaded FASTA: {path}")
        normalized, contigs, bases = normalized_content(path)
        writer.writerow({
            "label": row["label"], "assembly": row["assembly"], "fasta": path,
            "bytes": path.stat().st_size, "sha256": file_hash(path),
            "normalized_sequence_sha256": normalized, "contigs": contigs, "bases": bases,
        })
print(f"[OK] Wrote reference provenance: {output}")
PY
