#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}"
[[ -n "${YACHIDA_ENV:-}" ]] && source "$YACHIDA_ENV"
: "${SAMPLE_ID:?provided by stream_sample.py}"
: "${RAW_R1:?provided by stream_sample.py}"
: "${RAW_R2:?provided by stream_sample.py}"
: "${SAMPLE_WORK:?provided by stream_sample.py}"
: "${RECEIPT:?provided by stream_sample.py}"
: "${UPSTREAM_SIF:?set in YACHIDA_ENV}"
: "${HOST_INDEX:?set in YACHIDA_ENV}"

threads="${PREPROCESS_THREADS:-8}"
audit="$SAMPLE_WORK/host_policy_audit"
results="${HOST_POLICY_AUDIT_ROOT:-$ROOT/work/private_host_policy_audit/results}/$SAMPLE_ID"
mkdir -p "$audit" "$results"

fastp_r1="$audit/post_fastp_1.fastq.gz"
fastp_r2="$audit/post_fastp_2.fastq.gz"
apptainer exec --cleanenv "$UPSTREAM_SIF" fastp \
  -i "$RAW_R1" -I "$RAW_R2" -o "$fastp_r1" -O "$fastp_r2" \
  --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 60 --trim_poly_g --low_complexity_filter \
  --complexity_threshold 30 --thread "$threads" --detect_adapter_for_pe \
  --html "$results/fastp.html" --json "$results/fastp.json" \
  >"$results/fastp.log" 2>&1

legacy_r1="$audit/un_conc_1.fastq.gz"
legacy_r2="$audit/un_conc_2.fastq.gz"
apptainer exec --cleanenv "$UPSTREAM_SIF" bowtie2 \
  -x "$HOST_INDEX" -1 "$fastp_r1" -2 "$fastp_r2" --very-sensitive \
  -p "$threads" --un-conc-gz "$audit/un_conc_%.fastq.gz" -S /dev/null \
  >"$results/bowtie2_un_conc.log" 2>&1

strict_r1="$audit/both_unmapped_1.fastq.gz"
strict_r2="$audit/both_unmapped_2.fastq.gz"
apptainer exec --cleanenv "$UPSTREAM_SIF" bowtie2 \
  -x "$HOST_INDEX" -1 "$fastp_r1" -2 "$fastp_r2" --very-sensitive \
  --reorder -p "$threads" \
  2>"$results/bowtie2_both_unmapped.log" |
  python3 "$ROOT/datasets/yachida/extract_strict_unmapped_pairs.py" \
    --r1 "$strict_r1" --r2 "$strict_r2" \
    2>>"$results/bowtie2_both_unmapped.log"

for policy in un_conc both_unmapped; do
  python3 "$ROOT/scripts/validate_paired_fastq.py" \
    --r1 "$audit/${policy}_1.fastq.gz" --r2 "$audit/${policy}_2.fastq.gz" \
    --minimum-pairs 1 --output "$results/${policy}.paired_integrity.tsv"
done

python3 - "$SAMPLE_ID" "$audit" "$results/host_policy_summary.tsv" <<'PY'
import gzip, pathlib, sys
sample, root, output = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
def pairs(path):
    with gzip.open(path, "rt", encoding="ascii") as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise SystemExit(f"[ERROR] Incomplete FASTQ: {path}")
    return lines // 4
post = pairs(root / "post_fastp_1.fastq.gz")
legacy = pairs(root / "un_conc_1.fastq.gz")
strict = pairs(root / "both_unmapped_1.fastq.gz")
with output.open("w", encoding="utf-8") as handle:
    handle.write("sample_id\tpost_fastp_pairs\tun_conc_pairs\tboth_unmapped_pairs\textra_pairs_in_un_conc\textra_pct_of_post_fastp\n")
    handle.write(f"{sample}\t{post}\t{legacy}\t{strict}\t{legacy-strict}\t{100*(legacy-strict)/post:.8f}\n")
print(output.read_text(), end="")
PY

python3 - "$results" "$RECEIPT" <<'PY'
import hashlib, pathlib, sys
root, receipt = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2])
files = sorted(p for p in root.rglob("*") if p.is_file() and p.stat().st_size)
receipt.parent.mkdir(parents=True, exist_ok=True)
with receipt.open("w", encoding="utf-8") as out:
    out.write("path\tsha256\tbytes\n")
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        out.write(f"{path}\t{digest}\t{path.stat().st_size}\n")
PY
echo "[PASS] Host-policy retention audit completed: $SAMPLE_ID"
