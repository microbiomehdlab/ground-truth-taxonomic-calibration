#!/usr/bin/env bash
# Select several independently seeded paired subsets while reading each pool mate once.
set -euo pipefail
IFS=$'\n\t'

[[ $# -eq 6 ]] || { echo "Usage: $0 IMAGE POOL1 POOL2 REQUESTS OUTDIR TMPDIR" >&2; exit 2; }
IMG="$1"; POOL1="$2"; POOL2="$3"; REQUESTS="$4"; OUTDIR="$5"; TMPDIR="$6"
for f in "$IMG" "$POOL1" "$POOL2" "$REQUESTS"; do [[ -s "$f" ]] || { echo "[ERROR] Missing/empty: $f" >&2; exit 1; }; done
mkdir -p "$OUTDIR" "$TMPDIR"
OUTDIR="$(cd "$OUTDIR" && pwd -P)"; TMPDIR="$(cd "$TMPDIR" && pwd -P)"

declare -a tags=() fifos=() pids=()
while IFS=$'\t' read -r tag seed n; do
  [[ "$tag" == "tag" || -z "$tag" ]] && continue
  [[ "$tag" =~ ^[A-Za-z0-9._-]+$ && "$seed" =~ ^[0-9]+$ && "$n" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid request: $tag $seed $n" >&2; exit 1; }
  tags+=("$tag"); fifo="$TMPDIR/${tag}.r1.fifo"; rm -f "$fifo"; mkfifo "$fifo"; fifos+=("$fifo")
  apptainer exec --cleanenv -B "$TMPDIR:$TMPDIR" "$IMG" seqtk sample -s "$seed" "$fifo" "$n" > "$OUTDIR/${tag}_1.fq" & pids+=("$!")
done < "$REQUESTS"
(( ${#tags[@]} > 0 )) || { echo "[ERROR] No sampling requests" >&2; exit 1; }
tee "${fifos[@]}" < "$POOL1" > /dev/null
for p in "${pids[@]}"; do wait "$p"; done

fifos=(); pids=()
for tag in "${tags[@]}"; do
  ids1="$TMPDIR/${tag}.ids1"; ids2="$TMPDIR/${tag}.ids2"
  awk 'NR%4==1 {gsub(/^@/,"",$1); print $1}' "$OUTDIR/${tag}_1.fq" > "$ids1"
  awk '{o=$1; print o; id=o; if(id~/\/1$/){sub(/\/1$/,"/2",id); print id} else if(id~/_1$/){sub(/_1$/,"_2",id); print id} else if(id~/[.]1$/){sub(/[.]1$/,".2",id); print id}}' "$ids1" | awk '!seen[$0]++' > "$ids2"
  fifo="$TMPDIR/${tag}.r2.fifo"; rm -f "$fifo"; mkfifo "$fifo"; fifos+=("$fifo")
  apptainer exec --cleanenv -B "$TMPDIR:$TMPDIR" "$IMG" seqtk subseq "$fifo" "$ids2" > "$OUTDIR/${tag}_2.fq" & pids+=("$!")
done
tee "${fifos[@]}" < "$POOL2" > /dev/null
for p in "${pids[@]}"; do wait "$p"; done

while IFS=$'\t' read -r tag _seed n; do
  [[ "$tag" == "tag" || -z "$tag" ]] && continue
  n1=$(awk 'END{print int(NR/4)}' "$OUTDIR/${tag}_1.fq"); n2=$(awk 'END{print int(NR/4)}' "$OUTDIR/${tag}_2.fq")
  [[ "$n1" -eq "$n" && "$n2" -eq "$n" ]] || { echo "[ERROR] $tag requested=$n R1=$n1 R2=$n2" >&2; exit 1; }
done < "$REQUESTS"
rm -f -- "$TMPDIR"/*.fifo
echo "[OK] Single-pass paired sampling completed for ${#tags[@]} request(s)"
