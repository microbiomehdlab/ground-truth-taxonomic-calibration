#!/usr/bin/env bash
set -euo pipefail
DATASETS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin" "$tmp/state"
cat > "$tmp/jobs.tsv" <<'EOF'
generation	offset	tasks	download_job	compute_job
1	0	4	100	101
2	4	4	200	201
EOF
cat > "$tmp/bin/sacct" <<'EOF'
#!/usr/bin/env bash
case "$*" in
  *"-j 100,200 "*) cat <<'OUT'
100_1|COMPLETED|
100_2|COMPLETED|
100_3|COMPLETED|
100_4|COMPLETED|
200_1|PENDING|
200_2|PENDING|
200_3|PENDING|
200_4|PENDING|
OUT
  ;;
  *"-j 101,201 "*) cat <<'OUT'
101_1|RUNNING|
101_2|RUNNING|
101_3|RUNNING|
101_4|COMPLETED|
201_1|PENDING|
201_2|PENDING|
201_3|PENDING|
201_4|PENDING|
OUT
  ;;
  *"-j 900,901 "*) printf '900_5|RUNNING|\n901_5|PENDING|\n' ;;
  *) : ;;
esac
EOF
cat > "$tmp/bin/squeue" <<'EOF'
#!/usr/bin/env bash
case "$*" in
  *"-j 100,200 "*) printf '200|1|PENDING|Dependency\n200|2|PENDING|Dependency\n200|3|PENDING|Dependency\n200|4|PENDING|Dependency\n' ;;
  *"-j 101,201 "*) printf '101|1|RUNNING|None\n101|2|RUNNING|None\n101|3|RUNNING|None\n201|1|PENDING|Dependency\n201|2|PENDING|Dependency\n201|3|PENDING|Dependency\n201|4|PENDING|Dependency\n' ;;
  *"-j 900,901 "*) printf '900|5|RUNNING|None\n901|5|PENDING|Dependency\n' ;;
  *"-j 200 -t PD"*) printf '200|1|Dependency\n200|2|Dependency\n200|3|Dependency\n200|4|Dependency\n' ;;
esac
EOF
cat > "$tmp/bin/scontrol" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SCONTROL_CALLS"
EOF
chmod +x "$tmp/bin/"*
export PATH="$tmp/bin:$PATH" SCONTROL_CALLS="$tmp/calls"

# Three running computes occupy 3/4 slots: admit exactly one download.
bash "$DATASETS_DIR/supervise_crc_global_slots.sh" --ledger "$tmp/jobs.tsv" \
  --inflight-limit 4 --download-concurrent 2 --state-dir "$tmp/state" --once
test "$(wc -l < "$SCONTROL_CALLS")" -eq 1
grep -Fx 'update JobId=200_1 Dependency=' "$SCONTROL_CALLS"
grep -F $'2:1\t200_1\tDependency' "$tmp/state/admissions.tsv"

# Persistent admission occupies the fourth slot after restart: admit nothing.
: > "$SCONTROL_CALLS"
bash "$DATASETS_DIR/supervise_crc_global_slots.sh" --ledger "$tmp/jobs.tsv" \
  --inflight-limit 4 --download-concurrent 2 --state-dir "$tmp/state" --once
test ! -s "$SCONTROL_CALLS"

# With a 5-slot pool, the same observed state admits one additional sample.
bash "$DATASETS_DIR/supervise_crc_global_slots.sh" --ledger "$tmp/jobs.tsv" \
  --inflight-limit 5 --download-concurrent 2 --state-dir "$tmp/state" --once
grep -Fx 'update JobId=200_2 Dependency=' "$SCONTROL_CALLS"

# A completed download whose compute remains dependency-blocked does not occupy
# a runnable sample slot. Future/dependency-held work must not masquerade as an
# active lane.
rm -rf "$tmp/dependency-held-state"
: > "$SCONTROL_CALLS"
bash "$DATASETS_DIR/supervise_crc_global_slots.sh" --ledger "$tmp/jobs.tsv" \
  --inflight-limit 4 --download-concurrent 4 \
  --state-dir "$tmp/dependency-held-state" --once
test "$(wc -l < "$SCONTROL_CALLS")" -eq 1
grep -Fx 'update JobId=200_1 Dependency=' "$SCONTROL_CALLS"

# A retry pipeline occupies one sample slot and one download slot. With three
# original computes, the 4-slot pool is full and admits nothing.
: > "$SCONTROL_CALLS"
bash "$DATASETS_DIR/supervise_crc_global_slots.sh" --ledger "$tmp/jobs.tsv" \
  --inflight-limit 4 --download-concurrent 2 --external-pipeline 900:901 \
  --state-dir "$tmp/external-state" --once
test ! -s "$SCONTROL_CALLS"

echo '[PASS] global CRC sample-slot admission tests'
