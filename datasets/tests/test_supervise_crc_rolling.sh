#!/usr/bin/env bash
set -euo pipefail

DATASETS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin" "$tmp/state"

cat > "$tmp/jobs.tsv" <<'EOF'
generation	offset	tasks	download_job	compute_job
1	0	4	100	101
2	4	4	200	201
3	8	2	300	301
EOF

cat > "$tmp/bin/squeue" <<'EOF'
#!/usr/bin/env bash
case "$*" in
  *"-t R,CF,CG"*) printf '100_1\n100_2\n' ;;
  *"-j 200 -t PD"*)
    printf '200|1|Dependency\n200|2|DependencyNeverSatisfied\n200|3|Dependency\n'
    ;;
  *"-j 300 -t PD"*) printf '300|1|Dependency\n' ;;
esac
EOF

cat > "$tmp/bin/sacct" <<'EOF'
#!/usr/bin/env bash
case "$*" in
  *"-j 101_1 "*) printf 'COMPLETED|\n' ;;
  *"-j 101_2 "*) printf 'FAILED|\n' ;;
  *"-j 101_3 "*) printf 'RUNNING|\n' ;;
  *"-j 201_1 "*) printf 'COMPLETED|\n' ;;
  *) printf 'PENDING|\n' ;;
esac
EOF

cat > "$tmp/bin/scontrol" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SCONTROL_CALLS"
EOF
chmod +x "$tmp/bin/"*

export PATH="$tmp/bin:$PATH"
export SCONTROL_CALLS="$tmp/scontrol.calls"

bash "$DATASETS_DIR/supervise_crc_rolling.sh" \
  --ledger "$tmp/jobs.tsv" --download-concurrent 3 \
  --state-dir "$tmp/state" --once

grep -Fx 'update JobId=200_1 Dependency=' "$SCONTROL_CALLS"
if grep -Eq 'JobId=(200_2|200_3|300_1)' "$SCONTROL_CALLS"; then
  echo '[FAIL] Supervisor released a task without a completed direct predecessor' >&2
  exit 1
fi
grep -F $'200_1\t101_1\tCOMPLETED' "$tmp/state/releases.tsv"

# A restarted supervisor must not release a task already recorded in state.
: > "$SCONTROL_CALLS"
bash "$DATASETS_DIR/supervise_crc_rolling.sh" \
  --ledger "$tmp/jobs.tsv" --download-concurrent 3 \
  --state-dir "$tmp/state" --once
test ! -s "$SCONTROL_CALLS"

echo '[PASS] CRC rolling supervisor admission and restart test'
