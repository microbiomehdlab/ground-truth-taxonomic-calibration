# Automatic CRC rolling-lane admission

`supervise_crc_rolling.sh` fills otherwise idle download capacity across rolling
generations without manually editing Slurm dependencies. It implements the
same guarded recovery described in `CRC_ROLLING_PRODUCTION.md`:

- the global number of running or supervisor-admitted downloads is bounded;
- only download tasks are released;
- a task is eligible only after the preceding compute task in the same lane is
  exactly `COMPLETED`;
- compute dependencies are never changed;
- releases are persisted in `supervisor/releases.tsv` and survive restarts;
- a file lock prevents two supervisors from controlling the same ledger.

The supervisor deliberately does not rescue failed predecessors. Operational
failures remain visible and require a documented targeted retry.

## Start a new rolling cohort

Run the cohort preflight and submit the rolling DAG normally. Capture the
`jobs.tsv` path printed by `submit_crc_rolling.sh`, then perform a read-only
admission preview:

```bash
bash datasets/supervise_crc_rolling.sh \
  --ledger /absolute/path/to/jobs.tsv \
  --download-concurrent 8 \
  --once --dry-run
```

Inspect the proposed tasks. Every displayed predecessor must be the compute
array from the immediately preceding ledger generation and must be
`COMPLETED`.

For unattended operation on a login host where lightweight polling processes
are permitted:

```bash
ledger=/absolute/path/to/jobs.tsv
supervisor_dir="$(dirname "$ledger")/supervisor"
mkdir -p "$supervisor_dir"
nohup bash datasets/supervise_crc_rolling.sh \
  --ledger "$ledger" \
  --download-concurrent 8 \
  --poll-seconds 60 \
  > "$supervisor_dir/nohup.log" 2>&1 < /dev/null &
echo "$!" > "$supervisor_dir/pid"
```

If persistent processes are prohibited on login hosts, run the same foreground
command inside the cluster-approved persistent-session or utility-job mechanism.
The controller performs only brief `squeue`, `sacct`, and `scontrol` calls and
sleeps between checks; it does not process biological data.

## Monitor and stop

```bash
tail -f "$supervisor_dir/events.log"
column -t -s $'\t' "$supervisor_dir/releases.tsv"
```

After all download generations finish, stop the polling process and retain the
supervisor directory with the production evidence:

```bash
kill "$(cat "$supervisor_dir/pid")"
```

Restarting with the same ledger and state directory is idempotent: previously
released tasks remain accounted for and will not be released again.
