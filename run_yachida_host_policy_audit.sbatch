#!/usr/bin/env bash
#SBATCH --job-name=yachida_host_policy
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=yachida_host_policy.%j.out
#SBATCH --error=yachida_host_policy.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export absolute config/yachida.env path}"
: "${SAMPLE_ID:?export one sample ID}"
manifest="${YACHIDA_MANIFEST:-$ROOT/work/yachida_67x3/metadata/pilot_batched.tsv}"
python3 "$ROOT/datasets/yachida/stream_sample.py" \
  --manifest "$manifest" --sample-id "$SAMPLE_ID" \
  --scratch-root "${HOST_POLICY_AUDIT_SCRATCH:-$ROOT/work/private_host_policy_audit/scratch}" \
  --state-dir "${HOST_POLICY_AUDIT_STATE:-$ROOT/work/private_host_policy_audit/state}" \
  --runner "$ROOT/datasets/yachida/run_host_policy_audit.sh"
