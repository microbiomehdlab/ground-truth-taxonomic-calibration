#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"

git rev-parse --is-inside-work-tree >/dev/null

bad=0
tracked="$(git ls-files -co --exclude-standard)"

check_tracked_pattern() {
  local pattern="$1"
  local description="$2"
  local hits
  hits="$(grep -E "$pattern" <<<"$tracked" || true)"
  if [[ -n "$hits" ]]; then
    echo "[ERROR] Tracked $description:" >&2
    echo "$hits" >&2
    bad=1
  fi
}

check_tracked_pattern '(^|/)(global|spikein|dataset)[.]env$' "private configuration"
check_tracked_pattern '[.](fastq|fq)([.]gz)?$|[.]sra$' "raw sequencing data"
check_tracked_pattern '[.](sif|img)$' "container images"
check_tracked_pattern '(^|/)(results|profile_results|runs|logs)/' "generated output"
check_tracked_pattern '[.](out|err)$|(^|/)slurm-' "scheduler logs"
check_tracked_pattern '[.](tar|tar[.]gz|zip)$' "archives"

host_paths="$(
  git ls-files -co --exclude-standard -z |
    xargs -0 -r grep -IlE '^[[:space:]]*[^#[:space:]].*(/home/|/mnt/|/scratch/)' 2>/dev/null || true
)"
host_paths="$(grep -vE '^scripts/check_repository[.]' <<<"$host_paths" || true)"
if [[ -n "$host_paths" ]]; then
  echo "[ERROR] Tracked files contain host-specific absolute paths:" >&2
  echo "$host_paths" >&2
  bad=1
fi

large="$(
  git ls-files -co --exclude-standard -z |
    xargs -0 -r stat -c '%s	%n' |
    awk '$1 > 5000000'
)"
if [[ -n "$large" ]]; then
  echo "[ERROR] Tracked files exceed 5 MB:" >&2
  echo "$large" >&2
  bad=1
fi

stale_publication_refs="$(
  git grep --untracked -nE \
    'maaslin2_1180_v1[.]sif|TableA6_target_database_representation|Supplementary Table A6' \
    -- '*.md' '*.sh' '*.py' 2>/dev/null |
    grep -v '^scripts/check_repository[.]sh:' || true
)"
if [[ -n "$stale_publication_refs" ]]; then
  echo "[ERROR] Stale publication workflow references:" >&2
  echo "$stale_publication_refs" >&2
  bad=1
fi

if ((bad)); then
  echo "[FAIL] Repository publication checks failed." >&2
  exit 1
fi

echo "[PASS] Repository publication checks passed."
