#!/usr/bin/env bash
set -euo pipefail

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/opt/conda}"
if [[ -n "${MAMBA_BIN:-}" ]]; then
  [[ -x "$MAMBA_BIN" ]] || { echo "[ERROR] micromamba missing: $MAMBA_BIN" >&2; exit 1; }
else
  MAMBA_BIN="$(command -v micromamba || true)"
  if [[ -z "$MAMBA_BIN" ]]; then
    for candidate in /bin/micromamba /usr/bin/micromamba /usr/local/bin/micromamba; do
      if [[ -x "$candidate" ]]; then
        MAMBA_BIN="$candidate"
        break
      fi
    done
  fi
  [[ -n "$MAMBA_BIN" ]] || { echo "[ERROR] micromamba missing from the image" >&2; exit 1; }
fi

"$MAMBA_BIN" list -n taxonomic_tools --json |
  /opt/conda/envs/taxonomic_tools/bin/python -c '
import json, sys
expected = {
    "python": "3.11",
    "kraken2": "2.1.6",
    "bracken": "3.1",
    "metaphlan": "4.2.2",
    "bowtie2": "2.5.4",
    "art": "2016.06.05",
    "seqtk": "1.4",
    "fastqc": "0.12.1",
    "fastp": "0.23.4",
}
observed = {item["name"]: item["version"] for item in json.load(sys.stdin)}
errors = []
for name, version in expected.items():
    actual = observed.get(name)
    if actual is None or not actual.startswith(version):
        errors.append(f"{name}: expected {version}, found {actual}")
if errors:
    raise SystemExit("[ERROR] Pinned package verification failed:\n- " + "\n- ".join(errors))
print("[OK] Pinned package versions verified")
'

"$MAMBA_BIN" run -n taxonomic_tools bash -c '
  set -euo pipefail
  for tool in python3 fastqc fastp bowtie2 kraken2 bracken metaphlan art_illumina seqtk; do
    command -v "$tool" >/dev/null
    printf "[OK] %-12s %s\n" "$tool" "$(command -v "$tool")"
  done
'
