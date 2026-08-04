#!/usr/bin/env bash
set -euo pipefail

ENV_PREFIX="/opt/conda/envs/taxonomic_tools"
PYTHON="$ENV_PREFIX/bin/python"
[[ -x "$PYTHON" ]] || { echo "[ERROR] Environment Python missing: $PYTHON" >&2; exit 1; }

"$PYTHON" -c '
import glob, json
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
metadata = glob.glob("/opt/conda/envs/taxonomic_tools/conda-meta/*.json")
observed = {}
for filename in metadata:
    with open(filename, encoding="utf-8") as handle:
        item = json.load(handle)
    observed[item["name"]] = item["version"]
errors = []
for name, version in expected.items():
    actual = observed.get(name)
    if actual is None or not actual.startswith(version):
        errors.append(f"{name}: expected {version}, found {actual}")
if errors:
    raise SystemExit("[ERROR] Pinned package verification failed:\n- " + "\n- ".join(errors))
print("[OK] Pinned package versions verified")
'

export PATH="$ENV_PREFIX/bin:$PATH"
for tool in python3 fastqc fastp bowtie2 kraken2 bracken metaphlan art_illumina seqtk; do
  command -v "$tool" >/dev/null
  printf "[OK] %-12s %s\n" "$tool" "$(command -v "$tool")"
done
