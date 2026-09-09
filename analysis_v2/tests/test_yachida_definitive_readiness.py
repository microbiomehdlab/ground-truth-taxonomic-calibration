#!/usr/bin/env python3
import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


repo = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp); state = root / "state"; seal = state / "production_seal"
    (state / "samples").mkdir(parents=True); seal.mkdir()
    manifest = root / "manifest.tsv"
    manifest.write_text("sample_id\tTarget_Condition\nS1\tControl\nS2\tCRC\n", encoding="utf-8")
    for sample in ("S1", "S2"):
        (state / "samples" / f"{sample}.verified").write_text("status\tPASS\n")
        (state / "samples" / f"{sample}.retained_outputs.tsv").write_text("path\tsha256\tbytes\n")
    for name in ("dataset_completion.tsv", "dataset_completion.tsv.sha256", "pilot_batched.tsv", "SUCCESS"):
        (seal / name).write_text(f"{name}\n")
    sealed = [seal / name for name in ("dataset_completion.tsv", "dataset_completion.tsv.sha256", "pilot_batched.tsv", "SUCCESS")]
    (seal / "production_seal.sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in sealed))
    canonical = root / "canonical.tsv"
    fields = ["cohort", "sample_id", "analysis_population", "profiler", "include"]
    data = []
    for sample in ("S1", "S2"):
        for profiler in ("kraken2_bracken", "metaphlan4"):
            data.append(["yachida", sample, "community", profiler, "1"])
    for profiler in ("kraken2_bracken", "metaphlan4"):
        data.append(["yachida", "S1", "independent", profiler, "1"])
    with canonical.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(fields); writer.writerows(data)
    success = root / "canonical.SUCCESS"; success.write_text("status\tPASS\n")
    image = root / "analysis.sif"; image.write_text("fixture\n")
    sensitivity = root / "sensitivity.SUCCESS"; sensitivity.write_text("status\tPASS\n")
    out = root / "out"
    command = [sys.executable, str(repo / "analysis_v2/scripts/check_yachida_definitive_readiness.py"),
        "--manifest", str(manifest), "--state-dir", str(state), "--canonical", str(canonical),
        "--canonical-success", str(success), "--analysis-sif", str(image),
        "--assembly-sensitivity-success", str(sensitivity), "--expected-samples", "2",
        "--expected-independent", "1", "--outdir", str(out)]
    subprocess.run(command, check=True)
    assert (out / "SUCCESS").is_file() and (out / "readiness.sha256").is_file()
    canonical.write_text(canonical.read_text().replace("independent", "community"))
    failed = subprocess.run(command[:-1] + [str(root / "bad")], capture_output=True, text=True)
    assert failed.returncode != 0 and "expected community and independent" in failed.stderr
print("[PASS] Yachida definitive-readiness fixture")
