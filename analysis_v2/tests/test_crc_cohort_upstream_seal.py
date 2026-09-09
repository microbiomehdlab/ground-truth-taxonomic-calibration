#!/usr/bin/env python3
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

repo = Path(__file__).resolve().parents[2]


def sidecar(path: Path) -> None:
    path.with_suffix(path.suffix + ".sha256").write_text(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n", encoding="utf-8")


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp); state = root / "state/samples"; results = root / "results"; state.mkdir(parents=True)
    manifest = root / "production.tsv"; independent = root / "independent.tsv"
    header = "sample_id\tcondition\tstudy\tage\tsex\tbmi\tindependent_subset\n"
    manifest.write_text(header + "S1\tControl\tStudy\t50\tFemale\t20\t1\nS2\tCRC\tStudy\t60\tMale\t25\t0\n")
    independent.write_text(header + "S1\tControl\tStudy\t50\tFemale\t20\t1\n")
    sidecar(manifest); sidecar(independent)
    for sample, count in (("S1", 68), ("S2", 8)):
        for suffix in ("verified", "retained_outputs.tsv", "input_provenance.tsv"):
            (state / f"{sample}.{suffix}").write_text("status\tPASS\n")
        sample_root = results / "Study" / sample
        (sample_root / "profiles").mkdir(parents=True)
        (sample_root / "SUCCESS").write_text("status\tPASS\n")
        for index in range(count):
            marker = sample_root / "profiles" / str(index) / "SUCCESS"; marker.parent.mkdir(); marker.write_text("PASS\n")
    out = root / "seal"
    subprocess.run([sys.executable, str(repo / "analysis_v2/scripts/seal_crc_cohort_upstream.py"),
                    "--cohort", "feng", "--manifest", str(manifest),
                    "--independent-manifest", str(independent), "--state-dir", str(root / "state"),
                    "--results-root", str(results), "--outdir", str(out),
                    "--expected-samples", "2", "--expected-independent", "1"], check=True)
    assert (out / "SUCCESS").is_file() and (out / "production_seal.sha256").is_file()
    assert "S1" in (out / "sample_flow.tsv").read_text()
print("[PASS] CRC cohort upstream-seal fixture")
