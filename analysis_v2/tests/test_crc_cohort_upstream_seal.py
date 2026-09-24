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
    root = Path(tmp); state = root / "state/samples"; results = root / "results"; qc = root / "qc"; scratch = root / "scratch"
    state.mkdir(parents=True); qc.mkdir(); scratch.mkdir()
    manifest = root / "production.tsv"; independent = root / "independent.tsv"
    header = "sample_id\tcondition\tstudy\tage\tsex\tbmi\tindependent_subset\n"
    manifest.write_text(header + "S1\tControl\tStudy\t50\tFemale\t20\t1\nS2\tCRC\tStudy\t60\tMale\t25\t0\n")
    independent.write_text(header + "S1\tControl\tStudy\t50\tFemale\t20\t1\n")
    sidecar(manifest); sidecar(independent)
    payloads = {}
    for sample, count in (("S1", 68), ("S2", 8)):
        for suffix in ("verified", "input_provenance.tsv"):
            (state / f"{sample}.{suffix}").write_text("status\tPASS\n")
        sample_root = results / "Study" / sample
        (sample_root / "profiles").mkdir(parents=True)
        # Production sentinel markers are deliberately zero-byte files.
        (sample_root / "SUCCESS").touch()
        for index in range(count):
            marker = sample_root / "profiles" / str(index) / "SUCCESS"; marker.parent.mkdir(); marker.touch()
        payload = sample_root / "sample_completion.tsv"
        payload.write_text("field\tvalue\nstatus\tPASS\n")
        payloads[sample] = payload
        (state / f"{sample}.retained_outputs.tsv").write_text(
            "path\tsha256\tbytes\n"
            f"{payload}\t{hashlib.sha256(payload.read_bytes()).hexdigest()}\t{payload.stat().st_size}\n"
        )
    out = root / "seal"
    subprocess.run([sys.executable, str(repo / "analysis_v2/scripts/seal_crc_cohort_upstream.py"),
                    "--cohort", "feng", "--manifest", str(manifest),
                    "--independent-manifest", str(independent), "--state-dir", str(root / "state"),
                    "--results-root", str(results), "--qc-root", str(qc),
                    "--scratch-root", str(scratch), "--outdir", str(out),
                    "--expected-samples", "2", "--expected-independent", "1"], check=True)
    assert (out / "SUCCESS").is_file() and (out / "production_seal.sha256").is_file()
    assert "S1" in (out / "sample_flow.tsv").read_text()
    assert not (out / "AUDIT_IN_PROGRESS").exists()

    # A changed retained output must invalidate an existing seal and be visible
    # in the sample-flow ledger rather than being accepted on marker presence.
    payloads["S2"].write_text("field\tvalue\nstatus\tFAIL\n")
    failed = subprocess.run(
        [sys.executable, str(repo / "analysis_v2/scripts/seal_crc_cohort_upstream.py"),
         "--cohort", "feng", "--manifest", str(manifest),
         "--independent-manifest", str(independent), "--state-dir", str(root / "state"),
         "--results-root", str(results), "--qc-root", str(qc),
         "--scratch-root", str(scratch), "--outdir", str(out),
         "--expected-samples", "2", "--expected-independent", "1"],
        text=True, capture_output=True,
    )
    assert failed.returncode != 0
    assert not (out / "SUCCESS").exists()
    assert not (out / "production_seal.sha256").exists()
    assert (out / "AUDIT_IN_PROGRESS").is_file()
    assert "sha256_mismatch" in (out / "sample_flow.tsv").read_text()
print("[PASS] CRC cohort upstream-seal fixture")
