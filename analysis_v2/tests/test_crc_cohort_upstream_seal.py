#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

repo = Path(__file__).resolve().parents[2]


def sidecar(path: Path) -> None:
    path.with_suffix(path.suffix + ".sha256").write_text(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n", encoding="utf-8"
    )


def write_receipt(path: Path, payload: Path) -> None:
    path.write_text(
        "path\tsha256\tbytes\n"
        f"{payload}\t{hashlib.sha256(payload.read_bytes()).hexdigest()}\t{payload.stat().st_size}\n"
    )


def audit_command(root: Path, manifest: Path, independent: Path, out: Path) -> list[str]:
    return [
        sys.executable, str(repo / "analysis_v2/scripts/seal_crc_cohort_upstream.py"),
        "--cohort", "feng", "--manifest", str(manifest),
        "--independent-manifest", str(independent), "--state-dir", str(root / "state"),
        "--results-root", str(root / "results"), "--qc-root", str(root / "qc"),
        "--scratch-root", str(root / "scratch"), "--outdir", str(out),
        "--expected-samples", "2", "--expected-independent", "1",
    ]


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    state = root / "state/samples"
    results = root / "results"
    qc = root / "qc"
    scratch = root / "scratch"
    state.mkdir(parents=True)
    qc.mkdir()
    scratch.mkdir()
    manifest = root / "production.tsv"
    independent = root / "independent.tsv"
    header = "sample_id\tcondition\tstudy\tage\tsex\tbmi\tindependent_subset\n"
    manifest.write_text(
        header + "S1\tControl\tStudy\t50\tFemale\t20\t1\n"
        + "S2\tCRC\tStudy\t60\tMale\t25\t0\n"
    )
    independent.write_text(header + "S1\tControl\tStudy\t50\tFemale\t20\t1\n")
    sidecar(manifest)
    sidecar(independent)
    payloads: dict[str, Path] = {}
    community_markers: dict[str, list[Path]] = {}
    for sample, condition, independent_count in (("S1", "Control", 60), ("S2", "CRC", 0)):
        for suffix in ("verified", "input_provenance.tsv"):
            (state / f"{sample}.{suffix}").write_text("status\tPASS\n")
        sample_root = results / "Study" / sample
        sample_root.mkdir(parents=True)
        (sample_root / "SUCCESS").touch()
        baseline = sample_root / "profiles/baseline" / sample / "SUCCESS"
        baseline.parent.mkdir(parents=True)
        baseline.touch()
        for index in range(independent_count):
            marker = sample_root / "profiles/independent/target" / str(index) / "SUCCESS"
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.touch()
        community_markers[sample] = []
        for index in range(7):
            marker = sample_root / "profiles/community" / str(index) / "SUCCESS"
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.touch()
            community_markers[sample].append(marker)
        expected_profiles = 1 + independent_count + 7
        payload = sample_root / "sample_completion.tsv"
        payload.write_text(
            "field\tvalue\n" f"sample_id\t{sample}\n" "study\tStudy\n"
            f"condition\t{condition}\n" f"independent_subset\t{1 if independent_count else 0}\n"
            f"expected_profiles\t{expected_profiles}\n" f"observed_profiles\t{expected_profiles}\n"
            "community_design_rows\t7\n" f"independent_design_rows\t{independent_count}\n"
        )
        payloads[sample] = payload
        write_receipt(state / f"{sample}.retained_outputs.tsv", payload)

    out = root / "seal"
    subprocess.run(audit_command(root, manifest, independent, out), check=True)
    assert (out / "SUCCESS").is_file() and (out / "production_seal.sha256").is_file()
    flow = (out / "sample_flow.tsv").read_text()
    assert "observed_baseline_profiles" in flow
    assert "observed_independent_profiles" in flow
    assert "observed_community_profiles" in flow
    assert not (out / "AUDIT_IN_PROGRESS").exists()

    payloads["S2"].write_text(payloads["S2"].read_text().replace("CRC", "BAD"))
    failed = subprocess.run(audit_command(root, manifest, independent, out), text=True, capture_output=True)
    assert failed.returncode != 0
    assert not (out / "SUCCESS").exists()
    assert not (out / "production_seal.sha256").exists()
    assert (out / "AUDIT_IN_PROGRESS").is_file()
    assert "sha256_mismatch" in (out / "sample_flow.tsv").read_text()

    payloads["S2"].write_text(payloads["S2"].read_text().replace("BAD", "CRC"))
    write_receipt(state / "S2.retained_outputs.tsv", payloads["S2"])
    community_markers["S2"][0].unlink()
    wrong = results / "Study/S2/profiles/baseline/extra/SUCCESS"
    wrong.parent.mkdir()
    wrong.touch()
    topology_failure = subprocess.run(
        audit_command(root, manifest, independent, out), text=True, capture_output=True
    )
    assert topology_failure.returncode != 0
    flow = (out / "sample_flow.tsv").read_text()
    assert "baseline_profile_count;community_profile_count" in flow

print("[PASS] CRC cohort upstream-seal fixture")
