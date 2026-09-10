#!/usr/bin/env python3
import csv
import hashlib
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_yachida_production_audit_seals_verified_receipts():
    with tempfile.TemporaryDirectory() as name:
        root = pathlib.Path(name)
        state = root / "state"
        scratch = root / "scratch"
        retained = root / "results" / "result.tsv"
        retained.parent.mkdir(parents=True)
        retained.write_text("result\n", encoding="utf-8")
        samples = state / "samples"
        samples.mkdir(parents=True)
        (samples / "S1.verified").write_text("outputs\t1\n", encoding="utf-8")
        (samples / "S1.retained_outputs.tsv").write_text(
            f"path\tsha256\tbytes\n{retained}\t{sha256(retained)}\t{retained.stat().st_size}\n",
            encoding="utf-8",
        )
        manifest = root / "manifest.tsv"
        with manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["sample_id", "Target_Condition", "batch_id"],
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow({"sample_id": "S1", "Target_Condition": "Control", "batch_id": "batch_001"})
        independent = root / "independent.tsv"
        independent.write_text("sample_id\nS1\n", encoding="utf-8")
        sample_root = root / "results/YachidaS_2019/S1"
        for group, count in (("baseline", 1), ("independent", 60), ("community", 7)):
            for index in range(count):
                marker = sample_root / "profiles" / group / f"profile_{index}" / "SUCCESS"
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.touch()
        (sample_root / "sample_completion.tsv").write_text(
            "field\tvalue\n"
            "sample_id\tS1\n"
            "condition\tControl\n"
            "independent_subset\t1\n"
            "expected_profiles\t68\n"
            "observed_profiles\t68\n"
            "community_design_rows\t7\n"
            "independent_design_rows\t60\n",
            encoding="utf-8",
        )
        subprocess.run([
            "python3", str(ROOT / "datasets/yachida/audit_production.py"),
            "--manifest", str(manifest), "--scratch-root", str(scratch),
            "--results-root", str(root / "results"),
            "--independent-manifest", str(independent),
            "--state-dir", str(state), "--expected-samples", "1",
        ], check=True)
        assert (state / "batches/batch_001/SUCCESS").is_file()
        assert (state / "production_seal/SUCCESS").is_file()
        assert (state / "production_seal/production_seal.sha256").is_file()


def test_production_gate_shell_syntax():
    for path in [
        ROOT / "datasets/preflight_crc_production.sh",
        ROOT / "run_yachida_production_audit.sbatch",
    ]:
        subprocess.run(["bash", "-n", str(path)], check=True)
