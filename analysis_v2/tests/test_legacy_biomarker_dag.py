#!/usr/bin/env python3
import os
import pathlib
import subprocess
import tempfile

root = pathlib.Path(__file__).resolve().parents[2]
submitter = root / "analysis_v2/submit_legacy_biomarker_development_dag.sh"
worker = root / "analysis_v2/run_legacy_biomarker_development_stage.sbatch"

subprocess.run(["bash", "-n", str(submitter)], check=True)
subprocess.run(["bash", "-n", str(worker)], check=True)

with tempfile.TemporaryDirectory() as temporary:
    temporary = pathlib.Path(temporary)
    fake_sbatch = temporary / "sbatch"
    counter = temporary / "counter"
    calls = temporary / "calls"
    fake_sbatch.write_text(
        "#!/usr/bin/env bash\n"
        f"n=$(cat '{counter}' 2>/dev/null || echo 1000)\n"
        "n=$((n + 1))\n"
        f"printf '%s' \"$n\" > '{counter}'\n"
        f"printf '%s\\n' \"$*\" >> '{calls}'\n"
        "printf '%s\\n' \"$n\"\n",
        encoding="utf-8",
    )
    fake_sbatch.chmod(0o755)
    output = temporary / "output"
    env = os.environ | {
        "DEV_INPUT_ROOT": str(temporary / "input"),
        "ANALYSIS_SIF": str(temporary / "analysis.sif"),
        "OUTDIR": str(output),
        "SBATCH_BIN": str(fake_sbatch),
    }
    result = subprocess.run(
        ["bash", str(submitter)], cwd=root, env=env, text=True,
        capture_output=True, check=True,
    )
    submitted = calls.read_text(encoding="utf-8").splitlines()
    assert len(submitted) == 8
    assert "--dependency=afterok:1001" in submitted[1]
    assert "--dependency=afterok:1001" in submitted[2]
    assert "--dependency=afterok:1001" in submitted[3]
    assert "--dependency=afterok:1002:1003" in submitted[4]
    assert "--dependency=afterok:1004" in submitted[5]
    assert "--dependency=afterok:1002:1003" in submitted[6]
    assert "--dependency=afterok:1005:1006:1007" in submitted[7]
    assert "seal=1008" in result.stdout

print("[PASS] legacy biomarker-development DAG fixture")
