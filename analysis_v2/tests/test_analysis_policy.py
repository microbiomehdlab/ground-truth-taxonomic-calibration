#!/usr/bin/env python3
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "policy"
    subprocess.run(["python3", str(ROOT / "analysis_v2/scripts/validate_analysis_policy.py"),
                    "--policy", str(ROOT / "analysis_v2/ANALYSIS_POLICY.tsv"),
                    "--outdir", str(out)], check=True)
    assert (out / "SUCCESS").is_file()
    subprocess.run(["sha256sum", "-c", "analysis_policy.sha256"], cwd=out, check=True)
print("[PASS] frozen analysis-policy fixture")
