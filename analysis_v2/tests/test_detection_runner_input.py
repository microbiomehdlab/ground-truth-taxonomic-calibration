#!/usr/bin/env python3
"""Prevent definitive runners from dropping zero-dose detection rows."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
RUNNERS = (
    ROOT / "analysis_v2/run_crc_cohort_definitive_analysis.sh",
    ROOT / "analysis_v2/run_yachida_definitive_analysis.sh",
)


for runner in RUNNERS:
    text = runner.read_text(encoding="utf-8")
    calls = re.findall(
        r"fit_detection_dose_response\.R\s+\\\n\s+--input\s+([^ ]+)", text
    )
    assert calls, f"No detection-model invocation found in {runner.name}"
    assert calls == ['"$CANONICAL_INPUT"'], (
        f"{runner.name} must pass canonical input (including dose zero) "
        f"to the detection model; found {calls}"
    )

print("[PASS] definitive detection runners retain canonical zero-dose rows")
