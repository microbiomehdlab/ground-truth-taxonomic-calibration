#!/usr/bin/env python3
"""Fixture tests for deterministic paired endpoint derivation."""

import csv
import math
import subprocess
import tempfile
from pathlib import Path

from test_canonical_input import HEADER, row


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    script = repo / "analysis_v2/scripts/derive_paired_endpoints.py"
    with tempfile.TemporaryDirectory(prefix="paired_endpoints.") as name:
        root = Path(name)
        table = root / "canonical.tsv"
        # Community example: total spike is 10%, target contributes 1%.
        baseline = row("kraken2_bracken", "community-profile", "0", "0.2", "0")
        perturbed = row("kraken2_bracken", "community-profile", "0.01", "0.19", "1000")
        total_index = HEADER.index("spike_fraction_total")
        perturbed[total_index] = "0.1"
        with table.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows([baseline, perturbed])

        result = subprocess.run([
            "python3", str(script), "--input", str(table),
            "--outdir", str(root / "derived"),
        ], text=True, capture_output=True)
        assert result.returncode == 0, result.stderr
        with (root / "derived/paired_endpoints.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        result_row = rows[0]
        assert math.isclose(float(result_row["baseline_retained_after_dilution"]), 0.18)
        assert math.isclose(float(result_row["read_proportional_reference"]), 0.19)
        assert math.isclose(float(result_row["recovered_spike_signal"]), 0.01)
        assert math.isclose(float(result_row["response_ratio"]), 1.0)
        assert math.isclose(float(result_row["signed_reference_error"]), 0.0, abs_tol=1e-12)
        assert (root / "derived/input_validation/SUCCESS").is_file()
        assert (root / "derived/SUCCESS").is_file()

    print("[PASS] paired-endpoint fixture tests")


if __name__ == "__main__":
    main()
