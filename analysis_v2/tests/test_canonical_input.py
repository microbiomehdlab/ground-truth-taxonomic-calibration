#!/usr/bin/env python3
"""Fixture tests for the canonical v2 analytical input contract."""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path


HEADER = [
    "schema_version", "cohort", "study", "sample_id", "condition",
    "analysis_population", "target_label", "target_taxon", "assembly_arm",
    "profiler", "profile_id", "baseline_profile_id", "spike_fraction_total",
    "spike_fraction_target", "implanted_read_pairs_target", "native_abundance",
    "native_unit", "abundance_fraction", "detected_native_nonzero",
    "source_profile", "source_design", "include", "exclusion_reason",
]


def row(profiler: str, profile: str, dose: str, native: str, pairs: str) -> list[str]:
    baseline = f"base-{profiler}"
    unit = "fraction_total_reads" if profiler == "kraken2_bracken" else "relative_abundance_pct"
    fraction = native if profiler == "kraken2_bracken" else str(float(native) / 100)
    is_baseline = dose == "0"
    return [
        "paired-dose-response-v2.1", "yachida", "YachidaS_2019", "S1", "CRC",
        "independent", "Fnuc", "Fusobacterium nucleatum", "original", profiler,
        baseline if is_baseline else profile, baseline, dose, dose, pairs, native,
        unit, fraction, "1" if float(fraction) > 0 else "0",
        f"profiles/{profile}.tsv", "BASELINE" if is_baseline else "design.tsv", "1", "",
    ]


def run_validator(script: Path, table: Path, outdir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["python3", str(script), "--input", str(table), "--outdir", str(outdir)],
        text=True, capture_output=True,
    )


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    script = repo / "analysis_v2/scripts/validate_canonical_input.py"
    with tempfile.TemporaryDirectory(prefix="canonical_input.") as name:
        root = Path(name)
        valid = root / "valid.tsv"
        rows = []
        for profiler, base, spike in [
            ("kraken2_bracken", "0.01", "0.011"),
            ("metaphlan4", "1.0", "1.2"),
        ]:
            rows.extend([
                row(profiler, f"spike-{profiler}", "0", base, "0"),
                row(profiler, f"spike-{profiler}", "0.001", spike, "1000"),
            ])
        with valid.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(rows)
        result = run_validator(script, valid, root / "valid_audit")
        assert result.returncode == 0, result.stderr
        assert (root / "valid_audit/SUCCESS").is_file()

        invalid = root / "invalid.tsv"
        bad_rows = [rows[1]]
        with invalid.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(bad_rows)
        result = run_validator(script, invalid, root / "invalid_audit")
        assert result.returncode != 0
        assert "lacks unique baseline" in result.stderr

    print("[PASS] canonical-input fixture tests")


if __name__ == "__main__":
    main()
