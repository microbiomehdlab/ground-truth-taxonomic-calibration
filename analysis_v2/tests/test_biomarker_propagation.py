#!/usr/bin/env python3
"""Fixture test for biomarker-propagation metrics."""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path


HEADER = ["cohort", "study", "analysis_population", "target_label", "assembly_arm",
          "profiler", "contrast", "spike_fraction_target", "feature", "effect",
          "p_value", "q_value", "include", "exclusion_reason"]


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    script = repo / "analysis_v2/scripts/evaluate_biomarker_propagation.py"
    aliases = repo / "examples/spike_taxon_aliases.csv"
    panel = repo / "spikes/spike_panel.tsv"
    with tempfile.TemporaryDirectory(prefix="biomarker_propagation.") as name:
        root = Path(name)
        calls = root / "calls.tsv"
        outdir = root / "output"
        prefix = ["yachida", "Study1", "independent", "Pana", "original",
                  "kraken2_bracken", "CRC_vs_Control"]
        rows = [
            prefix + ["0", "Peptostreptococcus anaerobius", "0.1", "0.2", "0.2", "1", ""],
            prefix + ["0", "Other species", "0.2", "0.01", "0.04", "1", ""],
            prefix + ["0.01", "Peptostreptococcus anaerobius", "0.7", "0.001", "0.01", "1", ""],
            prefix + ["0.01", "Other species", "0.3", "0.01", "0.04", "1", ""],
        ]
        with calls.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(rows)
        result = subprocess.run([
            "python3", str(script), "--calls", str(calls), "--aliases", str(aliases),
            "--spike-panel", str(panel),
            "--outdir", str(outdir),
        ], text=True, capture_output=True)
        assert result.returncode == 0, result.stderr
        with (outdir / "biomarker_propagation_metrics.tsv").open(
                newline="", encoding="utf-8") as handle:
            metrics = list(csv.DictReader(handle, delimiter="\t"))
        assert len(metrics) == 2
        q05 = next(row for row in metrics if float(row["q_threshold"]) == 0.05)
        assert q05["target_called"] == "1"
        assert q05["enriched_calls"] == "2"
        assert q05["off_target_enriched_calls"] == "1"
        assert abs(float(q05["precision"]) - 0.5) < 1e-12
        assert abs(float(q05["target_effect_change_from_baseline"]) - 0.6) < 1e-12
        assert abs(float(q05["biomarker_set_jaccard_vs_baseline"]) - 0.5) < 1e-12
        assert q05["baseline_reference_kind"] == "observed_calls"
        assert (outdir / "SUCCESS").is_file()
    print("[PASS] biomarker-propagation fixture")


if __name__ == "__main__":
    main()
