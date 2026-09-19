#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import duckdb


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/build_perturbation_response_input.py"


def write(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ResponseInputTest(unittest.TestCase):
    def test_community_is_one_physical_observation_with_vector_expectation(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            manifest = root / "manifest.tsv"
            endpoints = root / "endpoints.tsv"
            abundance = root / "abundance.tsv"
            aliases = root / "aliases.tsv"
            out = root / "out"
            targets = [("A", "Species A", 0.03), ("B", "Species B", 0.02)]
            m_fields = [
                "cohort", "study", "analysis_population", "sample_id", "condition",
                "target_label", "assembly_arm", "profiler", "profile_id",
                "baseline_profile_id", "spike_fraction_target", "dose_level",
                "source_profile", "target_feature",
            ]
            e_fields = [
                "cohort", "study", "analysis_population", "sample_id", "condition",
                "target_label", "assembly_arm", "profiler", "profile_id",
                "baseline_profile_id", "spike_fraction_total", "spike_fraction_target",
                "source_baseline_profile", "source_profile",
            ]
            m_rows = []
            e_rows = []
            for label, feature, fraction in targets:
                m_rows.append({
                    "cohort": "yachida", "study": "S", "analysis_population": "community",
                    "sample_id": "X", "condition": "Control", "target_label": label,
                    "assembly_arm": "original", "profiler": "kraken2_bracken",
                    "profile_id": "X_mix", "baseline_profile_id": "X",
                    "spike_fraction_target": fraction, "dose_level": "dose_06",
                    "source_profile": "/profiles/mix", "target_feature": feature,
                })
                e_rows.append({
                    "cohort": "yachida", "study": "S", "analysis_population": "community",
                    "sample_id": "X", "condition": "Control", "target_label": label,
                    "assembly_arm": "original", "profiler": "kraken2_bracken",
                    "profile_id": "X_mix", "baseline_profile_id": "X",
                    "spike_fraction_total": 0.05, "spike_fraction_target": fraction,
                    "source_baseline_profile": "/profiles/base",
                    "source_profile": "/profiles/mix",
                })
            write(manifest, m_fields, m_rows)
            write(endpoints, e_fields, e_rows)
            write(abundance, ["profiler", "source_profile", "feature", "abundance_fraction"], [
                {"profiler": "kraken2_bracken", "source_profile": "/profiles/base",
                 "feature": "Bystander", "abundance_fraction": 0.10},
                {"profiler": "kraken2_bracken", "source_profile": "/profiles/mix",
                 "feature": "Bystander", "abundance_fraction": 0.095},
                {"profiler": "kraken2_bracken", "source_profile": "/profiles/mix",
                 "feature": "Species A", "abundance_fraction": 0.02},
                {"profiler": "kraken2_bracken", "source_profile": "/profiles/mix",
                 "feature": "Species A_A", "abundance_fraction": 0.01},
            ])
            write(aliases, ["profiler", "source_feature", "canonical_feature", "rationale"], [{
                "profiler": "kraken2_bracken", "source_feature": "Species A_A",
                "canonical_feature": "Species A", "rationale": "fixture split label",
            }])
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "--profile-manifest", str(manifest),
                 "--endpoints", str(endpoints), "--abundance", str(abundance),
                 "--feature-aliases", str(aliases),
                 "--outdir", str(out), "--threads", "1", "--memory-limit", "1GB",
                 "--expected-community-targets", "2"],
                text=True, capture_output=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            con = duckdb.connect()
            obs = con.execute(
                "SELECT implanted_target_count,effective_total_fraction FROM "
                f"read_csv('{out / 'observations.tsv'}', delim='\\t', header=true)"
            ).fetchall()
            self.assertEqual(obs, [(2, 0.05)])
            rows = con.execute(
                "SELECT feature,expected_abundance_fraction,observed_abundance_fraction,"
                "is_direct_target FROM read_parquet(?) ORDER BY feature",
                [str(out / "paired_feature_responses.parquet")],
            ).fetchall()
            self.assertEqual([row[0] for row in rows], ["Bystander", "Species A", "Species B"])
            by_feature = {row[0]: row[1:] for row in rows}
            self.assertAlmostEqual(by_feature["Bystander"][0], 0.095)
            self.assertAlmostEqual(by_feature["Species A"][0], 0.03)
            self.assertAlmostEqual(by_feature["Species B"][0], 0.02)
            self.assertEqual(by_feature["Species B"][1:], (0.0, 1))
            self.assertTrue((out / "SUCCESS").is_file())


if __name__ == "__main__":
    unittest.main()
