import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import duckdb

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/plot_fnuc_spike_biomarker_audit.py"


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class FnucSpikeBiomarkerAuditTest(unittest.TestCase):
    def test_shared_spike_call_and_response_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger, endpoints = [], []
            for cohort in ("feng", "yachida", "zeller"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for contrast in ("CRC_vs_Control", "Adenoma_vs_Control"):
                        common = dict(cohort=cohort, profiler=profiler, contrast=contrast,
                                      analysis_population="community", assembly_arm="original",
                                      target_label="CRCpanel", spike_fraction_target=".00001",
                                      q_threshold=".05")
                        ledger.append(dict(common, feature="Fusobacterium nucleatum",
                                           feature_role="implanted_target", baseline_called="0",
                                           dose_called="1", baseline_effect=".3", dose_effect="1",
                                           baseline_q_value=".8", dose_q_value=".01", transition="gained"))
                        ledger.append(dict(common, feature="Fusobacterium animalis",
                                           feature_role="bystander", baseline_called="0",
                                           dose_called="1", baseline_effect=".1", dose_effect=".5",
                                           baseline_q_value=".8", dose_q_value=".02", transition="gained"))
                    for condition in ("Control", "Adenoma", "CRC"):
                        for sample in range(3):
                            endpoints.append(dict(cohort=cohort, profiler=profiler, condition=condition,
                                sample_id=f"s{sample}", analysis_population="community",
                                assembly_arm="original", target_label="Fnuc",
                                spike_fraction_target=".00001",
                                reference_type="genome_equivalent" if profiler == "metaphlan4"
                                               else "read_proportional",
                                baseline_abundance_fraction="0", observed_abundance_fraction=".00001",
                                implanted_signal_profiler_scale=".00001",
                                recovered_spike_signal_profiler_scale=".00001"))
            write(root / "ledger.tsv", ledger)
            write(root / "endpoints.tsv", endpoints)
            db = duckdb.connect()
            db.execute("CREATE TABLE response (cohort VARCHAR, profiler VARCHAR, condition VARCHAR, "
                       "analysis_population VARCHAR, dose_index INTEGER, feature VARCHAR, "
                       "is_direct_target BOOLEAN, response_delta DOUBLE)")
            for cohort in ("feng", "yachida", "zeller"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for condition in ("Control", "Adenoma", "CRC"):
                        db.execute("INSERT INTO response VALUES (?, ?, ?, 'community', 1, "
                                   "'Fusobacterium animalis', false, 0.000002)",
                                   [cohort, profiler, condition])
            parquet_path = root / "response.parquet"
            db.execute("COPY response TO '" + str(parquet_path).replace("'", "''") +
                       "' (FORMAT PARQUET)")
            db.close()
            done = subprocess.run(["python3", str(SCRIPT), "--ledger", str(root / "ledger.tsv"),
                                   "--endpoints", str(root / "endpoints.tsv"),
                                   "--paired-features", str(parquet_path),
                                   "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            ET.parse(root / "out/fnuc_spike_biomarker_audit.svg")
            with (root / "out/fnuc_call_transition_summary.tsv").open() as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(summary), 12)
            self.assertEqual({r["target_dose_called"] for r in summary}, {"1"})
            self.assertEqual({r["related_gained"] for r in summary}, {"1"})
            with (root / "out/fnuc_target_condition_response.tsv").open() as handle:
                response = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(response), 18)
            self.assertEqual({r["zero_rescued"] for r in response}, {"3"})
            with (root / "out/fnuc_related_species_response.tsv").open() as handle:
                bystander = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(bystander), 18)
            self.assertTrue((root / "out/SUCCESS").is_file())


if __name__ == "__main__":
    unittest.main()
