import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/diagnose_crc_candidate_recovery.py"


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TargetedDiagnostic(unittest.TestCase):
    def fixture(self, root, bad_alias=False):
        alias_rows = [
            ("Parvimonas micra", "Parvimonas micra", "kraken2_bracken"),
            ("Parvimonas micra", "Parvimonas micra", "metaphlan4"),
            ("Dialister pneumosintes", "Dialister pneumosintes" if bad_alias else
             "Allisonella pneumosintes", "kraken2_bracken"),
            ("Dialister pneumosintes", "Dialister pneumosintes", "metaphlan4"),
        ]
        with (root / "aliases.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("canonical", "alias", "tool"))
            writer.writerows(alias_rows)
        calls, endpoints = [], []
        for cohort, target, canonical in (("zeller", "Pmic", "Parvimonas micra"),
                                          ("yachida", "Dpne", "Dialister pneumosintes")):
            for profiler in ("kraken2_bracken", "metaphlan4"):
                feature = next(a for c, a, p in alias_rows if c == canonical and p == profiler)
                for label in ("Bfrag", "Fnuc"):
                    calls.append(dict(cohort=cohort, analysis_population="community",
                        assembly_arm="original", profiler=profiler,
                        contrast="CRC_vs_Control", include="1", model_spec="primary_age_sex",
                        spike_fraction_target="0", target_label=label, feature=feature,
                        effect="1.5", standard_error="0.5", lower_95="0.5", upper_95="2.5",
                        q_value="0.01", n_samples="8"))
                for condition in ("Control", "CRC"):
                    for dose in (.0001, .0005, .001):
                        for sample in range(4):
                            baseline = .001 if sample % 2 else 0
                            signal = dose * (2 if profiler == "metaphlan4" else 1)
                            recovered = signal * (0.5 if profiler == "kraken2_bracken" else 1)
                            endpoints.append(dict(cohort=cohort, sample_id="s" + str(sample),
                                condition=condition, analysis_population="independent",
                                assembly_arm="original", profiler=profiler, target_label=target,
                                spike_fraction_target=dose,
                                reference_type="genome_equivalent" if profiler == "metaphlan4"
                                               else "read_proportional",
                                baseline_abundance_fraction=baseline,
                                observed_abundance_fraction=baseline + recovered,
                                implanted_signal_profiler_scale=signal,
                                recovered_spike_signal_profiler_scale=recovered,
                                implanted_read_pairs_target="10"))
        write(root / "calls.tsv", calls)
        write(root / "endpoints.tsv", endpoints)

    def test_paired_recovery_and_alias_are_explicit(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root)
            done = subprocess.run(["python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                "--endpoints", str(root / "endpoints.tsv"), "--aliases", str(root / "aliases.csv"),
                "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            for stem in ("zeller_Pmic", "yachida_Dpne"):
                ET.parse(root / "out" / (stem + "_diagnostic.svg"))
            with (root / "out/baseline_crc_model.tsv").open() as handle:
                model = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(next(r for r in model if r["target_label"] == "Dpne" and
                                  r["profiler"] == "kraken2_bracken")["reported_feature"],
                             "Allisonella pneumosintes")
            with (root / "out/paired_recovery_summary.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 24)
            self.assertEqual({int(r["n_paired"]) for r in rows}, {4})
            self.assertEqual({float(r["recovery_ratio_median"]) for r in rows if
                              r["profiler"] == "kraken2_bracken"}, {0.5})
            self.assertEqual({float(r["recovery_ratio_median"]) for r in rows if
                              r["profiler"] == "metaphlan4"}, {1.0})

    def test_alias_contract_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root, bad_alias=True)
            done = subprocess.run(["python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                "--endpoints", str(root / "endpoints.tsv"), "--aliases", str(root / "aliases.csv"),
                "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("alias contract changed", done.stderr)


if __name__ == "__main__":
    unittest.main()
