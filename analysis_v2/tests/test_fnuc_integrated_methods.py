import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/plot_fnuc_integrated_methods.py"


def write(path, records):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


class IntegratedMethodsTest(unittest.TestCase):
    def test_plot_and_missing_context_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            matched, audit, high = root / "matched", root / "audit", root / "high"
            matched.mkdir()
            audit.mkdir()
            high.mkdir()
            models, baselines, spikes, transitions, response = [], [], [], [], []
            for profiler in ("kraken2_bracken", "metaphlan4"):
                for cohort in ("feng", "yachida", "zeller"):
                    common = dict(profiler=profiler, cohort=cohort)
                    for contrast in ("CRC_vs_Control", "Adenoma_vs_Control"):
                        models.append(dict(common, feature="Fusobacterium nucleatum",
                                           contrast=contrast, effect="1", lower_95="0",
                                           upper_95="2", q_value=".4"))
                        transitions.append(dict(common, contrast=contrast,
                            target_baseline_effect="1", target_dose_effect=".5",
                            target_baseline_q=".4", target_dose_q=".8",
                            target_baseline_called="0", target_dose_called="0",
                            related_gained="0", related_lost="0",
                            other_bystander_gained="0", other_bystander_lost="1"))
                    for condition in ("Control", "Adenoma", "CRC"):
                        baselines.append(dict(common, feature="Fusobacterium nucleatum",
                            condition=condition, positive="2", n="10",
                            positive_abundance_median_percent=".002",
                            all_sample_iqr_percent=".001"))
                        spikes.append(dict(common, feature="Fusobacterium nucleatum",
                            condition=condition, dose_percent=".01", observed_positive="9",
                            n="10", recovery_median="1", recovery_q1=".8", recovery_q3="1.2"))
                        response.append(dict(common, condition=condition, baseline_zero="8",
                            zero_rescued="7", n="10", recovery_median=".9",
                            recovery_q1=".7", recovery_q3="1.1", below_half="1"))
            write(matched / "matched_species_models.tsv", models)
            write(matched / "matched_species_baseline_summary.tsv", baselines)
            write(matched / "matched_species_spike_summary.tsv", spikes)
            write(audit / "fnuc_call_transition_summary.tsv", transitions)
            write(audit / "fnuc_target_condition_response.tsv", response)
            write(high / "fnuc_call_transition_summary.tsv", transitions)
            write(high / "fnuc_target_condition_response.tsv", response)
            (high / "dose_metadata.tsv").write_text(
                "member_dose_percent\tapprox_total_mix_percent\n0.1\t1\n", encoding="utf-8")
            (audit / "fnuc_related_species_response.tsv").write_text(
                "cohort\tprofiler\tcondition\tfeature\tresponse_median\n", encoding="utf-8")
            command = ["python3", str(SCRIPT), "--matched-dir", str(matched),
                       "--audit-dir", str(audit), "--high-audit-dir", str(high),
                       "--outdir", str(root / "out")]
            done = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            ET.parse(root / "out/fnuc_integrated_methods.svg")
            figure = (root / "out/fnuc_integrated_methods.svg").read_text(encoding="utf-8")
            self.assertIn("0.001% adenoma effect", figure)
            self.assertIn("0.1% adenoma call", figure)
            self.assertIn("0.1% other calls", figure)
            self.assertNotIn("paired response not supplied", figure)
            self.assertNotIn("0.01% independent recovery", figure)
            self.assertTrue((root / "out/SUCCESS").is_file())
            write(matched / "matched_species_models.tsv", models[:-1])
            failed = subprocess.run(command[:-1] + [str(root / "bad")], capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertFalse((root / "bad/SUCCESS").exists())


if __name__ == "__main__":
    unittest.main()
