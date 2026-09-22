import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/plot_fnuc_adenoma_diagnostic.py"


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class FnucDiagnosticTest(unittest.TestCase):
    def test_six_contexts_and_readable_svg(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            matched = root / "matched"
            matched.mkdir()
            models, baselines, spikes, endpoints, operators = [], [], [], [], []
            for cohort in ("feng", "yachida", "zeller"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for contrast in ("CRC_vs_Control", "Adenoma_vs_Control"):
                        models.append(dict(cohort=cohort, profiler=profiler,
                                           feature="Fusobacterium nucleatum", contrast=contrast,
                                           effect="2" if contrast == "CRC_vs_Control" else ".2",
                                           lower_95="1" if contrast == "CRC_vs_Control" else "-.5",
                                           upper_95="3" if contrast == "CRC_vs_Control" else ".9",
                                           q_value=".01" if contrast == "CRC_vs_Control" else ".8"))
                    baselines.append(dict(cohort=cohort, profiler=profiler,
                                          feature="Fusobacterium nucleatum", condition="Adenoma",
                                          positive="1", n="4", positive_abundance_median_percent=".002",
                                          all_sample_iqr_percent=".001"))
                    spikes.append(dict(cohort=cohort, profiler=profiler,
                                       feature="Fusobacterium nucleatum", condition="Adenoma",
                                       dose_percent=".01", n="4", recovery_q1=".8",
                                       recovery_median="1", recovery_q3="1.1",
                                       fraction_below_half=".25"))
                    operators.append(dict(cohort=cohort, profiler=profiler,
                                          holdout_cohort="NONE", analysis_population="independent",
                                          target_label="Fnuc", feature_role="bystander",
                                          feature="Fusobacterium animalis", operator_slope=".2",
                                          eligible_contexts="12"))
                    for sample in range(4):
                        endpoints.append(dict(cohort=cohort, profiler=profiler,
                                              sample_id=f"s{sample}", condition="Adenoma",
                                              analysis_population="community", assembly_arm="original",
                                              target_label="Fnuc", spike_fraction_target=".00001",
                                              reference_type="genome_equivalent" if profiler == "metaphlan4"
                                                             else "read_proportional",
                                              baseline_abundance_fraction="0" if sample else ".001",
                                              observed_abundance_fraction=".00001" if sample else ".00101",
                                              implanted_signal_profiler_scale=".00001",
                                              recovered_spike_signal_profiler_scale=".00001"))
            write(matched / "matched_species_models.tsv", models)
            write(matched / "matched_species_baseline_summary.tsv", baselines)
            write(matched / "matched_species_spike_summary.tsv", spikes)
            write(root / "endpoints.tsv", endpoints)
            write(root / "operator.tsv", operators)
            done = subprocess.run(["python3", str(SCRIPT), "--matched-dir", str(matched),
                                   "--endpoints", str(root / "endpoints.tsv"),
                                   "--operator", str(root / "operator.tsv"),
                                   "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            ET.parse(root / "out/fnuc_adenoma_diagnostic.svg")
            with (root / "out/fnuc_adenoma_diagnostic.tsv").open() as handle:
                result = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(result), 6)
            self.assertEqual({r["community_zero_rescued"] for r in result}, {"3"})
            self.assertTrue((root / "out/SUCCESS").is_file())


if __name__ == "__main__":
    unittest.main()
