import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/plot_shared_crc_adenoma_bridge.py"
COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")
CONDITIONS = ("Control", "Adenoma", "CRC")


def write(path, rows, delimiter="\t"):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


class SharedCRCAdenomaTest(unittest.TestCase):
    def test_crc_anchor_keeps_nonsignificant_adenoma_and_marks_unspiked_species(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            calls, manifest, abundance, endpoints = [], [], [], []
            for cohort in COHORTS:
                for profiler in PROFILERS:
                    for feature in ("Species A", "Species B"):
                        for contrast in ("CRC_vs_Control", "Adenoma_vs_Control"):
                            for target in ("A", "B"):
                                calls.append(dict(cohort=cohort, analysis_population="community",
                                    assembly_arm="original", profiler=profiler, contrast=contrast,
                                    include="1", model_spec="primary_age_sex",
                                    spike_fraction_target="0", feature=feature,
                                    effect="2" if contrast=="CRC_vs_Control" else "0.2",
                                    lower_95="1" if contrast=="CRC_vs_Control" else "-1",
                                    upper_95="3" if contrast=="CRC_vs_Control" else "1",
                                    q_value="0.01" if contrast=="CRC_vs_Control" else "0.8",
                                    n_samples="12", target_label=target))
                    for condition in CONDITIONS:
                        for sample in range(4):
                            source = str(root/f"{cohort}_{profiler}_{condition}_{sample}.tsv")
                            for target in ("A", "B"):
                                manifest.append(dict(cohort=cohort, sample_id=f"s{sample}",
                                    condition=condition, analysis_population="community",
                                    assembly_arm="original", profiler=profiler, dose_level="baseline",
                                    include="1", source_profile=source, target_label=target))
                            abundance.append(dict(profiler=profiler, source_profile=source,
                                                  feature="Species A", abundance_fraction=.001))
                            if sample % 2:
                                abundance.append(dict(profiler=profiler, source_profile=source,
                                                      feature="Species B", abundance_fraction=.002))
                        for dose in (.0001, .0005, .001):
                            for sample in range(3):
                                endpoints.append(dict(cohort=cohort, sample_id=f"s{sample}",
                                    condition=condition, analysis_population="independent",
                                    assembly_arm="original", profiler=profiler, target_label="A",
                                    spike_fraction_target=dose,
                                    reference_type="genome_equivalent" if profiler=="metaphlan4"
                                                   else "read_proportional",
                                    baseline_abundance_fraction=.001, observed_abundance_fraction=.001+dose,
                                    implanted_signal_profiler_scale=dose,
                                    recovered_spike_signal_profiler_scale=dose))
            write(root/"calls.tsv", calls)
            write(root/"manifest.tsv", manifest)
            write(root/"abundance.tsv", abundance)
            write(root/"endpoints.tsv", endpoints)
            write(root/"panel.tsv", [dict(label="A", taxon_name="Species A")])
            write(root/"aliases.csv", [dict(canonical="Species A", alias="Species A",
                                            tool=profiler) for profiler in PROFILERS], ",")
            cmd = ["python3", str(SCRIPT), "--calls", str(root/"calls.tsv"),
                   "--manifest", str(root/"manifest.tsv"),
                   "--abundance", str(root/"abundance.tsv"),
                   "--endpoints", str(root/"endpoints.tsv"),
                   "--panel", str(root/"panel.tsv"),
                   "--aliases", str(root/"aliases.csv"), "--outdir", str(root/"out")]
            done = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout+done.stderr)
            for profiler in PROFILERS:
                ET.parse(root/"out"/(profiler+"_crc_anchored_adenoma.svg"))
                ET.parse(root/"out"/(profiler+"_adenoma_measurement.svg"))
            with (root/"out/shared_crc_adenoma_model.tsv").open() as handle:
                model = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(model), 24)
            self.assertEqual({r["feature"] for r in model}, {"Species A", "Species B"})
            self.assertEqual({r["significant"] for r in model if r["contrast"]=="Adenoma_vs_Control"}, {"0"})
            self.assertEqual({r["direct_spike_target"] for r in model if r["feature"]=="Species B"}, {""})
            with (root/"out/shared_crc_baseline_summary.tsv").open() as handle:
                baseline = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(baseline), 36)
            self.assertEqual({r["n"] for r in baseline}, {"4"})
            with (root/"out/shared_crc_direct_spike_summary.tsv").open() as handle:
                spike = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(spike), 54)
            self.assertEqual({r["feature"] for r in spike}, {"Species A"})
            self.assertTrue((root/"out/SUCCESS").is_file())


if __name__ == "__main__":
    unittest.main()
