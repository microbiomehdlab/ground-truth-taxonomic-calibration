#!/usr/bin/env python3
"""Small end-to-end gate for the three-cohort baseline figure inputs."""
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "analysis_v2/scripts/build_three_cohort_baseline_figure_input.py"
PLOT = ROOT / "analysis_v2/scripts/plot_three_cohort_baseline_discordance.R"


def write(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class BaselineFigureTest(unittest.TestCase):
    def test_three_cohorts_both_profilers_and_zero_filling(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = root / "manifest.tsv"
            abundance = root / "abundance.tsv"
            manifest_rows = []
            abundance_rows = []
            for cohort in ("feng", "yachida", "zeller"):
                for profiler in ("kraken2_bracken", "metaphlan4"):
                    for sample, condition in (("A", "Control"), ("B", "CRC")):
                        source = "/profiles/{}/{}/{}".format(cohort, profiler, sample)
                        for label in ("Bfrag", "Fnuc"):
                            manifest_rows.append(dict(cohort=cohort, sample_id=sample,
                                                      profiler=profiler, condition=condition,
                                                      source_profile=source, dose_level="baseline",
                                                      include="1", target_label=label))
                        abundance_rows.append(dict(profiler=profiler, source_profile=source,
                                                   feature="Bacteroides fragilis",
                                                   abundance_fraction="0.001"))
            write(manifest, list(manifest_rows[0]), manifest_rows)
            write(abundance, list(abundance_rows[0]), abundance_rows)
            out = root / "input"
            result = subprocess.run([
                "python3", str(BUILDER), "--manifest", str(manifest),
                "--abundance", str(abundance), "--aliases",
                str(ROOT / "examples/spike_taxon_aliases.csv"),
                "--panel", str(ROOT / "spikes/spike_panel.tsv"),
                "--outdir", str(out)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (out / "baseline_four_taxa_samples.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3 * 2 * 2 * 4)
            self.assertEqual(sum(float(r["abundance_fraction"]) > 0 for r in rows), 12)
            plot = subprocess.run(["Rscript", str(PLOT), "--input",
                                   str(out / "baseline_four_taxa_samples.tsv"),
                                   "--outdir", str(root / "figures")],
                                  capture_output=True, text=True)
            self.assertEqual(plot.returncode, 0, plot.stdout + plot.stderr)
            self.assertTrue((root / "figures/baseline_three_cohort_combined.pdf").stat().st_size > 0)


if __name__ == "__main__":
    unittest.main()
