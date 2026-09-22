import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from analysis_v2.tests.test_diagnose_crc_candidate_recovery import TargetedDiagnostic

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/plot_crc_candidate_mechanisms.py"


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class MechanismAtlasTest(unittest.TestCase):
    def test_full_baselines_and_paired_spikes_are_separate(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            TargetedDiagnostic().fixture(root)
            manifest, abundance = [], []
            for cohort in ("zeller", "yachida"):
                for condition in ("Control", "CRC"):
                    for profiler in ("kraken2_bracken", "metaphlan4"):
                        feature = "Parvimonas micra" if cohort == "zeller" else (
                            "Allisonella pneumosintes" if profiler == "kraken2_bracken"
                            else "Dialister pneumosintes")
                        for i in range(4):
                            source = str(root / f"{cohort}_{condition}_{profiler}_{i}.tsv")
                            for target in ("Bfrag", "Fnuc"):
                                manifest.append(dict(cohort=cohort, sample_id=f"s{i}",
                                                     condition=condition, analysis_population="community",
                                                     assembly_arm="original", profiler=profiler,
                                                     dose_level="baseline", include="1",
                                                     target_label=target, source_profile=source))
                            if i % 2:
                                abundance.append(dict(profiler=profiler, source_profile=source,
                                                      feature=feature, abundance_fraction=.001))
            write(root / "manifest.tsv", manifest)
            write(root / "abundance.tsv", abundance)
            done = subprocess.run([
                "python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                "--endpoints", str(root / "endpoints.tsv"),
                "--manifest", str(root / "manifest.tsv"),
                "--abundance", str(root / "abundance.tsv"),
                "--aliases", str(root / "aliases.csv"),
                "--outdir", str(root / "out"),
            ], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            for stem in ("zeller_Pmic", "yachida_Dpne"):
                ET.parse(root / "out" / (stem + "_mechanism_atlas.svg"))
            with (root / "out/full_baseline_and_low_dose_summary.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 8)
            self.assertEqual({int(r["n_full_baseline"]) for r in rows}, {4})
            self.assertEqual({float(r["baseline_prevalence"]) for r in rows}, {.5})
            self.assertEqual({int(r["paired_spike_n_at_0p01"]) for r in rows}, {4})
            self.assertTrue((root / "out/SUCCESS").is_file())


if __name__ == "__main__":
    unittest.main()
