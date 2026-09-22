import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/plot_three_cohort_metaphlan_crc_candidates.py"


class CandidatePlot(unittest.TestCase):
    def fixture(self, path, contradiction=False):
        rows = []
        values = {
            "A": ((1.2, .01), (.8, .03), (1.0, .02)),
            "B": ((-1.1, .01), (-.4, .2), (-.9, .04)),
            "C": ((.6, .03), (-.5, .01), (.1, .8)),
            "D": ((.1, .8), (.2, .7), (.3, .6)),
        }
        for feature, fits in values.items():
            for cohort, (effect, q) in zip(("feng", "yachida", "zeller"), fits):
                for target in ("Bfrag", "Fnuc"):
                    rows.append(dict(cohort=cohort, analysis_population="community",
                                     assembly_arm="original", profiler="metaphlan4",
                                     contrast="CRC_vs_Control", spike_fraction_target="0",
                                     include="1", model_spec="primary_age_sex", feature=feature,
                                     effect=effect + (1 if contradiction and feature == "A" and
                                                     cohort == "feng" and target == "Fnuc" else 0),
                                     q_value=q, n_samples="50"))
        rows.append(dict(rows[0], analysis_population="independent", feature="Ignore", q_value="0.001"))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_shared_candidates_and_all_cohort_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "calls.tsv")
            result = subprocess.run(["python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                                     "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            ET.parse(root / "out/metaphlan_crc_candidates.svg")
            with (root / "out/metaphlan_crc_candidates_all.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({r["feature"] for r in rows}, {"A", "B", "C"})
            self.assertEqual(len(rows), 9)
            self.assertEqual({r["same_direction_significant_cohorts"] for r in rows if r["feature"] == "A"}, {"3"})
            self.assertEqual({r["same_direction_significant_cohorts"] for r in rows if r["feature"] == "B"}, {"2"})
            self.assertEqual({r["same_direction_significant_cohorts"] for r in rows if r["feature"] == "C"}, {"1"})

    def test_contradictory_repeated_baseline_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root / "calls.tsv", contradiction=True)
            result = subprocess.run(["python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                                     "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("contradictory baseline fit", result.stderr)


if __name__ == "__main__":
    unittest.main()
