import csv
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from plot_matched_crc_adenoma_species import LABELS, matched_calls, matched_features


class MatchedSpeciesTest(unittest.TestCase):
    def test_aliases_include_nonshared_features_in_both_profilers(self):
        aliases = {
            (p, feature): (label, feature)
            for p in ("kraken2_bracken", "metaphlan4")
            for label, feature in (("Fnuc", "Fusobacterium nucleatum"),
                                   ("Pmic", "Parvimonas micra"),
                                   ("Dpne", "Allisonella pneumosintes" if p == "kraken2_bracken"
                                    else "Dialister pneumosintes"))
        }
        selected = matched_features(aliases)
        self.assertEqual(selected["kraken2_bracken"][2], "Allisonella pneumosintes")
        self.assertEqual(selected["metaphlan4"][2], "Dialister pneumosintes")
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "calls.tsv"
            fields = ("cohort", "analysis_population", "assembly_arm", "profiler", "contrast",
                      "include", "model_spec", "spike_fraction_target", "feature", "effect",
                      "lower_95", "upper_95", "q_value", "n_samples")
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for profiler, features in selected.items():
                    for feature in features:
                        for cohort in ("feng", "yachida", "zeller"):
                            for contrast in ("CRC_vs_Control", "Adenoma_vs_Control"):
                                writer.writerow(dict(cohort=cohort, analysis_population="community",
                                    assembly_arm="original", profiler=profiler, contrast=contrast,
                                    include=1, model_spec="primary_age_sex", spike_fraction_target=0,
                                    feature=feature, effect=1, lower_95=0, upper_95=2,
                                    q_value=.8 if feature == "Parvimonas micra" else .01,
                                    n_samples=100))
            rows = matched_calls(path, selected)
        self.assertEqual(len(rows), 2*len(LABELS)*3*2)
        self.assertEqual({r["profiler"] for r in rows if r["feature"] == "Parvimonas micra"},
                         {"kraken2_bracken", "metaphlan4"})
        self.assertEqual({r["significant"] for r in rows if r["feature"] == "Parvimonas micra"}, {0})

    def test_missing_alias_rejected(self):
        with self.assertRaisesRegex(ValueError, "expected one"):
            matched_features({})


if __name__ == "__main__":
    unittest.main()
