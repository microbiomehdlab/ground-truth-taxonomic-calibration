import csv
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis_v2/scripts/plot_three_cohort_crc_spike_impact.py"
COHORTS = ("feng", "yachida", "zeller")
PROFILERS = ("kraken2_bracken", "metaphlan4")


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class SpikeImpact(unittest.TestCase):
    def fixture(self, root, missing=False):
        calls, ledger = [], []
        for profiler in PROFILERS:
            for cohort in COHORTS:
                for feature, effect, q in (("Species A", 1, .01),
                                           ("Species B", -1, .02 if cohort != "zeller" else .4)):
                    for target in ("A", "B"):
                        calls.append(dict(cohort=cohort, analysis_population="community",
                                          assembly_arm="original", profiler=profiler,
                                          contrast="CRC_vs_Control", include="1",
                                          model_spec="primary_age_sex", spike_fraction_target="0",
                                          feature=feature, effect=effect, q_value=q,
                                          n_samples="30", target_label=target))
                    if missing and profiler == "metaphlan4" and cohort == "zeller" and feature == "Species A":
                        continue
                    called = q <= .05
                    is_direct = profiler == "metaphlan4" and cohort == "feng" and feature == "Species A"
                    dose_called = feature == "Species A" or cohort == "feng"
                    ledger.append(dict(cohort=cohort, analysis_population="community",
                                       assembly_arm="original", profiler=profiler,
                                       contrast="CRC_vs_Control", target_label="CRCpanel",
                                       spike_fraction_target="0.001",
                                       spike_fraction_total={"feng": "0.008", "yachida": "0.01", "zeller": "0.012"}[cohort],
                                       q_threshold="0.05", feature=feature,
                                       feature_role="implanted_target" if is_direct else "bystander",
                                       baseline_called=str(int(called)), dose_called=str(int(dose_called)),
                                       baseline_q_value=q,
                                       dose_q_value="0.02" if dose_called else "0.3",
                                       effect_sign_changed="0"))
        write(root / "calls.tsv", calls)
        write(root / "ledger.tsv", ledger)

    def test_both_profilers_all_cohorts_and_spike_fates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root)
            done = subprocess.run(["python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                                   "--ledger", str(root / "ledger.tsv"), "--top", "2",
                                   "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            for profiler in PROFILERS:
                ET.parse(root / "out" / (profiler + "_top_candidates.svg"))
                ET.parse(root / "out" / (profiler + "_shared_candidates.svg"))
            with (root / "out/top_candidate_spike_fates.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 12)
            self.assertEqual({r["profiler"] for r in rows}, set(PROFILERS))
            self.assertEqual({r["cohort"] for r in rows}, set(COHORTS))
            self.assertEqual({r["feature"] for r in rows}, {"Species A", "Species B"})
            self.assertEqual(next(r for r in rows if r["profiler"] == "metaphlan4" and
                                  r["cohort"] == "feng" and r["feature"] == "Species A")["spike_fate"],
                             "direct_target_excluded")
            self.assertEqual(next(r for r in rows if r["profiler"] == "kraken2_bracken" and
                                  r["cohort"] == "yachida" and r["feature"] == "Species B")["spike_fate"],
                             "lost_significance")
            self.assertEqual({float(r["total_mixture_percent"]) for r in rows}, {0.8, 1.0, 1.2})
            with (root / "out/shared_candidate_spike_fates.tsv").open() as handle:
                shared = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(shared), 6)
            self.assertEqual({r["feature"] for r in shared}, {"Species A"})
            self.assertEqual(next(r for r in shared if r["profiler"] == "metaphlan4" and
                                  r["cohort"] == "feng")["post_spike_q"], "")
            self.assertEqual(next(r for r in shared if r["profiler"] == "kraken2_bracken" and
                                  r["cohort"] == "feng")["post_spike_q"], "0.02")

    def test_missing_significant_fate_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.fixture(root, missing=True)
            done = subprocess.run(["python3", str(SCRIPT), "--calls", str(root / "calls.tsv"),
                                   "--ledger", str(root / "ledger.tsv"),
                                   "--outdir", str(root / "out")], capture_output=True, text=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("missing spike fate", done.stderr)


if __name__ == "__main__":
    unittest.main()
