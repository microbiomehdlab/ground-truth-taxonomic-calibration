import csv
import importlib.util
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/bridge_crc_replication_reliability.py"
ANALYZER = Path(__file__).resolve().parents[1] / "scripts/analyze_perturbation_response.py"
spec = importlib.util.spec_from_file_location("bridge", SCRIPT)
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)


def write(path, fields, rows, delimiter="\t"):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


class BridgeTest(unittest.TestCase):
    def test_existing_analyzer_accepts_spike_panel_alias_csv(self):
        analyzer_spec = importlib.util.spec_from_file_location("response_analyzer", ANALYZER)
        analyzer = importlib.util.module_from_spec(analyzer_spec)
        import sys
        sys.modules[analyzer_spec.name] = analyzer
        analyzer_spec.loader.exec_module(analyzer)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            disease = [dict(cohort="feng", analysis_population="community",
                            profiler="kraken2_bracken", feature="Alias A", effect=2,
                            q_value=.01, include=1, dose_level="baseline",
                            contrast="CRC_vs_Control"),
                       dict(cohort="yachida", analysis_population="community",
                            profiler="kraken2_bracken", feature="Alias A", effect=1,
                            q_value=.02, include=1, dose_level="baseline",
                            contrast="CRC_vs_Control")]
            write(root / "disease.tsv", tuple(disease[0]), disease)
            cert = [dict(holdout_cohort="feng", training_cohorts="yachida",
                         profiler="kraken2_bracken", feature="Species A",
                         overall_reliability=.8)]
            write(root / "cert.tsv", tuple(cert[0]), cert)
            write(root / "alias.csv", ("canonical", "alias", "tool"),
                  [dict(canonical="Species A", alias="Alias A", tool="kraken2_bracken")], ",")
            rows = analyzer.build_biomarker_replication(root / "disease.tsv", root / "cert.tsv",
                                                       {"feng", "yachida"}, root / "alias.csv")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["feature"], "Species A")

    def test_source_destination_and_alias(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            calls = []
            for cohort, effect, q in (("feng", 2, .01), ("yachida", 1, .2), ("zeller", -1, .01)):
                for target in ("A", "B"):
                    calls.append(dict(cohort=cohort, analysis_population="community",
                                      assembly_arm="original", profiler="kraken2_bracken",
                                      contrast="CRC_vs_Control", include="1",
                                      model_spec="primary_age_sex", spike_fraction_target="0",
                                      feature="Alias A", effect=effect, q_value=q,
                                      n_samples=50, target_label=target))
            write(root / "calls.tsv", tuple(calls[0]), calls)
            certs = [dict(holdout_cohort=cohort, profiler="kraken2_bracken", feature="Alias A",
                          overall_reliability=.8, eligible_contexts=100,
                          evaluation_overall_reliability=score, evaluation_contexts=20)
                     for cohort, score in (("feng", .9), ("yachida", .4), ("zeller", .2))]
            write(root / "cert.tsv", tuple(certs[0]), certs)
            write(root / "alias.csv", ("canonical", "alias", "tool"),
                  [dict(canonical="Species A", alias="Alias A", tool="kraken2_bracken")], ",")
            rows = bridge.make_rows(root / "calls.tsv", root / "cert.tsv", root / "alias.csv")
            self.assertEqual(len(rows), 4)
            feng_to_yachida = next(r for r in rows if r["source_cohort"] == "feng"
                                   and r["destination_cohort"] == "yachida")
            self.assertEqual(feng_to_yachida["feature"], "Species A")
            self.assertEqual(feng_to_yachida["replication_status"], "SAME_DIRECTION_NOT_SIGNIFICANT")
            self.assertEqual(feng_to_yachida["destination_reliability"], .4)
            feng_to_zeller = next(r for r in rows if r["source_cohort"] == "feng"
                                  and r["destination_cohort"] == "zeller")
            self.assertEqual(feng_to_zeller["replication_status"], "DIRECTION_REVERSED")
            self.assertEqual(feng_to_zeller["replicated"], 0)
            summary = bridge.make_summary(rows)
            bridge.plot(summary, root / "figure.svg")
            self.assertGreater((root / "figure.svg").stat().st_size, 1000)
            ET.parse(root / "figure.svg")

    def test_missing_destination_is_not_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            calls = [dict(cohort="feng", analysis_population="community", assembly_arm="original",
                          profiler="metaphlan4", contrast="CRC_vs_Control", include="1",
                          model_spec="primary_age_sex", spike_fraction_target="0",
                          feature="Species B", effect=2, q_value=.01, n_samples=40)]
            write(root / "calls.tsv", tuple(calls[0]), calls)
            certs = [dict(holdout_cohort="yachida", profiler="metaphlan4", feature="Species B",
                          overall_reliability=.7, eligible_contexts=50,
                          evaluation_overall_reliability="", evaluation_contexts=0)]
            write(root / "cert.tsv", tuple(certs[0]), certs)
            write(root / "alias.csv", ("canonical", "alias", "tool"), [], ",")
            rows = bridge.make_rows(root / "calls.tsv", root / "cert.tsv", root / "alias.csv")
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(r["replication_status"] == "NOT_EVALUABLE" for r in rows))
            self.assertEqual(bridge.make_summary(rows)[3]["n"], 0)


if __name__ == "__main__":
    unittest.main()
