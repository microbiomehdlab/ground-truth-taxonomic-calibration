#!/usr/bin/env python3
import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

import duckdb


ROOT = Path(__file__).resolve().parents[2]


class TargetRecoverySummaryTest(unittest.TestCase):
    def test_recovery_classes_and_threshold(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            parquet = tmp / "responses.parquet"
            con = duckdb.connect()
            con.execute(f"""
              COPY (SELECT * FROM (VALUES
                -- Community rows. fnuc_baseline_prevalence.tsv deliberately
                -- restricts to analysis_population='community', so the fixture
                -- must exercise that population; c1 is detected at baseline and
                -- c2 is not, giving exactly one detected sample.
                ('feng','study','c1','CRC','community','kraken2_bracken','cb1',1,
                 0.001,0.001,'Fusobacterium nucleatum',0.002,0.003,0.001,1,1,'read_proportional',0.001,0.001,0.0,0.0,0.0,NULL,0.0,0.0),
                ('feng','study','c2','Control','community','kraken2_bracken','cb2',1,
                 0.001,0.001,'Fusobacterium nucleatum',0.0,0.001,0.001,1,1,'read_proportional',0.001,0.001,0.0,0.0,0.0,NULL,0.0,0.0),
                ('feng','study','s1','CRC','independent','kraken2_bracken','b1',1,
                 0.001,0.001,'Fusobacterium nucleatum',0.0,0.001,0.001,1,1,'read_proportional',0.001,0.001,0.0,0.0,0.0,NULL,0.0,0.0),
                ('feng','study','s2','Control','independent','kraken2_bracken','b2',1,
                 0.001,0.001,'Fusobacterium nucleatum',0.002,0.0022,0.0002,1,1,'read_proportional',0.0002,0.001,0.0,0.0,0.0,NULL,0.0,0.0),
                ('feng','study','s3','CRC','independent','kraken2_bracken','b3',1,
                 0.001,0.001,'Bacteroides fragilis',0.0,0.0,0.0,0,1,'read_proportional',0.0,0.001,0.0,0.0,0.0,NULL,0.0,0.0)
              ) AS t(cohort,study,sample_id,condition,analysis_population,profiler,
                baseline_id,dose_index,nominal_target_fraction,target_fraction_for_feature,
                feature,baseline_abundance_fraction,observed_abundance_fraction,
                response_signal,observed_detected,is_direct_target,reference_type,
                response_signal_profiler_scale,implanted_signal_profiler_scale,
                expected_abundance_profiler_scale,
                retained_baseline_profiler_scale,response_delta_profiler_scale,
                quantitative_log2_ratio_profiler_scale,
                signed_bounded_error_profiler_scale,
                absolute_bounded_error_profiler_scale))
              TO '{str(parquet).replace("'", "''")}' (FORMAT PARQUET)
            """)
            metrics = tmp / "metrics.tsv"
            metrics.write_text(
                "cohort\tanalysis_population\ttarget_label\tprofiler\tcontrast\t"
                "spike_fraction_target\tq_threshold\ttarget_called\ttarget_effect\ttarget_q_value\n"
                "feng\tindependent\tFnuc\tkraken2_bracken\t"
                "spiked_vs_matched_baseline__pooled_conditions\t0.001\t0.05\t1\t2\t0.01\n",
                encoding="utf-8",
            )
            out = tmp / "out"
            subprocess.run([
                "python3", str(ROOT / "analysis_v2/scripts/summarize_target_recovery.py"),
                "--responses", str(parquet), "--biomarker-metrics", str(metrics),
                "--outdir", str(out), "--reference-scale", "profiler_scale",
            ], check=True)
            with (out / "target_recovery_observations.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            # Ordered by cohort, analysis_population, profiler, feature,
            # sample_id: the two community Fnuc rows sort before the
            # independent ones.
            self.assertEqual([row["recovery_class"] for row in rows],
                             ["Good", "Good",
                              "Poor / missed", "Good", "Poor / missed"])
            with (out / "fnuc_baseline_prevalence.tsv").open(newline="") as handle:
                prevalence = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(prevalence, "community population must reach prevalence output")
            self.assertTrue(all(row["analysis_population"] == "community"
                                for row in prevalence))
            self.assertEqual(sum(int(row["detected_samples"]) for row in prevalence), 1)
            self.assertTrue((out / "biomarker_detection_thresholds.tsv").is_file())
            self.assertTrue((out / "SUCCESS").is_file())


if __name__ == "__main__":
    unittest.main()
