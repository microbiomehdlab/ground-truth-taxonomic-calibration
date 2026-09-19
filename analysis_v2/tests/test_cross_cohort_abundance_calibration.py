#!/usr/bin/env python3
import csv
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
M = ["cohort", "study", "analysis_population", "sample_id", "condition", "target_label",
     "assembly_arm", "profiler", "profile_id", "baseline_profile_id", "spike_fraction_target",
     "dose_level", "source_profile", "target_taxon", "target_feature", "age", "sex", "bmi",
     "include", "exclusion_reason"]
E = ["cohort", "study", "sample_id", "condition", "analysis_population", "target_label",
     "assembly_arm", "profiler", "profile_id", "baseline_profile_id", "spike_fraction_total",
     "spike_fraction_target", "source_baseline_profile", "source_profile"]
A = ["profiler", "source_profile", "feature", "abundance_fraction"]

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp); manifest = tmp / "manifest.tsv"; endpoints = tmp / "endpoints.tsv"
    abundance = tmp / "abundance.tsv"; out = tmp / "out"
    mr, er, ar = [], [], []
    for cohort in ("feng", "zeller"):
        for i in (1, 2):
            sample = f"{cohort}{i}"; base = str(tmp / f"{sample}.base"); spike = str(tmp / f"{sample}.spike")
            common = [cohort, cohort + "_study", "independent", sample, "CRC" if i == 1 else "Control",
                      "Target", "original", "metaphlan4"]
            mr.append(common + [sample + ".base", sample + ".base", "0", "baseline", base,
                                "Target species", "Target species", "60", "Female", "25", "1", ""])
            mr.append(common + [sample + ".spike", sample + ".base", "0.1", "dose_01", spike,
                                "Target species", "Target species", "60", "Female", "25", "1", ""])
            er.append([cohort, cohort + "_study", sample, common[4], "independent", "Target",
                       "original", "metaphlan4", sample + ".spike", sample + ".base",
                       "0.1", "0.1", base, spike])
            ar.extend([["metaphlan4", base, "Bystander", "0.2"],
                       ["metaphlan4", base, "Target species", "0.05"],
                       ["metaphlan4", spike, "Bystander", "0.19"], # expected .18 + .01 inflation
                       ["metaphlan4", spike, "Target species", "0.145"]])
    for path, fields, rows in ((manifest, M, mr), (endpoints, E, er), (abundance, A, ar)):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields); writer.writerows(rows)
    subprocess.run(["python3", str(ROOT / "analysis_v2/scripts/calibrate_abundance_cross_cohort.py"),
                    "--profile-manifest", str(manifest), "--abundance-long", str(abundance),
                    "--paired-endpoints", str(endpoints), "--outdir", str(out),
                    "--min-training-samples", "2", "--min-positive-residual-fraction", "0.5"], check=True)
    transfer = out / "transfers/feng__to__zeller"
    with (transfer / "calibration_coefficients.tsv").open(newline="", encoding="utf-8") as handle:
        coef = [r for r in csv.DictReader(handle, delimiter="\t") if r["feature"] == "Bystander"][0]
    assert coef["correction_active"] == "1"
    assert abs(float(coef["slope_fraction_per_implanted_fraction"]) - .1) < 1e-10
    with (transfer / "corrected_abundance_long.tsv").open(newline="", encoding="utf-8") as handle:
        corrected = [r for r in csv.DictReader(handle, delimiter="\t")
                     if r["feature"] == "Bystander" and "virtual_profiles" in r["source_profile"]]
    assert corrected and all(abs(float(r["abundance_fraction"]) - .18) < 1e-10 for r in corrected)
    with (transfer / "corrected_abundance_long.tsv").open(newline="", encoding="utf-8") as handle:
        protected = [r for r in csv.DictReader(handle, delimiter="\t")
                     if r["feature"] == "Target species" and "virtual_profiles" in r["source_profile"]]
    assert protected and all(abs(float(r["abundance_fraction"]) - .145) < 1e-10 for r in protected)
    assert (transfer / "SUCCESS").is_file() and (out / "SUCCESS").is_file()
print("[PASS] cross-cohort abundance-calibration fixture")
