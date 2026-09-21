#!/usr/bin/env python3
"""Build three-cohort biomarker-recoverability figure sources on profiler scale."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import statistics
from collections import defaultdict
from pathlib import Path

COHORTS = ("feng", "yachida", "zeller")
CONDITIONS = ("Control", "Adenoma", "CRC")
PROFILERS = ("kraken2_bracken", "metaphlan4")
DOSES = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
FOCUS = 0.0001
DRIVERS = ("baseline_log10", "detectability", "recovery_error", "recovery_iqr")


def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("empty table: {}".format(path))
        yield from reader


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def number(row, field):
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError("nonfinite {}".format(field))
    return value


def condition(value):
    key = value.strip().lower()
    names = {"control": "Control", "healthy": "Control", "adenoma": "Adenoma",
             "adenomas": "Adenoma", "crc": "CRC", "colorectal carcinoma": "CRC"}
    if key not in names:
        raise ValueError("unknown condition: {}".format(value))
    return names[key]


def nominal(value):
    dose = min(DOSES, key=lambda d: abs(value - d) / d)
    if abs(value - dose) / dose > 0.05:
        raise ValueError("independent dose does not match frozen grid: {}".format(value))
    return dose


def tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(biomarkers, endpoints, panel, outdir):
    labels = [row["label"] for row in read(panel)]
    if len(labels) != 10 or len(set(labels)) != 10:
        raise ValueError("expected ten unique frozen panel labels")
    expected_contexts = {(c, g, p, label) for c in COHORTS for g in CONDITIONS
                         for p in PROFILERS for label in labels}
    biomarker = {}
    seen_doses = defaultdict(set)
    for row in read(biomarkers):
        if row["analysis_population"] != "independent" or row["assembly_arm"] != "original":
            continue
        prefix = "spiked_vs_matched_baseline__background_"
        if not row["contrast"].startswith(prefix) or abs(number(row, "q_threshold") - .05) > 1e-12:
            continue
        key = (row["cohort"], condition(row["contrast"][len(prefix):]),
               row["profiler"], row["target_label"])
        if key not in expected_contexts:
            raise ValueError("unexpected biomarker context: {}".format(key))
        dose = nominal(number(row, "spike_fraction_target"))
        identity = key + (dose,)
        if identity in biomarker:
            raise ValueError("duplicate biomarker context/dose: {}".format(identity))
        if row["target_called"] not in {"0", "1"}:
            raise ValueError("invalid target_called")
        q = number(row, "target_q_value")
        if not 0 <= q <= 1:
            raise ValueError("target q outside [0,1]")
        biomarker[identity] = (int(row["target_called"]), q)
        seen_doses[key].add(dose)
    if set(seen_doses) != expected_contexts or any(v != set(DOSES) for v in seen_doses.values()):
        raise ValueError("stratified biomarker table lacks a complete three-cohort six-dose grid")

    observed = defaultdict(dict)
    for row in read(endpoints):
        if row["analysis_population"] != "independent" or row["assembly_arm"] != "original":
            continue
        key = (row["cohort"], condition(row["condition"]), row["profiler"], row["target_label"])
        if key not in expected_contexts:
            raise ValueError("unexpected paired-endpoint context: {}".format(key))
        expected_type = "genome_equivalent" if row["profiler"] == "metaphlan4" else "read_proportional"
        if row["reference_type"] != expected_type:
            raise ValueError("wrong row reference type for {}".format(key))
        dose = nominal(number(row, "spike_fraction_target"))
        identity = key + (dose,)
        sample = row["sample_id"]
        if sample in observed[identity]:
            raise ValueError("duplicate endpoint sample: {} {}".format(identity, sample))
        baseline = number(row, "baseline_abundance_fraction")
        implanted = number(row, "implanted_signal_profiler_scale")
        recovered = number(row, "recovered_spike_signal_profiler_scale")
        expected = number(row, "expected_abundance_profiler_scale")
        actual = number(row, "observed_abundance_fraction")
        if baseline < 0 or implanted <= 0 or expected <= 0 or actual < 0:
            raise ValueError("invalid profiler-scale endpoint")
        detected = row["observed_detected_native_nonzero"]
        if detected not in {"0", "1"} or int(detected) != int(actual > 0):
            raise ValueError("inconsistent observed detection")
        observed[identity][sample] = (baseline, implanted, recovered, expected,
                                      actual, int(detected))
    if {key[:-1] for key in observed} != expected_contexts or any(
            {key[-1] for key in observed if key[:-1] == context} != set(DOSES)
            for context in expected_contexts):
        raise ValueError("paired endpoints lack a complete three-cohort six-dose grid")

    thresholds = []
    drivers = []
    for key in sorted(expected_contexts):
        calls = [dose for dose in DOSES if biomarker[key + (dose,)][0] == 1]
        thresholds.append(dict(cohort=key[0], condition=key[1], profiler=key[2],
                               target_label=key[3], minimum_fraction=(min(calls) if calls else "NR"),
                               tested_doses=len(DOSES)))
        sample_rows = list(observed[key + (FOCUS,)].values())
        if not sample_rows:
            raise ValueError("empty focus-dose endpoint context")
        baseline_log10 = math.log10(statistics.median(r[0] for r in sample_rows) + 1e-6)
        detectability = statistics.mean(r[5] for r in sample_rows)
        recovery_error = statistics.median(abs(r[2] - r[1]) / r[1] for r in sample_rows)
        ratios = sorted(r[2] / r[1] for r in sample_rows)
        # Inclusive quartiles are defined for small independent subsets too.
        q1, _, q3 = statistics.quantiles(ratios, n=4, method="inclusive")
        q = biomarker[key + (FOCUS,)][1]
        values = (baseline_log10, detectability, recovery_error, q3 - q1)
        for name, value in zip(DRIVERS, values):
            drivers.append(dict(cohort=key[0], condition=key[1], profiler=key[2],
                                target_label=key[3], focus_fraction=FOCUS,
                                driver=name, driver_value=format(value, ".17g"),
                                biomarker_q=format(q, ".17g"),
                                biomarker_strength=format(-math.log10(max(q, 1e-300)), ".17g"),
                                samples=len(sample_rows)))
    outdir.mkdir(parents=True, exist_ok=False)
    threshold_path = outdir / "minimum_biomarker_fraction.tsv"
    driver_path = outdir / "biomarker_recovery_drivers.tsv"
    tsv(threshold_path, ("cohort", "condition", "profiler", "target_label",
                         "minimum_fraction", "tested_doses"), thresholds)
    tsv(driver_path, ("cohort", "condition", "profiler", "target_label",
                      "focus_fraction", "driver", "driver_value", "biomarker_q",
                      "biomarker_strength", "samples"), drivers)
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "source.sha256").write_text("".join("{}  {}\n".format(digest(p), p.resolve())
                                               for p in (biomarkers, endpoints, panel,
                                                         threshold_path, driver_path)), encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] {} threshold contexts and {} driver rows".format(len(thresholds), len(drivers)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("biomarkers", "endpoints", "panel", "outdir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.biomarkers, args.endpoints, args.panel, args.outdir)
    except (OSError, ValueError, KeyError, statistics.StatisticsError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error


if __name__ == "__main__":
    main()
