#!/usr/bin/env python3
"""Make four-taxon baseline figure data from a sealed three-cohort native input."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path

LABELS = ("Bfrag", "Fnuc", "Pint", "Pmic")
COHORTS = {"feng", "yachida", "zeller"}
PROFILERS = {"kraken2_bracken", "metaphlan4"}
FIELDS = ("cohort", "sample_id", "condition", "profiler", "source_profile",
          "target_label", "abundance_fraction")


def table(path, delimiter="\t"):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("empty table: {}".format(path))
        yield from reader


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build(manifest, abundance, aliases, panel, outdir):
    taxa = {row["label"]: row["taxon_name"] for row in table(panel)}
    if not set(LABELS) <= taxa.keys():
        raise ValueError("spike panel lacks a requested taxon")
    alias_map = {(row["canonical"], row["tool"]): row["alias"]
                 for row in table(aliases, ",")}
    feature_to_label = {}
    for label in LABELS:
        for profiler in PROFILERS:
            feature = alias_map.get((taxa[label], profiler))
            if not feature:
                raise ValueError("missing alias for {} {}".format(label, profiler))
            key = (profiler, feature)
            if key in feature_to_label:
                raise ValueError("two selected taxa share profiler feature {}".format(key))
            feature_to_label[key] = label

    samples = {}
    profile_owners = defaultdict(set)
    for row in table(manifest):
        if row["dose_level"] != "baseline" or row["include"] != "1":
            continue
        key = (row["cohort"], row["sample_id"], row["profiler"])
        value = (row["condition"], str(Path(row["source_profile"]).resolve()))
        previous = samples.setdefault(key, value)
        if previous != value:
            raise ValueError("conflicting baseline for {}".format(key))
        profile_owners[(row["profiler"], value[1])].add(key)
    if {key[0] for key in samples} != COHORTS:
        raise ValueError("baseline input must contain exactly Feng, Yachida, Zeller")
    if {key[2] for key in samples} != PROFILERS:
        raise ValueError("baseline input must contain both profilers")
    if {(key[0], key[2]) for key in samples} != {
            (cohort, profiler) for cohort in COHORTS for profiler in PROFILERS}:
        raise ValueError("a cohort-profiler baseline stratum is missing")
    if not samples:
        raise ValueError("no included baseline profiles")
    abundance_by_sample = {}
    for row in table(abundance):
        label = feature_to_label.get((row["profiler"], row["feature"]))
        if label is None:
            continue
        value = float(row["abundance_fraction"])
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("invalid native abundance")
        resolved_source = str(Path(row["source_profile"]).resolve())
        for sample_key in profile_owners.get((row["profiler"], resolved_source), ()):
            key = sample_key + (label,)
            if key in abundance_by_sample:
                raise ValueError("duplicate selected feature for {}".format(key))
            abundance_by_sample[key] = value

    outdir.mkdir(parents=True, exist_ok=False)
    output = outdir / "baseline_four_taxa_samples.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (cohort, sample_id, profiler), (condition, source) in sorted(samples.items()):
            for label in LABELS:
                writer.writerow(dict(cohort=cohort, sample_id=sample_id,
                                     condition=condition, profiler=profiler,
                                     source_profile=source, target_label=label,
                                     abundance_fraction=format(abundance_by_sample.get(
                                         (cohort, sample_id, profiler, label), 0.0), ".17g")))
    counts = Counter((key[0], key[2]) for key in samples)
    summary = outdir / "baseline_four_taxa_summary.tsv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("cohort", "profiler", "baseline_samples"))
        for (cohort, profiler), count in sorted(counts.items()):
            writer.writerow((cohort, profiler, count))
    checksums = outdir / "baseline_four_taxa.sha256"
    checksums.write_text("".join("{}  {}\n".format(sha256(p), p.resolve())
                                 for p in (manifest, abundance, aliases, panel, output, summary)),
                         encoding="utf-8")
    (outdir / "DEVELOPMENT_ONLY.txt").write_text("status\tDEVELOPMENT_ONLY\n", encoding="utf-8")
    (outdir / "SUCCESS").write_text("status\tPASS\n", encoding="utf-8")
    print("[PASS] {} baseline sample-profiler profiles; {} figure rows".format(
        len(samples), len(samples) * len(LABELS)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("manifest", "abundance", "aliases", "panel", "outdir"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.manifest, args.abundance, args.aliases, args.panel, args.outdir)
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error


if __name__ == "__main__":
    main()
