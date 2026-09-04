#!/usr/bin/env python3
"""Evaluate target and off-target propagation from standardized DA results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path
from typing import Optional


KEY = ["cohort", "study", "analysis_population", "target_label",
       "assembly_arm", "profiler", "contrast"]
REQUIRED = KEY + ["spike_fraction_target", "feature", "effect", "p_value",
                  "q_value", "include", "exclusion_reason"]
OUTPUT = KEY + ["spike_fraction_target", "q_threshold", "target_alias",
                "target_called", "enriched_calls", "off_target_enriched_calls",
                "precision", "recall", "f1", "target_effect", "target_q_value",
                "target_effect_change_from_baseline", "biomarker_set_jaccard_vs_baseline"]


def read(path: Path, delimiter: str = "\t",
         required: Optional[list[str]] = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError("empty table: {}".format(path))
        missing = set(required or REQUIRED) - set(reader.fieldnames)
        if missing:
            raise ValueError("{} missing columns: {}".format(path, ", ".join(sorted(missing))))
        return list(reader)


def number(row: dict[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except ValueError as error:
        raise ValueError("invalid {} for feature {}".format(field, row.get("feature", "?"))) from error
    if not math.isfinite(value):
        raise ValueError("non-finite {} for feature {}".format(field, row.get("feature", "?")))
    return value


def render(value: float) -> str:
    return format(value, ".17g")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--aliases", type=Path, required=True)
    parser.add_argument("--spike-panel", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--q-thresholds", default="0.05,0.10")
    args = parser.parse_args()
    try:
        rows = read(args.calls)
        aliases = read(args.aliases, delimiter=",",
                       required=["canonical", "alias", "tool"])
        panel = read(args.spike_panel, required=["label", "taxon_name"])
        thresholds = [float(value) for value in args.q_thresholds.split(",")]
        if not thresholds or any(not 0 < value < 1 for value in thresholds) or len(set(thresholds)) != len(thresholds):
            raise ValueError("q thresholds must be distinct values between zero and one")
        taxa = {row["label"]: row["taxon_name"] for row in panel}
        alias_by_taxon = {(row["canonical"], row["tool"]): row["alias"]
                          for row in aliases}
        alias_map = {}
        for label, taxon in taxa.items():
            for profiler in ("kraken2_bracken", "metaphlan4"):
                if (taxon, profiler) in alias_by_taxon:
                    alias_map[(label, profiler)] = alias_by_taxon[(taxon, profiler)]
        included = []
        excluded = []
        identities = set()
        for row in rows:
            if row["include"] not in {"0", "1"}:
                raise ValueError("include must be 0 or 1")
            if (row["include"] == "0") != bool(row["exclusion_reason"].strip()):
                raise ValueError("include/exclusion_reason mismatch")
            dose = number(row, "spike_fraction_target")
            effect = number(row, "effect"); p_value = number(row, "p_value"); q_value = number(row, "q_value")
            if dose < 0 or not 0 <= p_value <= 1 or not 0 <= q_value <= 1:
                raise ValueError("dose or probability outside valid range")
            identity = tuple(row[field] for field in KEY) + (render(dose), row["feature"])
            if identity in identities:
                raise ValueError("duplicate analytical feature row: {}".format(identity))
            identities.add(identity)
            row = dict(row, _dose=dose, _effect=effect, _q=q_value)
            (included if row["include"] == "1" else excluded).append(row)

        groups = defaultdict(list)
        for row in included:
            groups[tuple(row[field] for field in KEY) + (row["_dose"],)].append(row)
        base = {key[:-1]: values for key, values in groups.items() if key[-1] == 0}
        outputs = []
        for group_key, values in sorted(groups.items(), key=lambda item: tuple(map(str, item[0]))):
            dose = group_key[-1]
            if dose == 0:
                continue
            context = group_key[:-1]
            baseline = base.get(context)
            if baseline is None:
                raise ValueError("positive-dose context lacks zero-dose DA results: {}".format(context))
            target_label = context[KEY.index("target_label")]
            profiler = context[KEY.index("profiler")]
            target = alias_map.get((target_label, profiler))
            if target is None:
                raise ValueError("missing target alias for {} {}".format(target_label, profiler))
            target_rows = [row for row in values if row["feature"] == target]
            baseline_targets = [row for row in baseline if row["feature"] == target]
            if len(target_rows) != 1 or len(baseline_targets) != 1:
                raise ValueError("target feature must occur exactly once at baseline and dose")
            for threshold in thresholds:
                calls = {row["feature"] for row in values if row["_effect"] > 0 and row["_q"] <= threshold}
                baseline_calls = {row["feature"] for row in baseline if row["_effect"] > 0 and row["_q"] <= threshold}
                target_called = int(target in calls)
                precision = target_called / len(calls) if calls else 0.0
                recall = float(target_called)
                f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
                union = calls | baseline_calls
                jaccard = len(calls & baseline_calls) / len(union) if union else 1.0
                output = {field: value for field, value in zip(KEY, context)}
                output.update({
                    "spike_fraction_target": render(dose), "q_threshold": render(threshold),
                    "target_alias": target, "target_called": str(target_called),
                    "enriched_calls": str(len(calls)),
                    "off_target_enriched_calls": str(len(calls - {target})),
                    "precision": render(precision), "recall": render(recall), "f1": render(f1),
                    "target_effect": render(target_rows[0]["_effect"]),
                    "target_q_value": render(target_rows[0]["_q"]),
                    "target_effect_change_from_baseline": render(target_rows[0]["_effect"] - baseline_targets[0]["_effect"]),
                    "biomarker_set_jaccard_vs_baseline": render(jaccard),
                })
                outputs.append(output)
        if not outputs:
            raise ValueError("no positive-dose contexts were evaluated")
        args.outdir.mkdir(parents=True, exist_ok=True)
        output_path = args.outdir / "biomarker_propagation_metrics.tsv"
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(outputs)
        exclusion_path = args.outdir / "excluded_da_rows.tsv"
        with exclusion_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=REQUIRED, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows([{field: row[field] for field in REQUIRED} for row in excluded])
        summary = args.outdir / "biomarker_propagation_summary.tsv"
        summary.write_text("metric\tvalue\ninput_rows\t{}\nincluded_rows\t{}\nexcluded_rows\t{}\nevaluated_rows\t{}\nstatus\tPASS\n".format(
            len(rows), len(included), len(excluded), len(outputs)), encoding="utf-8")
        checksum = args.outdir / "biomarker_propagation.sha256"
        checksum.write_text("".join("{}  {}\n".format(sha256(path), path.resolve()) for path in
                                    (args.calls, args.aliases, args.spike_panel, output_path,
                                     exclusion_path, summary)), encoding="utf-8")
        (args.outdir / "SUCCESS").write_text("evaluated_rows\t{}\nstatus\tPASS\n".format(len(outputs)), encoding="utf-8")
        print("[PASS] Evaluated {} biomarker-propagation rows".format(len(outputs)))
    except (ValueError, FileNotFoundError, KeyError) as error:
        raise SystemExit("[ERROR] {}".format(error)) from error


if __name__ == "__main__":
    main()
