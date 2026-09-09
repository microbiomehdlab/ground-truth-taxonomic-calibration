#!/usr/bin/env python3
"""Compare original and clean assembly arms for detection and biomarker outcomes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def render(value: float) -> str:
    return format(value, ".17g")


def exact_mcnemar(gains: int, losses: int) -> float:
    discordant = gains + losses
    if discordant == 0:
        return 1.0
    tail = sum(math.comb(discordant, k) for k in range(min(gains, losses) + 1)) / (2 ** discordant)
    return min(1.0, 2 * tail)


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values); running = 1.0; count = len(values)
    for rank_index in range(count - 1, -1, -1):
        index = order[rank_index]
        running = min(running, values[index] * count / (rank_index + 1))
        adjusted[index] = running
    return adjusted


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--biomarker-metrics", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    try:
        canonical = read(args.canonical)
        metrics = read(args.biomarker_metrics)
        needed = {"sample_id", "target_label", "assembly_arm", "profiler",
                  "spike_fraction_target", "detected_native_nonzero", "include"}
        if not canonical or not needed.issubset(canonical[0]):
            raise ValueError("canonical input lacks assembly-detection columns")
        canonical = [row for row in canonical if row["include"] == "1" and
                     row["analysis_population"] == "independent" and
                     float(row["spike_fraction_target"]) > 0]
        cell: dict[tuple[str, ...], dict[str, int]] = defaultdict(dict)
        for row in canonical:
            key = (row["sample_id"], row["target_label"], row["profiler"],
                   render(float(row["spike_fraction_target"])))
            arm = row["assembly_arm"]
            if arm in cell[key]:
                raise ValueError(f"duplicate detection arm: {key} {arm}")
            cell[key][arm] = int(row["detected_native_nonzero"])
        if not cell or any(set(arms) != {"original", "clean"} for arms in cell.values()):
            raise ValueError("every detection cell must contain original and clean arms")
        detection_pairs = []
        for key, arms in sorted(cell.items()):
            original, clean = arms["original"], arms["clean"]
            detection_pairs.append({
                "sample_id": key[0], "target_label": key[1], "profiler": key[2],
                "spike_fraction_target": key[3], "original_detected": original,
                "clean_detected": clean, "clean_minus_original": clean - original,
            })
        grouped: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
        for row in detection_pairs:
            grouped[(str(row["target_label"]), str(row["profiler"]),
                     str(row["spike_fraction_target"]))].append(row)
        detection_summary = []
        for key, values in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1], float(item[0][2]))):
            original = [int(row["original_detected"]) for row in values]
            clean = [int(row["clean_detected"]) for row in values]
            gains = sum(c > o for o, c in zip(original, clean)); losses = sum(c < o for o, c in zip(original, clean))
            detection_summary.append({
                "target_label": key[0], "profiler": key[1], "spike_fraction_target": key[2],
                "samples": len(values), "original_detection_rate": render(sum(original) / len(values)),
                "clean_detection_rate": render(sum(clean) / len(values)),
                "clean_minus_original_detection_rate": render((sum(clean) - sum(original)) / len(values)),
                "clean_only_detections": gains, "original_only_detections": losses,
                "exact_mcnemar_p_value": render(exact_mcnemar(gains, losses)),
            })
        adjusted = bh([float(row["exact_mcnemar_p_value"]) for row in detection_summary])
        for row, value in zip(detection_summary, adjusted):
            row["q_value_bh_across_target_profiler_dose"] = render(value)

        metric_needed = {"target_label", "assembly_arm", "profiler", "contrast",
                         "spike_fraction_target", "q_threshold", "target_called",
                         "target_effect", "off_target_enriched_calls", "precision"}
        if not metrics or not metric_needed.issubset(metrics[0]):
            raise ValueError("biomarker metrics lack assembly-comparison columns")
        metric_cell: dict[tuple[str, ...], dict[str, dict[str, str]]] = defaultdict(dict)
        for row in metrics:
            key = tuple(row[field] for field in
                        ("target_label", "profiler", "contrast", "spike_fraction_target", "q_threshold"))
            if row["assembly_arm"] in metric_cell[key]:
                raise ValueError(f"duplicate biomarker arm: {key}")
            metric_cell[key][row["assembly_arm"]] = row
        if not metric_cell or any(set(arms) != {"original", "clean"} for arms in metric_cell.values()):
            raise ValueError("every biomarker context must contain original and clean arms")
        biomarker = []
        for key, arms in sorted(metric_cell.items(), key=lambda item: tuple(map(str, item[0]))):
            original, clean = arms["original"], arms["clean"]
            record: dict[str, object] = dict(zip(
                ("target_label", "profiler", "contrast", "spike_fraction_target", "q_threshold"), key))
            for field in ("target_called", "target_effect", "off_target_enriched_calls", "precision"):
                old, new = float(original[field]), float(clean[field])
                record[f"original_{field}"] = render(old)
                record[f"clean_{field}"] = render(new)
                record[f"clean_minus_original_{field}"] = render(new - old)
            biomarker.append(record)

        args.outdir.mkdir(parents=True, exist_ok=True)
        detection_pair_fields = list(detection_pairs[0])
        detection_summary_fields = list(detection_summary[0])
        biomarker_fields = list(biomarker[0])
        paths = [args.outdir / "paired_detection_outcomes.tsv",
                 args.outdir / "detection_arm_comparison.tsv",
                 args.outdir / "biomarker_arm_comparison.tsv"]
        write(paths[0], detection_pairs, detection_pair_fields)
        write(paths[1], detection_summary, detection_summary_fields)
        write(paths[2], biomarker, biomarker_fields)
        summary = args.outdir / "comparison_summary.tsv"
        summary.write_text(
            "metric\tvalue\n"
            f"paired_detection_rows\t{len(detection_pairs)}\n"
            f"detection_summary_rows\t{len(detection_summary)}\n"
            f"biomarker_comparison_rows\t{len(biomarker)}\n"
            "inference\tdetection_exact_mcnemar;biomarker_descriptive\n"
            "status\tPASS\n", encoding="utf-8")
        seal = args.outdir / "assembly_outcome_comparison.sha256"
        seal.write_text("".join(f"{digest(path)}  {path.resolve()}\n" for path in
                                (args.canonical, args.biomarker_metrics, *paths, summary)), encoding="utf-8")
        (args.outdir / "SUCCESS").write_text("analysis\tassembly_outcome_comparison\nstatus\tPASS\n", encoding="utf-8")
        print(f"[PASS] Assembly outcome comparison: {len(detection_pairs)} paired detection cells")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
