#!/usr/bin/env python3
"""Compare implanted-target abundances between two Yachida result trees."""

import argparse
import csv
import math
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--aliases", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tolerance-pp", type=float, default=0.001)
    return parser.parse_args()


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_targets(panel_path, alias_path):
    targets = {}
    for row in read_tsv(panel_path):
        targets[row["label"]] = {"canonical": row["taxon_name"], "aliases": {}}
    with alias_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for target in targets.values():
                if row["canonical"] == target["canonical"]:
                    target["aliases"].setdefault(row["tool"], set()).add(row["alias"])
    for target in targets.values():
        for tool in ("kraken2_bracken", "metaphlan4"):
            target["aliases"].setdefault(tool, set()).add(target["canonical"])
    return targets


def normalize_name(value):
    value = value.strip().replace("_", " ")
    if "|" in value:
        value = value.rsplit("|", 1)[-1]
    if value.startswith("s  "):
        value = value[3:]
    elif value.startswith("s__"):
        value = value[3:]
    return " ".join(value.split()).casefold()


def bracken_abundances(path):
    values = {}
    for row in read_tsv(path):
        values[normalize_name(row["name"])] = float(row["fraction_total_reads"]) * 100.0
    return values


def metaphlan_abundances(path):
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and "s__" in fields[0]:
                values[normalize_name(fields[0])] = float(fields[2])
    return values


def lookup(values, names):
    matches = [(name, values[normalize_name(name)]) for name in names if normalize_name(name) in values]
    # Profiler tables conventionally omit zero-abundance features. Absence is
    # therefore a valid zero, whereas more than one matching alias is ambiguous.
    if not matches:
        return "", 0.0, 0
    if len(matches) > 1:
        return "", math.nan, len(matches)
    return matches[0][0], matches[0][1], 1


def profile_metadata(relative):
    parts = relative.parts
    if "independent" in parts:
        index = parts.index("independent")
        return "independent", [parts[index + 1]], parts[index + 2]
    if "community" in parts:
        index = parts.index("community")
        return "community", None, parts[index + 1]
    return None, None, None


def fraction_from_id(profile_id):
    tag = profile_id.rsplit("_f", 1)[-1]
    return float("0." + tag[2:]) if tag.startswith("0p") else float(tag.replace("p", "."))


def main():
    args = arguments()
    targets = load_targets(args.panel, args.aliases)
    rows = []
    reference_profiles = args.reference_root / "profiles"
    candidate_profiles = args.candidate_root / "profiles"
    patterns = (("kraken2_bracken", "*.bracken.S.tsv", bracken_abundances),
                ("metaphlan4", "*.metaphlan.tsv", metaphlan_abundances))
    for tool, pattern, reader in patterns:
        for reference_file in sorted(reference_profiles.rglob(pattern)):
            relative = reference_file.relative_to(reference_profiles)
            design, labels, profile_id = profile_metadata(relative)
            if design is None:
                continue
            candidate_file = candidate_profiles / relative
            if not candidate_file.is_file():
                raise SystemExit(f"[ERROR] Missing candidate profile: {candidate_file}")
            reference_values, candidate_values = reader(reference_file), reader(candidate_file)
            for label in labels or list(targets):
                target = targets[label]
                names = target["aliases"][tool]
                reference_name, reference_value, reference_matches = lookup(reference_values, names)
                candidate_name, candidate_value, candidate_matches = lookup(candidate_values, names)
                unambiguous = reference_matches <= 1 and candidate_matches <= 1
                difference = abs(candidate_value - reference_value) if unambiguous else math.nan
                relative_difference = (difference / abs(reference_value) * 100.0
                                       if unambiguous and reference_value != 0 else math.nan)
                rows.append({
                    "design": design, "profile_id": profile_id,
                    "spike_fraction": f"{fraction_from_id(profile_id):.10g}",
                    "target_label": label, "canonical_target": target["canonical"], "profiler": tool,
                    "reference_feature": reference_name, "candidate_feature": candidate_name,
                    "reference_abundance_pct": reference_value, "candidate_abundance_pct": candidate_value,
                    "absolute_difference_pp": difference, "relative_difference_pct": relative_difference,
                    "reference_matches": reference_matches, "candidate_matches": candidate_matches,
                    "within_tolerance": unambiguous and difference <= args.tolerance_pp,
                })
    if not rows:
        raise SystemExit("[ERROR] No spiked profiles were found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    failures = [row for row in rows if not row["within_tolerance"]]
    maximum = max((row["absolute_difference_pp"] for row in rows if math.isfinite(row["absolute_difference_pp"])), default=math.nan)
    print(f"[INFO] Compared {len(rows)} implanted-target abundance values")
    print(f"[INFO] Maximum absolute difference: {maximum:.12g} percentage points")
    print(f"[INFO] Tolerance: {args.tolerance_pp:.12g} percentage points")
    print(f"[INFO] Output: {args.output}")
    if failures:
        raise SystemExit(f"[FAIL] {len(failures)} target values failed matching or tolerance")
    print("[PASS] Every implanted-target abundance is within tolerance")


if __name__ == "__main__":
    main()
