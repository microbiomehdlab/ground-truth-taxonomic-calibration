#!/usr/bin/env python3
"""Audit native Bracken and MetaPhlAn profile totals without renormalizing them."""

from __future__ import annotations

import argparse
import csv
import hashlib
from decimal import Decimal, InvalidOperation
from pathlib import Path


OUTPUT_FIELDS = [
    "profiler", "profile", "sample_id", "rows", "species_rows",
    "reported_total", "species_total", "unclassified_total",
    "non_species_total", "excess_over_composition", "rounding_tolerance",
    "native_unit", "status",
]
PERCENT_TOLERANCE = 0.001  # Printed MetaPhlAn rounding, in percentage points.


def number(value: str, field: str, path: Path) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid {field}: {value!r}") from error


def sample_id(path: Path, suffix: str) -> str:
    name = path.name
    return name[:-len(suffix)] if suffix and name.endswith(suffix) else name


def audit_bracken(path: Path) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"taxonomy_lvl", "new_est_reads", "fraction_total_reads"}
    fields = set(rows[0].keys()) if rows else set()
    missing = required - fields
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    species = [row for row in rows if row["taxonomy_lvl"].strip().upper() == "S"]
    fraction_text = [row["fraction_total_reads"].strip() for row in species]
    fractions = [number(value, "fraction_total_reads", path)
                 for value in fraction_text]
    estimates = [number(row["new_est_reads"], "new_est_reads", path)
                 for row in species]
    if any(value < 0 for value in fractions + estimates):
        raise ValueError(f"{path}: negative native abundance")
    fraction_sum = sum(fractions)
    try:
        # Each displayed value may differ from its unrounded value by at most
        # half a unit in its final printed decimal place. Summing those bounds
        # gives a conservative, profile-specific rounding tolerance.
        rounding_tolerance = sum(
            0.5 * 10.0 ** min(0, Decimal(value).as_tuple().exponent)
            for value in fraction_text
        )
    except InvalidOperation as error:
        raise ValueError(f"{path}: invalid printed Bracken fraction") from error
    excess = max(0.0, fraction_sum - 1.0)
    if excess <= 1e-12:
        status = "PASS"
    elif excess <= rounding_tolerance + 1e-12:
        status = "PASS_ROUNDING_TOLERANCE"
    else:
        status = "FAIL_FRACTION_GT_ONE"
    return {
        "profiler": "kraken2_bracken", "profile": str(path.resolve()),
        "sample_id": sample_id(path, ".bracken.S.tsv"), "rows": len(rows),
        "species_rows": len(species), "reported_total": fraction_sum,
        "species_total": fraction_sum, "unclassified_total": "",
        "non_species_total": "", "excess_over_composition": excess,
        "rounding_tolerance": rounding_tolerance,
        "native_unit": "fraction_total_reads",
        "status": status,
    }


def metaphlan_rows(path: Path) -> list[list[str]]:
    result = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 3:
                raise ValueError(f"{path}: expected at least three columns")
            result.append(row)
    return result


def audit_metaphlan(path: Path) -> dict[str, object]:
    rows = metaphlan_rows(path)
    abundances = [number(row[2], "relative_abundance", path) for row in rows]
    if any(value < 0 for value in abundances):
        raise ValueError(f"{path}: negative native abundance")
    species = [value for row, value in zip(rows, abundances)
               if row[0].split("|")[-1].startswith("s__")]
    unclassified = sum(value for row, value in zip(rows, abundances)
                       if row[0].upper() == "UNCLASSIFIED")
    # MetaPhlAn can emit the same community mass at several taxonomic ranks.
    # Therefore the sum of every row is not a meaningful compositional total.
    # Identify terminal leaves and separately validate every reported rank.
    clades = {row[0] for row in rows if row[0].upper() != "UNCLASSIFIED"}
    nonleaves = set()
    for clade in clades:
        parts = clade.split("|")
        nonleaves.update("|".join(parts[:index]) for index in range(1, len(parts)))
    leaf_values = [
        value for row, value in zip(rows, abundances)
        if row[0].upper() == "UNCLASSIFIED" or row[0] not in nonleaves
    ]
    non_species_leaf_values = [
        value for row, value in zip(rows, abundances)
        if row[0].upper() != "UNCLASSIFIED"
        and row[0] not in nonleaves
        and not row[0].split("|")[-1].startswith("s__")
        and not row[0].split("|")[-1].startswith("t__")
    ]
    rank_totals = {}
    for row, value in zip(rows, abundances):
        clade = row[0]
        if clade.upper() == "UNCLASSIFIED":
            continue
        rank = clade.split("|")[-1].split("__", 1)[0]
        rank_totals[rank] = rank_totals.get(rank, 0.0) + value
    total = sum(leaf_values)
    species_total = sum(species)
    non_species = sum(non_species_leaf_values)
    upper_bound = 100.0 + PERCENT_TOLERANCE
    bad_ranks = [rank for rank, value in rank_totals.items() if value > upper_bound]
    if total > upper_bound:
        status = "FAIL_LEAF_PERCENT_GT_100"
    elif unclassified > upper_bound:
        status = "FAIL_UNCLASSIFIED_GT_100"
    elif bad_ranks:
        status = "FAIL_RANK_PERCENT_GT_100:" + ",".join(sorted(bad_ranks))
    else:
        status = "PASS"
    return {
        "profiler": "metaphlan4", "profile": str(path.resolve()),
        "sample_id": sample_id(path, ".metaphlan.tsv"), "rows": len(rows),
        "species_rows": len(species), "reported_total": total,
        "species_total": species_total, "unclassified_total": unclassified,
        "non_species_total": non_species,
        "excess_over_composition": max(0.0, total - 100.0),
        "rounding_tolerance": PERCENT_TOLERANCE,
        "native_unit": "relative_abundance_pct", "status": status,
    }


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bracken", type=Path, action="append", default=[])
    parser.add_argument("--metaphlan", type=Path, action="append", default=[])
    parser.add_argument(
        "--canonical", type=Path,
        help="Read unique included native profile paths from a canonical input TSV.",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def profiles_from_canonical(path: Path) -> tuple[list[Path], list[Path]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"profiler", "source_profile", "include"}
    fields = set(rows[0]) if rows else set()
    missing = required - fields
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    profiles: dict[str, set[Path]] = {
        "kraken2_bracken": set(), "metaphlan4": set()
    }
    for row in rows:
        if row["include"] != "1":
            continue
        profiler = row["profiler"]
        if profiler not in profiles:
            raise ValueError(f"{path}: unsupported profiler: {profiler!r}")
        profiles[profiler].add(Path(row["source_profile"]))
    return (sorted(profiles["kraken2_bracken"]),
            sorted(profiles["metaphlan4"]))


def main() -> None:
    args = parse_args()
    try:
        if args.canonical:
            canonical_bracken, canonical_metaphlan = profiles_from_canonical(args.canonical)
            args.bracken.extend(canonical_bracken)
            args.metaphlan.extend(canonical_metaphlan)
        args.bracken = sorted(set(args.bracken))
        args.metaphlan = sorted(set(args.metaphlan))
    except (OSError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
    if not args.bracken and not args.metaphlan:
        raise SystemExit("[ERROR] supply at least one native profile")
    for path in args.bracken + args.metaphlan:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"[ERROR] missing or empty profile: {path}")

    records = ([audit_bracken(path) for path in args.bracken] +
               [audit_metaphlan(path) for path in args.metaphlan])
    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / "profile_semantics_audit.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    failures = [record for record in records
                if not str(record["status"]).startswith("PASS")]
    if failures:
        for record in failures[:20]:
            print("[ERROR] {status}: {profile}".format(**record))
        if len(failures) > 20:
            print("[ERROR] ... and {} additional failures".format(len(failures) - 20))
        print("[INFO] Complete failure report: {}".format(output))
        raise SystemExit("[ERROR] {} native profiles failed validation".format(len(failures)))

    checksums = args.outdir / "profile_semantics_inputs.sha256"
    with checksums.open("w", encoding="utf-8") as handle:
        if args.canonical:
            handle.write(f"{digest(args.canonical)}  {args.canonical.resolve()}\n")
        for path in args.bracken + args.metaphlan:
            handle.write(f"{digest(path)}  {path.resolve()}\n")
        handle.write(f"{digest(output)}  {output.resolve()}\n")
    (args.outdir / "SUCCESS").write_text(
        f"profiles\t{len(records)}\nstatus\tPASS\n", encoding="utf-8"
    )
    print(f"[PASS] Audited {len(records)} native profiles")
    print(f"[INFO] {output}")


if __name__ == "__main__":
    main()
