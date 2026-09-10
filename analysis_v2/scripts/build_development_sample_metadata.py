#!/usr/bin/env python3
"""Build canonical-matched Feng/Zeller metadata for development analyses."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical", required=True, type=Path)
    parser.add_argument("--feng-manifest", required=True, type=Path)
    parser.add_argument("--zeller-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    canonical = read_tsv(args.canonical)
    if not canonical or not {"cohort", "sample_id"}.issubset(canonical[0]):
        raise SystemExit("[ERROR] canonical input lacks cohort/sample_id")

    wanted = {(row["cohort"], row["sample_id"]) for row in canonical if row.get("include") == "1"}
    if not wanted:
        raise SystemExit("[ERROR] canonical input contains no included samples")

    records = {}
    for cohort, path in (("feng", args.feng_manifest), ("zeller", args.zeller_manifest)):
        rows = read_tsv(path)
        required = {"sample_id", "condition", "study", "age", "sex", "bmi"}
        if not rows or not required.issubset(rows[0]):
            raise SystemExit(f"[ERROR] incomplete metadata manifest: {path}")
        for row in rows:
            key = (cohort, row["sample_id"])
            if key in records:
                raise SystemExit(f"[ERROR] duplicate metadata key: {key}")
            records[key] = {field: row[field] for field in required}
            records[key]["cohort"] = cohort

    missing = sorted(wanted - set(records))
    if missing:
        raise SystemExit(
            "[ERROR] metadata missing canonical samples: " + ", ".join(map(str, missing))
        )

    selected = [records[key] for key in sorted(wanted)]
    sample_ids = [row["sample_id"] for row in selected]
    if len(sample_ids) != len(set(sample_ids)):
        raise SystemExit("[ERROR] sample_id is not globally unique across cohorts")
    if any(not row["age"] or not row["sex"] for row in selected):
        raise SystemExit("[ERROR] primary age/sex metadata is incomplete")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "condition", "study", "age", "sex", "bmi", "cohort"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(selected)

    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        f"{digest(args.output)}  {args.output.name}\n", encoding="utf-8"
    )
    print(f"[PASS] Development metadata: {len(selected)} canonical-matched samples")


if __name__ == "__main__":
    main()
