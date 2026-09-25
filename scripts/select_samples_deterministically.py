#!/usr/bin/env python3
"""Select a balanced subset using stable accession hashes, not row order."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib

SELECTION_FIELDS = ("selection_rank", "selection_hash", "selection_seed")


def prefixed(prefix: str, field: str) -> str:
    return "%s_%s" % (prefix, field) if prefix else field


def resolve_selection_names(fields, existing_prefix, new_prefix):
    """Decide the output names for the selection-provenance columns.

    Without a collision nothing changes, so existing callers keep their exact
    historical output. With a collision the run fails unless the caller states
    explicit provenance prefixes: appending a second `selection_rank` beside the
    first is what produced the ambiguous historical Yachida header, and it is
    never done silently again.
    """
    colliding = [field for field in SELECTION_FIELDS if field in fields]
    if colliding and not (existing_prefix or new_prefix):
        raise SystemExit(
            "[ERROR] the input manifest already carries "
            + ", ".join(colliding)
            + "; rerunning would append a duplicate column of the same name. "
            "Pass --existing-selection-prefix (for example pilot) and "
            "--new-selection-prefix (for example independent) to keep both "
            "sets of provenance under unambiguous names.")
    rename = {}
    if existing_prefix:
        for field in colliding:
            rename[field] = prefixed(existing_prefix, field)
    output_fields = [rename.get(field, field) for field in fields]
    new_names = [prefixed(new_prefix, field) for field in SELECTION_FIELDS]
    duplicates = sorted({
        name for name in output_fields + new_names
        if (output_fields + new_names).count(name) > 1})
    if duplicates:
        raise SystemExit(
            "[ERROR] the requested output header is ambiguous; duplicate "
            "column(s): " + ", ".join(duplicates))
    return output_fields, rename, dict(zip(SELECTION_FIELDS, new_names))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--per-condition", type=int, required=True)
    parser.add_argument("--selection-seed", default="ground-truth-taxonomic-calibration-yachida-v1")
    parser.add_argument("--id-column", default="sample_id")
    parser.add_argument("--condition-column", default="Target_Condition")
    parser.add_argument(
        "--existing-selection-prefix", default="",
        help="Rename selection_rank/hash/seed already present in the input to "
             "<PREFIX>_selection_* (for example pilot), preserving them.")
    parser.add_argument(
        "--new-selection-prefix", default="",
        help="Write this run's provenance as <PREFIX>_selection_* (for example "
             "independent) instead of the bare selection_* names.")
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    duplicated = sorted({name for name in fields if fields.count(name) > 1})
    if duplicated:
        raise SystemExit(
            "[ERROR] input manifest has duplicate column(s): "
            + ", ".join(duplicated))
    for field in (args.id_column, args.condition_column):
        if field not in fields:
            raise SystemExit(f"[ERROR] missing column {field!r}")
    output_fields, rename, new_names = resolve_selection_names(
        fields, args.existing_selection_prefix, args.new_selection_prefix)

    seen = set()
    grouped: dict[str, list[tuple[str, dict[str, str]]]] = {}
    for row in rows:
        sample = row[args.id_column]
        if not sample or sample in seen:
            raise SystemExit(f"[ERROR] blank or duplicate sample ID: {sample!r}")
        seen.add(sample)
        key = hashlib.sha256(
            f"stable-selection-v1\0{args.selection_seed}\0{sample}".encode("utf-8")
        ).hexdigest()
        grouped.setdefault(row[args.condition_column], []).append((key, row))

    selected = []
    for condition in sorted(grouped):
        candidates = sorted(grouped[condition], key=lambda item: (item[0], item[1][args.id_column]))
        if len(candidates) < args.per_condition:
            raise SystemExit(f"[ERROR] {condition}: need {args.per_condition}, found {len(candidates)}")
        for rank, (key, row) in enumerate(candidates[: args.per_condition], 1):
            record = {rename.get(name, name): value for name, value in row.items()}
            record[new_names["selection_rank"]] = str(rank)
            record[new_names["selection_hash"]] = key
            record[new_names["selection_seed"]] = args.selection_seed
            selected.append(record)

    output_fields = output_fields + [new_names[field] for field in SELECTION_FIELDS]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        f"{digest}  {args.output.name}\n", encoding="utf-8"
    )
    print(f"[OK] selected {len(selected)} samples across {len(grouped)} conditions -> {args.output}")


if __name__ == "__main__":
    main()
