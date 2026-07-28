#!/usr/bin/env python3
"""Resolve before/after FastQC archives for metadata-selected MetaPrep samples."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metaprep-root", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-list", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--sample-column", default=None)
    return parser.parse_args()


def read_metadata(path: Path, requested_column: str | None) -> tuple[str, list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        dialect = csv.Sniffer().sniff(handle.read(8192), delimiters="\t,")
        handle.seek(0)
        reader = csv.DictReader(handle, dialect=dialect)
        fields = reader.fieldnames or []
        candidates = [
            requested_column,
            "sample_id",
            "Name",
            "Sample Source ID",
            "sample",
            "SampleID",
        ]
        sample_column = next((x for x in candidates if x and x in fields), None)
        if sample_column is None:
            raise ValueError(
                f"Could not identify the sample column in {path}. Columns: {fields}"
            )
        samples = []
        seen = set()
        for row in reader:
            sample = (row.get(sample_column) or "").strip()
            if sample and sample not in seen:
                seen.add(sample)
                samples.append(sample)
    if not samples:
        raise ValueError(f"No sample identifiers found in {path}:{sample_column}")
    return sample_column, samples


def index_sample_directories(
    root: Path,
) -> tuple[dict[str, list[Path]], dict[tuple[str, str, int], list[Path]]]:
    candidates = list(root.glob("sequencing/*"))
    candidates.extend(root.glob("*/sequencing/*"))
    directory_index: dict[str, list[Path]] = {}
    accession_index: dict[tuple[str, str, int], list[Path]] = {}
    for path in candidates:
        if path.is_dir():
            resolved = path.resolve()
            directory_index.setdefault(path.name, []).append(resolved)
            for stage in ("qc_before", "qc_after"):
                for archive in sorted((path / stage).glob("*_fastqc.zip")):
                    match = re.fullmatch(r"(.+)_([12])_fastqc[.]zip", archive.name)
                    if match:
                        key = (match.group(1), stage, int(match.group(2)))
                        accession_index.setdefault(key, []).append(archive.resolve())
    return directory_index, accession_index


def resolve_one(stage_dir: Path, mate: int) -> Path:
    hits = sorted(stage_dir.glob(f"*_{mate}_fastqc.zip"))
    if len(hits) != 1:
        raise ValueError(
            f"Expected exactly one mate-{mate} FastQC ZIP in {stage_dir}; "
            f"found {len(hits)}"
        )
    if hits[0].stat().st_size == 0:
        raise ValueError(f"FastQC ZIP is empty: {hits[0]}")
    return hits[0].resolve()


def main() -> int:
    args = parse_args()
    root = args.metaprep_root.resolve()
    if not root.is_dir():
        raise ValueError(f"MetaPrep root does not exist: {root}")
    if not args.metadata.is_file():
        raise ValueError(f"Metadata file does not exist: {args.metadata}")

    sample_column, samples = read_metadata(args.metadata, args.sample_column)
    sample_dirs, accession_files = index_sample_directories(root)
    rows = []
    errors = []

    for sample in samples:
        locations = sample_dirs.get(sample, [])
        try:
            if len(locations) == 1:
                sample_dir = locations[0]
                resolved_files = {
                    (stage, mate): resolve_one(sample_dir / stage, mate)
                    for stage in ("qc_before", "qc_after")
                    for mate in (1, 2)
                }
                matched_by = "sample_directory"
            elif len(locations) > 1:
                raise ValueError(
                    f"expected at most one MetaPrep sample directory named "
                    f"{sample}; found {len(locations)}"
                )
            else:
                resolved_files = {}
                for stage in ("qc_before", "qc_after"):
                    for mate in (1, 2):
                        hits = accession_files.get((sample, stage, mate), [])
                        if len(hits) != 1:
                            raise ValueError(
                                f"no directory named {sample}, and expected exactly "
                                f"one archive {sample}_{mate}_fastqc.zip under "
                                f"{stage}; found {len(hits)}"
                            )
                        resolved_files[(stage, mate)] = hits[0]
                biological_dirs = {
                    path.parent.parent.resolve() for path in resolved_files.values()
                }
                if len(biological_dirs) != 1:
                    raise ValueError(
                        "accession-matched archives span multiple biological "
                        f"sample directories: {sorted(map(str, biological_dirs))}"
                    )
                matched_by = "fastqc_accession"

            for stage in ("qc_before", "qc_after"):
                for mate in (1, 2):
                    archive = resolved_files[(stage, mate)]
                    if archive.stat().st_size == 0:
                        raise ValueError(f"FastQC ZIP is empty: {archive}")
                    rows.append(
                        {
                            "sample_id": sample,
                            "matched_by": matched_by,
                            "metaprep_sample_directory": str(
                                archive.parent.parent.resolve()
                            ),
                            "stage": stage,
                            "mate": str(mate),
                            "fastqc_zip": str(archive),
                        }
                    )
        except ValueError as exc:
            errors.append(f"{sample}: {exc}")

    if errors:
        preview = "\n".join(f"  - {x}" for x in errors[:30])
        suffix = (
            f"\n  ... and {len(errors) - 30} more"
            if len(errors) > 30
            else ""
        )
        raise ValueError(
            f"MetaPrep QC validation failed for {len(errors)} of "
            f"{len(samples)} metadata samples:\n{preview}{suffix}"
        )

    args.output_list.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_list.write_text(
        "".join(f"{row['fastqc_zip']}\n" for row in rows), encoding="utf-8"
    )
    with args.output_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "matched_by",
                "metaprep_sample_directory",
                "stage",
                "mate",
                "fastqc_zip",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    expected = len(samples) * 4
    if len(rows) != expected:
        raise RuntimeError(f"Internal error: expected {expected} files; found {len(rows)}")
    print(
        f"[PASS] Resolved {len(rows)} MetaPrep FastQC ZIPs for "
        f"{len(samples)} samples "
        f"from metadata column '{sample_column}'."
    )
    print(f"[INFO] List:     {args.output_list}")
    print(f"[INFO] Manifest: {args.output_manifest}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
