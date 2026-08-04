#!/usr/bin/env python3
"""Join curated Yachida metadata to the complete PRJDB4176 FASTQ inventory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import urllib.request


ENA_REPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB4176"
    "&result=read_run&fields=run_accession,sample_accession,fastq_ftp,"
    "fastq_md5,fastq_bytes&format=tsv"
)
EXPECTED_COUNTS = {"Control": 290, "Adenoma": 67, "CRC": 258}


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(url: str, path: pathlib.Path) -> None:
    if path.is_file() and path.stat().st_size:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with urllib.request.urlopen(url) as response, partial.open("wb") as output:
        while block := response.read(1024 * 1024):
            output.write(block)
    partial.replace(path)


def paired(values: str, md5s: str, sizes: str) -> tuple[tuple[str, str, str], tuple[str, str, str]]:
    records = list(zip(values.split(";"), md5s.split(";"), sizes.split(";")))
    r1 = [record for record in records if record[0].endswith("_1.fastq.gz")]
    r2 = [record for record in records if record[0].endswith("_2.fastq.gz")]
    if len(r1) != 1 or len(r2) != 1:
        raise ValueError(f"expected exactly one paired FASTQ set, found R1={len(r1)} R2={len(r2)}")
    return r1[0], r2[0]


def read_unique(path: pathlib.Path, key: str) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    if key not in fields:
        raise ValueError(f"{path} lacks required column {key!r}")
    values = [row[key].strip() for row in rows]
    if any(not value for value in values) or len(values) != len(set(values)):
        raise ValueError(f"{path}: {key!r} contains blank or duplicate values")
    return rows, fields


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-metadata", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("cohort_manifest.tsv"))
    parser.add_argument("--cache-dir", type=pathlib.Path, default=pathlib.Path("metadata_cache"))
    args = parser.parse_args()

    ena_path = args.cache_dir / "PRJDB4176.ena.tsv"
    fetch(ENA_REPORT, ena_path)
    metadata, _ = read_unique(args.study_metadata, "Sample Source ID")
    if len(metadata) != 615:
        raise ValueError(f"expected 615 curated metadata samples, found {len(metadata)}")

    with ena_path.open(newline="", encoding="utf-8") as handle:
        ena_rows = list(csv.DictReader(handle, delimiter="\t"))
    by_sample: dict[str, list[dict[str, str]]] = {}
    for row in ena_rows:
        by_sample.setdefault(row["sample_accession"], []).append(row)

    rows = []
    for source in metadata:
        biosample = source["Sample Source ID"].strip()
        records = by_sample.get(biosample, [])
        if len(records) != 1:
            raise ValueError(f"{biosample}: expected exactly one PRJDB4176 run, found {len(records)}")
        record = records[0]
        r1, r2 = paired(record["fastq_ftp"], record["fastq_md5"], record["fastq_bytes"])
        condition = source["Study condition"].strip()
        if condition not in EXPECTED_COUNTS:
            raise ValueError(f"{biosample}: unrecognized Study condition {condition!r}")
        rows.append({
            "sample_id": biosample,
            "run_accession": record["run_accession"],
            "biosample_accession": biosample,
            "Study": "YachidaS_2019",
            "Target_Condition": condition,
            "age": source["Age"].strip(),
            "sex": source["Sex"].strip(),
            "bmi": source["BMI"].strip(),
            "tumor_stage_ajcc": source["Tumor Staging AJCC"].strip(),
            "primary_tumor_location": source["Primary Tumor Location"].strip(),
            "fastq1_url": "https://" + r1[0],
            "fastq2_url": "https://" + r2[0],
            "fastq1_md5": r1[1],
            "fastq2_md5": r2[1],
            "fastq1_bytes": r1[2],
            "fastq2_bytes": r2[2],
        })

    rows.sort(key=lambda row: row["sample_id"])
    counts = {condition: sum(row["Target_Condition"] == condition for row in rows)
              for condition in EXPECTED_COUNTS}
    if counts != EXPECTED_COUNTS:
        raise ValueError(f"unexpected curated condition counts: {counts}")

    included = {row["biosample_accession"] for row in rows}
    excluded = sorted(
        (row for row in ena_rows if row["sample_accession"] not in included),
        key=lambda row: (row["sample_accession"], row["run_accession"]),
    )
    if len(excluded) != 30:
        raise ValueError(f"expected 30 project runs outside curated metadata, found {len(excluded)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    excluded_path = args.output.with_name(args.output.stem + ".excluded_project_runs.tsv")
    with excluded_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(excluded[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(excluded)
    for path in (args.output, excluded_path):
        digest = sha256(path)
        path.with_suffix(path.suffix + ".sha256").write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    provenance_path = args.output.with_name(args.output.stem + ".provenance.tsv")
    provenance = [
        ("manifest_schema", "yachida-prjdb4176-manifest-v1"),
        ("study_metadata_file", args.study_metadata.name),
        ("study_metadata_sha256", sha256(args.study_metadata)),
        ("ena_inventory_url", ENA_REPORT),
        ("ena_inventory_sha256", sha256(ena_path)),
        ("curated_samples", str(len(rows))),
        ("excluded_project_runs", str(len(excluded))),
        *[(f"condition_{condition}", str(counts[condition])) for condition in EXPECTED_COUNTS],
    ]
    with provenance_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("field", "value")); writer.writerows(provenance)
    print(f"[OK] {args.output}: 615 curated samples ({counts})")
    print(f"[OK] Excluded project runs not in curated metadata: {len(excluded)} -> {excluded_path}")
    print(f"[OK] Provenance: {provenance_path}")


if __name__ == "__main__":
    main()
