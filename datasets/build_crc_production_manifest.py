#!/usr/bin/env python3
"""Join eligible samples and legacy URL inventory to an official ENA run report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import re
import urllib.parse


RUN_RE = re.compile(r"^([SED]RR\d+)_([12])\.fastq\.gz$")


def read_table(path: pathlib.Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(path: pathlib.Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalized_url(url: str) -> str:
    url = url.strip()
    if url.startswith("ftp://") or url.startswith("http://") or url.startswith("https://"):
        parsed = urllib.parse.urlparse(url)
        return f"https://{parsed.netloc}{parsed.path}"
    return f"https://{url}"


def parse_inventory(path: pathlib.Path) -> dict[str, dict[str, dict[int, str]]]:
    inventory: dict[str, dict[str, dict[int, str]]] = {}
    sample = ""
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith("# Sample: "):
            sample = line.removeprefix("# Sample: ").strip()
            if not sample or sample in inventory:
                raise SystemExit(f"[ERROR] blank or duplicate sample block: {sample!r}")
            inventory[sample] = {}
        elif line.startswith("wget "):
            if not sample:
                raise SystemExit("[ERROR] wget URL appears before a sample block")
            url = line.split()[-1]
            filename = pathlib.PurePosixPath(urllib.parse.urlparse(url).path).name
            match = RUN_RE.match(filename)
            if not match:
                raise SystemExit(f"[ERROR] unsupported paired FASTQ filename: {filename}")
            run, mate_text = match.groups()
            mate = int(mate_text)
            if mate in inventory[sample].setdefault(run, {}):
                raise SystemExit(f"[ERROR] duplicate mate URL: {sample} {run} mate {mate}")
            inventory[sample][run][mate] = normalized_url(url)
    return inventory


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eligible-manifest", required=True, type=pathlib.Path)
    parser.add_argument("--url-inventory", required=True, type=pathlib.Path)
    parser.add_argument("--ena-report", required=True, type=pathlib.Path)
    parser.add_argument("--independent-selection", type=pathlib.Path)
    parser.add_argument("--ena-study-accession", required=True)
    parser.add_argument("--ena-report-retrieved-date", required=True,
                        help="UTC retrieval date in YYYY-MM-DD form")
    parser.add_argument("--output", required=True, type=pathlib.Path,
                        help="One-row-per-biological-sample production manifest")
    args = parser.parse_args()

    eligible, eligible_fields = read_table(args.eligible_manifest)
    for field in ("Name", "Study condition", "Study name"):
        if field not in eligible_fields:
            raise SystemExit(f"[ERROR] eligible manifest lacks {field!r}")
    eligible_by_id = {row["Name"]: row for row in eligible}
    if len(eligible_by_id) != len(eligible):
        raise SystemExit("[ERROR] eligible sample IDs are not unique")

    selected_by_id: dict[str, dict[str, str]] = {}
    if args.independent_selection:
        selected, selection_fields = read_table(args.independent_selection)
        required_selection = {
            "Name", "Study condition", "selection_rank", "selection_hash", "selection_seed",
        }
        if missing := required_selection.difference(selection_fields):
            raise SystemExit(f"[ERROR] independent selection lacks fields: {sorted(missing)}")
        selected_by_id = {row["Name"]: row for row in selected}
        if len(selected_by_id) != len(selected):
            raise SystemExit("[ERROR] independent selection IDs are not unique")
        if not set(selected_by_id) <= set(eligible_by_id):
            raise SystemExit("[ERROR] independent selection is not nested in eligibility")

    inventory = parse_inventory(args.url_inventory)
    if set(inventory) != set(eligible_by_id):
        missing = sorted(set(eligible_by_id).difference(inventory))
        extra = sorted(set(inventory).difference(eligible_by_id))
        raise SystemExit(f"[ERROR] inventory/eligibility mismatch; missing={missing}, extra={extra}")

    ena_rows, ena_fields = read_table(args.ena_report)
    required_ena = {
        "run_accession", "sample_accession", "secondary_sample_accession",
        "fastq_ftp", "fastq_md5", "fastq_bytes",
    }
    if missing := required_ena.difference(ena_fields):
        raise SystemExit(f"[ERROR] ENA report lacks fields: {sorted(missing)}")
    ena_by_run = {row["run_accession"]: row for row in ena_rows}
    if len(ena_by_run) != len(ena_rows):
        raise SystemExit("[ERROR] ENA report contains duplicate run accessions")

    run_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for sample in sorted(inventory):
        metadata = eligible_by_id[sample]
        assembled = []
        for run in sorted(inventory[sample]):
            mates = inventory[sample][run]
            if set(mates) != {1, 2}:
                raise SystemExit(f"[ERROR] incomplete URL pair: {sample} {run}: {sorted(mates)}")
            if run not in ena_by_run:
                raise SystemExit(f"[ERROR] run absent from official ENA report: {run}")
            official = ena_by_run[run]
            official_urls = [normalized_url(value) for value in official["fastq_ftp"].split(";") if value]
            official_md5 = [value for value in official["fastq_md5"].split(";") if value]
            official_bytes = [value for value in official["fastq_bytes"].split(";") if value]
            if not (len(official_urls) == len(official_md5) == len(official_bytes) == 2):
                raise SystemExit(f"[ERROR] ENA run is not one paired FASTQ set: {run}")
            official_by_name = {
                pathlib.PurePosixPath(urllib.parse.urlparse(url).path).name: (url, md5, size)
                for url, md5, size in zip(official_urls, official_md5, official_bytes)
            }
            values = {}
            for mate in (1, 2):
                expected_name = f"{run}_{mate}.fastq.gz"
                if expected_name not in official_by_name:
                    raise SystemExit(f"[ERROR] ENA mate filename missing: {expected_name}")
                url, md5, size_text = official_by_name[expected_name]
                if pathlib.PurePosixPath(urllib.parse.urlparse(mates[mate]).path).name != expected_name:
                    raise SystemExit(f"[ERROR] legacy inventory disagrees with ENA: {sample} {run}")
                if not re.fullmatch(r"[0-9a-f]{32}", md5):
                    raise SystemExit(f"[ERROR] invalid ENA MD5: {sample} {run} mate {mate}")
                size = int(size_text)
                if size <= 0:
                    raise SystemExit(f"[ERROR] invalid ENA byte count: {sample} {run} mate {mate}")
                values[mate] = (url, md5, size)
            run_row = {
                "sample_id": sample,
                "condition": metadata["Study condition"],
                "study": metadata["Study name"],
                "run_accession": run,
                "ena_sample_accession": official["sample_accession"],
                "ena_secondary_sample_accession": official["secondary_sample_accession"],
                "fastq1_url": values[1][0], "fastq2_url": values[2][0],
                "fastq1_md5": values[1][1], "fastq2_md5": values[2][1],
                "fastq1_bytes": values[1][2], "fastq2_bytes": values[2][2],
            }
            run_rows.append(run_row)
            assembled.append(run_row)

        sample_rows.append({
            "sample_id": sample,
            "condition": metadata["Study condition"],
            "study": metadata["Study name"],
            "age": metadata.get("Age", ""),
            "sex": metadata.get("Sex", ""),
            "bmi": metadata.get("BMI", ""),
            "independent_subset": "1" if sample in selected_by_id else "0",
            "independent_selection_rank": selected_by_id.get(sample, {}).get("selection_rank", ""),
            "independent_selection_hash": selected_by_id.get(sample, {}).get("selection_hash", ""),
            "independent_selection_seed": selected_by_id.get(sample, {}).get("selection_seed", ""),
            "run_count": len(assembled),
            "run_accessions": ";".join(str(row["run_accession"]) for row in assembled),
            "fastq1_urls": ";".join(str(row["fastq1_url"]) for row in assembled),
            "fastq2_urls": ";".join(str(row["fastq2_url"]) for row in assembled),
            "fastq1_md5s": ";".join(str(row["fastq1_md5"]) for row in assembled),
            "fastq2_md5s": ";".join(str(row["fastq2_md5"]) for row in assembled),
            "fastq1_bytes": ";".join(str(row["fastq1_bytes"]) for row in assembled),
            "fastq2_bytes": ";".join(str(row["fastq2_bytes"]) for row in assembled),
            "total_download_bytes": sum(int(row["fastq1_bytes"]) + int(row["fastq2_bytes"]) for row in assembled),
            "run_combination_rule": "concatenate mates separately in ascending run_accession order",
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    run_path = args.output.with_name(f"{args.output.stem}.runs.tsv")
    independent_path = args.output.with_name(f"{args.output.stem}.independent.tsv")
    provenance_path = args.output.with_name(f"{args.output.stem}.provenance.tsv")
    sample_fields = list(sample_rows[0])
    run_fields = list(run_rows[0])
    write_table(args.output, sample_fields, sample_rows)
    write_table(run_path, run_fields, run_rows)
    write_table(
        independent_path,
        sample_fields,
        [row for row in sample_rows if row["independent_subset"] == "1"],
    )
    provenance = [
        {"key": "eligible_manifest", "value": str(args.eligible_manifest)},
        {"key": "eligible_manifest_sha256", "value": sha256(args.eligible_manifest)},
        {"key": "url_inventory", "value": str(args.url_inventory)},
        {"key": "url_inventory_sha256", "value": sha256(args.url_inventory)},
        {"key": "ena_report", "value": str(args.ena_report)},
        {"key": "ena_report_sha256", "value": sha256(args.ena_report)},
        {"key": "ena_study_accession", "value": args.ena_study_accession},
        {"key": "ena_report_retrieved_utc", "value": args.ena_report_retrieved_date},
        {"key": "independent_selection", "value": str(args.independent_selection or "")},
        {"key": "independent_selection_sha256",
         "value": sha256(args.independent_selection) if args.independent_selection else ""},
        {"key": "biological_samples", "value": str(len(sample_rows))},
        {"key": "paired_runs", "value": str(len(run_rows))},
        {"key": "independent_samples", "value": str(len(selected_by_id))},
        {"key": "run_order", "value": "ascending run_accession"},
        {"key": "multi_run_rule", "value": "concatenate R1 runs and R2 runs separately as gzip members"},
        {"key": "outcome_use", "value": "none"},
    ]
    write_table(provenance_path, ["key", "value"], provenance)
    for path in (args.output, run_path, independent_path, provenance_path):
        path.with_suffix(path.suffix + ".sha256").write_text(
            f"{sha256(path)}  {path.name}\n", encoding="utf-8"
        )
    print(f"[PASS] production manifest: {len(sample_rows)} samples, {len(run_rows)} paired runs")
    print(f"[INFO] Sample manifest: {args.output}")
    print(f"[INFO] Run manifest: {run_path}")
    print(f"[INFO] Independent manifest: {independent_path}")


if __name__ == "__main__":
    main()
