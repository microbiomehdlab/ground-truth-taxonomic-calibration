#!/usr/bin/env python3
"""Build the frozen DRA006684 (Yachida 2019) paired-read manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import urllib.request
import xml.etree.ElementTree as ET


DDBJ_BASE = "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA006/DRA006684"
ENA_REPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=DRA006684"
    "&result=read_run&fields=run_accession,sample_accession,fastq_ftp,"
    "fastq_md5,fastq_bytes&format=tsv"
)


def fetch(url: str, path: pathlib.Path) -> None:
    if path.is_file() and path.stat().st_size:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with urllib.request.urlopen(url) as response, partial.open("wb") as output:
        while block := response.read(1024 * 1024):
            output.write(block)
    partial.replace(path)


def attributes(sample: ET.Element) -> dict[str, str]:
    result = {}
    for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        tag = item.findtext("TAG", "").strip()
        if tag:
            result[tag] = item.findtext("VALUE", "").strip()
    return result


def paired(values: str, md5s: str, sizes: str) -> tuple[tuple[str, str, str], tuple[str, str, str]]:
    records = list(zip(values.split(";"), md5s.split(";"), sizes.split(";")))
    r1 = [record for record in records if record[0].endswith("_1.fastq.gz")]
    r2 = [record for record in records if record[0].endswith("_2.fastq.gz")]
    if len(r1) != 1 or len(r2) != 1:
        raise ValueError(f"expected exactly one paired FASTQ set, found R1={len(r1)} R2={len(r2)}")
    return r1[0], r2[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("cohort_manifest.tsv"))
    parser.add_argument("--cache-dir", type=pathlib.Path, default=pathlib.Path("metadata_cache"))
    args = parser.parse_args()

    files = {
        "run": args.cache_dir / "DRA006684.run.xml",
        "experiment": args.cache_dir / "DRA006684.experiment.xml",
        "sample": args.cache_dir / "DRA006684.sample.xml",
        "ena": args.cache_dir / "DRA006684.ena.tsv",
    }
    for kind in ("run", "experiment", "sample"):
        fetch(f"{DDBJ_BASE}/DRA006684.{kind}.xml", files[kind])
    fetch(ENA_REPORT, files["ena"])

    samples = {}
    for node in ET.parse(files["sample"]).getroot():
        attrs = attributes(node)
        samples[node.attrib["accession"]] = {
            "biosample_accession": node.findtext("./IDENTIFIERS/PRIMARY_ID", "").strip(),
            "sample_name": attrs.get("sample_name", ""),
            "host_disease_stat": attrs.get("host_disease_stat", ""),
            "age": attrs.get("age", ""),
            "sex": attrs.get("sex", ""),
        }

    experiments = {}
    for node in ET.parse(files["experiment"]).getroot():
        descriptor = node.find("./DESIGN/SAMPLE_DESCRIPTOR")
        if descriptor is None:
            raise ValueError(f"experiment {node.attrib['accession']} lacks SAMPLE_DESCRIPTOR")
        experiments[node.attrib["accession"]] = descriptor.attrib["accession"]

    runs = {}
    for node in ET.parse(files["run"]).getroot():
        reference = node.find("EXPERIMENT_REF")
        if reference is None:
            raise ValueError(f"run {node.attrib['accession']} lacks EXPERIMENT_REF")
        runs[node.attrib["accession"]] = reference.attrib["accession"]

    with files["ena"].open(newline="", encoding="utf-8") as handle:
        ena = {row["run_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}

    condition_map = {"Healthy control": "Control", "CRC (Stage III/IV)": "CRC"}
    rows = []
    for run_accession in sorted(runs):
        experiment_accession = runs[run_accession]
        sample_accession = experiments[experiment_accession]
        sample = samples[sample_accession]
        source_condition = sample["host_disease_stat"]
        if source_condition not in condition_map:
            raise ValueError(f"unrecognized condition for {run_accession}: {source_condition!r}")
        record = ena.get(run_accession)
        if record is None:
            raise ValueError(f"ENA report lacks {run_accession}")
        r1, r2 = paired(record["fastq_ftp"], record["fastq_md5"], record["fastq_bytes"])
        rows.append({
            "sample_id": run_accession,
            "run_accession": run_accession,
            "experiment_accession": experiment_accession,
            "sample_accession": sample_accession,
            "biosample_accession": sample["biosample_accession"],
            "source_sample_name": sample["sample_name"],
            "Study": "Yachida_2019",
            "Target_Condition": condition_map[source_condition],
            "source_condition": source_condition,
            "age": sample["age"],
            "sex": sample["sex"],
            "fastq1_url": "https://" + r1[0],
            "fastq2_url": "https://" + r2[0],
            "fastq1_md5": r1[1],
            "fastq2_md5": r2[1],
            "fastq1_bytes": r1[2],
            "fastq2_bytes": r2[2],
        })

    if len(rows) != 80:
        raise ValueError(f"expected 80 DRA006684 runs, found {len(rows)}")
    counts = {condition: sum(row["Target_Condition"] == condition for row in rows) for condition in condition_map.values()}
    if counts != {"Control": 40, "CRC": 40}:
        raise ValueError(f"unexpected condition counts: {counts}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        f"{digest}  {args.output.name}\n", encoding="utf-8"
    )
    print(f"[OK] {args.output}: 80 samples ({counts})")
    print(f"[OK] SHA-256: {digest}")


if __name__ == "__main__":
    main()
