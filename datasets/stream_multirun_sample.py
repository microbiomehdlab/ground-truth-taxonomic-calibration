#!/usr/bin/env python3
"""Download, verify, assemble, and process one frozen multi-run sample."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import pathlib
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request


SENTINEL = ".ground_truth_stream_sample"


def checksum(path: pathlib.Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def positive_environment_integer(name: str, default: int) -> int:
    value = os.environ.get(name, str(default))
    try:
        parsed = int(value)
    except ValueError as error:
        raise SystemExit(f"[ERROR] {name} must be a positive integer; observed {value!r}") from error
    if parsed < 1:
        raise SystemExit(f"[ERROR] {name} must be a positive integer; observed {value!r}")
    return parsed


def download(url: str, destination: pathlib.Path, expected_md5: str, expected_bytes: int) -> None:
    if destination.is_file() and destination.stat().st_size == expected_bytes:
        if checksum(destination, "md5") == expected_md5:
            print(f"[SKIP] verified download: {destination}")
            return
    attempts = positive_environment_integer("CRC_DOWNLOAD_ATTEMPTS", 4)
    retry_seconds = positive_environment_integer("CRC_DOWNLOAD_RETRY_SECONDS", 15)
    timeout_seconds = positive_environment_integer("CRC_DOWNLOAD_TIMEOUT_SECONDS", 120)
    partial = destination.with_suffix(destination.suffix + ".partial")
    for attempt in range(1, attempts + 1):
        partial.unlink(missing_ok=True)
        try:
            with urllib.request.urlopen(url, timeout=timeout_seconds) as response, partial.open("wb") as output:
                shutil.copyfileobj(response, output, length=4 * 1024 * 1024)
            observed_bytes = partial.stat().st_size
            observed_md5 = checksum(partial, "md5")
            if observed_bytes != expected_bytes or observed_md5 != expected_md5:
                raise RuntimeError(
                    "size/MD5 verification failed "
                    f"(bytes={observed_bytes}/{expected_bytes}, md5={observed_md5}/{expected_md5})"
                )
            partial.replace(destination)
            return
        except (OSError, urllib.error.URLError, RuntimeError) as error:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise SystemExit(f"[ERROR] download failed after {attempts} attempts: {url}: {error}") from error
            delay = retry_seconds * attempt
            print(f"[WARN] Download attempt {attempt}/{attempts} failed: {url}: {error}; retrying in {delay}s",
                  file=sys.stderr)
            time.sleep(delay)


def split_field(row: dict[str, str], field: str) -> list[str]:
    values = [value for value in row[field].split(";") if value]
    if not values:
        raise SystemExit(f"[ERROR] empty manifest field: {field}")
    return values


def validate_manifest_lists(row: dict[str, str]) -> list[dict[str, object]]:
    fields = (
        "run_accessions", "fastq1_urls", "fastq2_urls", "fastq1_md5s", "fastq2_md5s",
        "fastq1_bytes", "fastq2_bytes",
    )
    values = {field: split_field(row, field) for field in fields}
    lengths = {len(items) for items in values.values()}
    if lengths != {int(row["run_count"])}:
        raise SystemExit(f"[ERROR] inconsistent multi-run manifest lengths: {sorted(lengths)}")
    result = []
    for index, run in enumerate(values["run_accessions"]):
        if index and run <= values["run_accessions"][index - 1]:
            raise SystemExit("[ERROR] run accessions are not strictly ascending")
        result.append({
            "run_order": index + 1,
            "run_accession": run,
            "fastq1_url": values["fastq1_urls"][index],
            "fastq2_url": values["fastq2_urls"][index],
            "fastq1_md5": values["fastq1_md5s"][index],
            "fastq2_md5": values["fastq2_md5s"][index],
            "fastq1_bytes": int(values["fastq1_bytes"][index]),
            "fastq2_bytes": int(values["fastq2_bytes"][index]),
        })
    return result


def concatenate_gzip_members(sources: list[pathlib.Path], destination: pathlib.Path) -> None:
    partial = destination.with_suffix(destination.suffix + ".partial")
    partial.unlink(missing_ok=True)
    with partial.open("wb") as output:
        for source in sources:
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, output, length=4 * 1024 * 1024)
    with gzip.open(partial, "rb") as handle:
        while handle.read(4 * 1024 * 1024):
            pass
    partial.replace(destination)


def verify_receipt(receipt: pathlib.Path, scratch: pathlib.Path) -> list[dict[str, str]]:
    if not receipt.is_file():
        raise SystemExit(f"[ERROR] runner did not create receipt: {receipt}")
    with receipt.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"path", "sha256", "bytes"}:
        raise SystemExit("[ERROR] receipt must contain path, sha256,bytes columns and at least one output")
    scratch = scratch.resolve()
    for row in rows:
        path = pathlib.Path(row["path"]).resolve()
        if path == scratch or scratch in path.parents:
            raise SystemExit(f"[ERROR] retained output is inside disposable scratch: {path}")
        if not path.is_file() or not path.stat().st_size:
            raise SystemExit(f"[ERROR] retained output missing/empty: {path}")
        if path.stat().st_size != int(row["bytes"]) or checksum(path, "sha256") != row["sha256"]:
            raise SystemExit(f"[ERROR] retained output checksum mismatch: {path}")
    return rows


def write_checksum_sidecar(path: pathlib.Path) -> pathlib.Path:
    sidecar = path.with_suffix(path.suffix + ".sha256")
    sidecar.write_text(f"{checksum(path, 'sha256')}  {path.name}\n", encoding="utf-8")
    return sidecar


def verify_checksum_sidecar(path: pathlib.Path) -> None:
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if not path.is_file() or not sidecar.is_file():
        raise SystemExit(f"[ERROR] input provenance or checksum missing: {path}")
    fields = sidecar.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1] != path.name or fields[0] != checksum(path, "sha256"):
        raise SystemExit(f"[ERROR] input provenance checksum mismatch: {path}")


def write_input_provenance(path: pathlib.Path, sample: str, runs: list[dict[str, object]],
                           assembled_r1: pathlib.Path, assembled_r2: pathlib.Path) -> None:
    fields = [
        "sample_id", "run_order", "run_accession", "fastq1_url", "fastq2_url",
        "fastq1_md5", "fastq2_md5", "fastq1_bytes", "fastq2_bytes",
        "assembled_r1", "assembled_r2", "assembled_r1_sha256", "assembled_r2_sha256",
        "assembled_r1_bytes", "assembled_r2_bytes", "assembly_rule",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        common = {
            "sample_id": sample,
            "assembled_r1": str(assembled_r1), "assembled_r2": str(assembled_r2),
            "assembled_r1_sha256": checksum(assembled_r1, "sha256"),
            "assembled_r2_sha256": checksum(assembled_r2, "sha256"),
            "assembled_r1_bytes": assembled_r1.stat().st_size,
            "assembled_r2_bytes": assembled_r2.stat().st_size,
            "assembly_rule": "concatenate mates separately in ascending run_accession order",
        }
        for run in runs:
            writer.writerow({**run, **common})
    write_checksum_sidecar(path)


def safe_cleanup(work: pathlib.Path) -> None:
    if not (work / SENTINEL).is_file():
        raise SystemExit(f"[ERROR] cleanup sentinel missing: {work / SENTINEL}")
    shutil.rmtree(work)
    print(f"[CLEANED] deleted verified disposable inputs: {work}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=pathlib.Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--scratch-root", required=True, type=pathlib.Path)
    parser.add_argument("--state-dir", required=True, type=pathlib.Path)
    parser.add_argument("--runner", required=True, type=pathlib.Path)
    parser.add_argument("--delete-inputs-after-verification", action="store_true")
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        matches = [row for row in csv.DictReader(handle, delimiter="\t") if row["sample_id"] == args.sample_id]
    if len(matches) != 1:
        raise SystemExit(f"[ERROR] expected one manifest row for {args.sample_id}, found {len(matches)}")
    row = matches[0]
    runs = validate_manifest_lists(row)
    if not args.runner.is_file() or not os.access(args.runner, os.X_OK):
        raise SystemExit(f"[ERROR] runner is missing or not executable: {args.runner}")

    root = args.scratch_root.resolve()
    if root == pathlib.Path("/"):
        raise SystemExit("[ERROR] scratch root cannot be /")
    work = (root / args.sample_id).resolve()
    if root not in work.parents:
        raise SystemExit("[ERROR] unsafe sample scratch path")
    state_dir = args.state_dir.resolve()
    if state_dir == root or root in state_dir.parents:
        raise SystemExit("[ERROR] persistent state directory cannot be inside disposable scratch")
    state_dir.mkdir(parents=True, exist_ok=True)
    receipt = state_dir / f"{args.sample_id}.retained_outputs.tsv"
    input_provenance = state_dir / f"{args.sample_id}.input_provenance.tsv"
    verified_marker = state_dir / f"{args.sample_id}.verified"
    if verified_marker.is_file() and receipt.is_file() and input_provenance.is_file():
        verify_checksum_sidecar(input_provenance)
        outputs = verify_receipt(receipt, work)
        print(f"[SKIP] already verified {len(outputs)} retained outputs for {args.sample_id}")
        if args.delete_inputs_after_verification and work.exists():
            safe_cleanup(work)
        return

    work.mkdir(parents=True, exist_ok=True)
    (work / SENTINEL).write_text("managed disposable sample directory\n", encoding="utf-8")
    run_root = work / "raw" / "runs"
    run_root.mkdir(parents=True, exist_ok=True)
    mate_paths: dict[int, list[pathlib.Path]] = {1: [], 2: []}
    for run in runs:
        accession = str(run["run_accession"])
        for mate in (1, 2):
            destination = run_root / f"{accession}_{mate}.fastq.gz"
            download(str(run[f"fastq{mate}_url"]), destination,
                     str(run[f"fastq{mate}_md5"]), int(run[f"fastq{mate}_bytes"]))
            mate_paths[mate].append(destination)

    assembled_root = work / "raw" / "assembled"
    assembled_root.mkdir(parents=True, exist_ok=True)
    assembled_r1 = assembled_root / f"{args.sample_id}_1.fastq.gz"
    assembled_r2 = assembled_root / f"{args.sample_id}_2.fastq.gz"
    concatenate_gzip_members(mate_paths[1], assembled_r1)
    concatenate_gzip_members(mate_paths[2], assembled_r2)
    if assembled_r1.stat().st_size != sum(path.stat().st_size for path in mate_paths[1]):
        raise SystemExit("[ERROR] assembled R1 byte count does not equal its gzip members")
    if assembled_r2.stat().st_size != sum(path.stat().st_size for path in mate_paths[2]):
        raise SystemExit("[ERROR] assembled R2 byte count does not equal its gzip members")
    write_input_provenance(input_provenance, args.sample_id, runs, assembled_r1, assembled_r2)
    print(f"[PASS] assembled {len(runs)} verified paired runs for {args.sample_id}")

    environment = os.environ.copy()
    environment.update({
        "SAMPLE_ID": args.sample_id,
        "RAW_R1": str(assembled_r1), "RAW_R2": str(assembled_r2),
        "SAMPLE_WORK": str(work), "RECEIPT": str(receipt),
        "TARGET_CONDITION": row["condition"], "STUDY": row["study"],
        "BATCH_ID": row.get("batch_id", ""), "BATCH_POSITION": row.get("batch_position", ""),
        "INPUT_PROVENANCE": str(input_provenance),
    })
    subprocess.run([str(args.runner.resolve())], check=True, env=environment)
    outputs = verify_receipt(receipt, work)
    verified_marker.write_text(
        f"outputs\t{len(outputs)}\nmanifest\t{args.manifest.resolve()}\n"
        f"input_provenance\t{input_provenance}\n",
        encoding="utf-8",
    )
    print(f"[OK] verified {len(outputs)} retained outputs for {args.sample_id}")
    if args.delete_inputs_after_verification:
        safe_cleanup(work)
    else:
        print(f"[KEEP] inputs retained: {work}")


if __name__ == "__main__":
    main()
