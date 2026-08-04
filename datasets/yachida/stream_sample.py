#!/usr/bin/env python3
"""Download, run, verify, and optionally clean one manifest sample.

The site-specific runner receives SAMPLE_ID, RAW_R1, RAW_R2, SAMPLE_WORK and
RECEIPT as environment variables.  On success it must write RECEIPT as a TSV
with header ``path<TAB>sha256<TAB>bytes`` and one row per retained output.
Only after every receipt entry has been independently verified can this driver
delete the sample's scratch directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import pathlib
import shutil
import subprocess
import urllib.request


SENTINEL = ".ground_truth_stream_sample"


def checksum(path: pathlib.Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: pathlib.Path, expected_md5: str, expected_bytes: int) -> None:
    if destination.is_file() and destination.stat().st_size == expected_bytes:
        if checksum(destination, "md5") == expected_md5:
            print(f"[SKIP] verified download: {destination}")
            return
    partial = destination.with_suffix(destination.suffix + ".partial")
    partial.unlink(missing_ok=True)
    with urllib.request.urlopen(url) as response, partial.open("wb") as output:
        while block := response.read(4 * 1024 * 1024):
            output.write(block)
    if partial.stat().st_size != expected_bytes or checksum(partial, "md5") != expected_md5:
        partial.unlink(missing_ok=True)
        raise SystemExit(f"[ERROR] size/MD5 verification failed: {url}")
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--scratch-root", type=pathlib.Path, required=True)
    parser.add_argument("--state-dir", type=pathlib.Path, required=True)
    parser.add_argument("--runner", type=pathlib.Path, required=True)
    parser.add_argument("--delete-inputs-after-verification", action="store_true")
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        matches = [row for row in csv.DictReader(handle, delimiter="\t") if row["sample_id"] == args.sample_id]
    if len(matches) != 1:
        raise SystemExit(f"[ERROR] expected one manifest row for {args.sample_id}, found {len(matches)}")
    row = matches[0]
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
    verified_marker = state_dir / f"{args.sample_id}.verified"
    if verified_marker.is_file() and receipt.is_file():
        outputs = verify_receipt(receipt, work)
        print(f"[SKIP] already verified {len(outputs)} retained outputs for {args.sample_id}")
        return

    work.mkdir(parents=True, exist_ok=True)
    (work / SENTINEL).write_text("managed disposable sample directory\n", encoding="utf-8")
    raw = work / "raw"
    raw.mkdir(exist_ok=True)
    r1 = raw / f"{args.sample_id}_1.fastq.gz"
    r2 = raw / f"{args.sample_id}_2.fastq.gz"
    download(row["fastq1_url"], r1, row["fastq1_md5"], int(row["fastq1_bytes"]))
    download(row["fastq2_url"], r2, row["fastq2_md5"], int(row["fastq2_bytes"]))

    environment = os.environ.copy()
    environment.update({
        "SAMPLE_ID": args.sample_id,
        "RAW_R1": str(r1),
        "RAW_R2": str(r2),
        "SAMPLE_WORK": str(work),
        "RECEIPT": str(receipt),
        "TARGET_CONDITION": row["Target_Condition"],
        "STUDY": row["Study"],
    })
    subprocess.run([str(args.runner.resolve())], check=True, env=environment)
    outputs = verify_receipt(receipt, work)
    verified_marker.write_text(
        f"outputs\t{len(outputs)}\nmanifest\t{args.manifest.resolve()}\n", encoding="utf-8"
    )
    print(f"[OK] verified {len(outputs)} retained outputs for {args.sample_id}")

    if args.delete_inputs_after_verification:
        if not (work / SENTINEL).is_file():
            raise SystemExit(f"[ERROR] cleanup sentinel missing: {work / SENTINEL}")
        shutil.rmtree(work)
        print(f"[CLEANED] deleted verified disposable inputs: {work}")
    else:
        print(f"[KEEP] inputs retained: {work}")


if __name__ == "__main__":
    main()
