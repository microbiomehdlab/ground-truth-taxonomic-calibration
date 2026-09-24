#!/usr/bin/env python3
"""Audit and seal a completed Feng or Zeller upstream cohort."""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import Counter
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def marker_exists(path: Path) -> bool:
    """Return whether a sentinel exists; production creates these with touch."""
    return path.is_file()


def verify_sidecar(path: Path) -> None:
    sidecar = path.with_suffix(path.suffix + ".sha256")
    if not nonempty_file(sidecar):
        raise ValueError(f"missing manifest checksum: {sidecar}")
    fields = sidecar.read_text(encoding="utf-8").split()
    if len(fields) != 2 or fields[0] != digest(path) or fields[1] != path.name:
        raise ValueError(f"manifest checksum mismatch: {path}")


def verify_receipt(
    receipt: Path, scratch_root: Path, persistent_roots: tuple[Path, ...]
) -> tuple[bool, int, str]:
    """Rehash every retained output and reject disposable-scratch paths."""
    if not nonempty_file(receipt):
        return False, 0, "missing_or_empty_receipt"
    try:
        with receipt.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["path", "sha256", "bytes"]:
                return False, 0, "invalid_receipt_header"
            rows = list(reader)
        if not rows:
            return False, 0, "empty_receipt"
        scratch = scratch_root.resolve()
        allowed = tuple(root.resolve() for root in persistent_roots)
        seen: set[Path] = set()
        for index, row in enumerate(rows, start=2):
            raw_path = row.get("path", "")
            if not raw_path or not row.get("sha256") or not row.get("bytes"):
                return False, len(rows), f"malformed_receipt_row_{index}"
            path = Path(raw_path).resolve()
            if path in seen:
                return False, len(rows), f"duplicate_receipt_path_{index}"
            seen.add(path)
            if path == scratch or scratch in path.parents:
                return False, len(rows), f"retained_output_inside_scratch_{index}"
            if not any(path == root or root in path.parents for root in allowed):
                return False, len(rows), f"retained_output_outside_sample_roots_{index}"
            if not nonempty_file(path):
                return False, len(rows), f"missing_or_empty_output_{index}"
            try:
                expected_bytes = int(row["bytes"])
            except ValueError:
                return False, len(rows), f"invalid_bytes_{index}"
            if expected_bytes <= 0 or path.stat().st_size != expected_bytes:
                return False, len(rows), f"size_mismatch_{index}"
            if digest(path) != row["sha256"]:
                return False, len(rows), f"sha256_mismatch_{index}"
        return True, len(rows), ""
    except (OSError, csv.Error) as error:
        detail = str(error).replace("\t", " ").replace("\n", " ")
        return False, 0, f"receipt_read_error:{detail}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True, choices=("feng", "zeller"))
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--independent-manifest", required=True, type=Path)
    parser.add_argument("--state-dir", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--qc-root", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--expected-samples", required=True, type=int)
    parser.add_argument("--expected-independent", type=int, default=30)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for stale in (args.outdir / "SUCCESS", args.outdir / "production_seal.sha256"):
        stale.unlink(missing_ok=True)
    in_progress = args.outdir / "AUDIT_IN_PROGRESS"
    in_progress.write_text(
        f"cohort\t{args.cohort}\nstatus\tIN_PROGRESS\nmanifest\t{args.manifest.resolve()}\n",
        encoding="utf-8",
    )
    try:
        verify_sidecar(args.manifest); verify_sidecar(args.independent_manifest)
        manifest = read(args.manifest); independent_rows = read(args.independent_manifest)
        required = {"sample_id", "condition", "study", "age", "sex", "bmi", "independent_subset"}
        if not manifest or not required.issubset(manifest[0]):
            raise ValueError("production manifest lacks required metadata columns")
        if len(manifest) != args.expected_samples:
            raise ValueError(f"expected {args.expected_samples} samples; found {len(manifest)}")
        identifiers = [row["sample_id"] for row in manifest]
        independent = {row["sample_id"] for row in independent_rows}
        if len(set(identifiers)) != len(identifiers) or len(independent) != args.expected_independent or not independent.issubset(identifiers):
            raise ValueError("production/independent sample identity contract failed")
        if args.expected_independent == 30 and Counter(row["condition"] for row in independent_rows) != {"Control": 10, "Adenoma": 10, "CRC": 10}:
            raise ValueError("independent subset is not balanced 10/10/10")
        flow = []
        for row in manifest:
            sample = row["sample_id"]; sample_root = args.results_root / row["study"] / sample
            expected_profiles = 68 if sample in independent else 8
            observed_profiles = sum(
                1 for path in (sample_root / "profiles").rglob("SUCCESS")
                if marker_exists(path)
            )
            receipt_ok, receipt_files, receipt_error = verify_receipt(
                args.state_dir / "samples" / f"{sample}.retained_outputs.tsv",
                args.scratch_root / sample,
                (sample_root, args.qc_root / row["study"] / sample),
            )
            checks = {
                "verified_marker": nonempty_file(args.state_dir / "samples" / f"{sample}.verified"),
                "retained_output_receipt": receipt_ok,
                "input_provenance": nonempty_file(args.state_dir / "samples" / f"{sample}.input_provenance.tsv"),
                "sample_success": marker_exists(sample_root / "SUCCESS"),
                "profile_count": observed_profiles == expected_profiles,
            }
            flow.append({"sample_id": sample, "study": row["study"], "condition": row["condition"],
                         "independent_subset": int(sample in independent), "expected_profiles": expected_profiles,
                         "observed_profiles": observed_profiles, "retained_output_files": receipt_files,
                         **{name: int(ok) for name, ok in checks.items()},
                         "status": "PASS" if all(checks.values()) else "FAIL",
                         "failure_reasons": ";".join(name for name, ok in checks.items() if not ok),
                         "receipt_error": receipt_error})
        flow_path = args.outdir / "sample_flow.tsv"
        with flow_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(flow[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(flow)
        covariates = args.outdir / "covariate_audit.tsv"
        with covariates.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["field", "missing", "present", "distinct_nonmissing"])
            for field in ("condition", "age", "sex", "bmi"):
                values = [row[field].strip() for row in manifest]
                writer.writerow([field, sum(not value for value in values), sum(bool(value) for value in values),
                                 len({value for value in values if value})])
        for source, name in ((args.manifest, "production_manifest.tsv"),
                             (args.independent_manifest, "production_manifest.independent.tsv")):
            shutil.copyfile(source, args.outdir / name)
        if any(row["status"] != "PASS" for row in flow):
            raise ValueError(f"{sum(row['status'] != 'PASS' for row in flow)} samples are incomplete; see {flow_path}")
        success = args.outdir / "SUCCESS"
        success.write_text(f"cohort\t{args.cohort}\nsamples\t{len(manifest)}\nindependent_samples\t{len(independent)}\nstatus\tPASS\n", encoding="utf-8")
        artifacts = [flow_path, covariates, args.outdir / "production_manifest.tsv",
                     args.outdir / "production_manifest.independent.tsv", success]
        (args.outdir / "production_seal.sha256").write_text(
            "".join(f"{digest(path)}  {path.name}\n" for path in artifacts), encoding="utf-8")
        in_progress.unlink()
        print(f"[PASS] Sealed {args.cohort} upstream: {len(manifest)} samples")
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}") from error


if __name__ == "__main__":
    main()
