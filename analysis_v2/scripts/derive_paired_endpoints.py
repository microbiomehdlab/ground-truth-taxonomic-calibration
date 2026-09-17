#!/usr/bin/env python3
"""Derive baseline-adjusted paired endpoints from canonical v2 input."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import sys
from pathlib import Path


KEY_FIELDS = [
    "cohort", "sample_id", "analysis_population", "target_label", "profiler",
    "assembly_arm",
]
OUTPUT_FIELDS = [
    "schema_version", "cohort", "study", "sample_id", "condition",
    "analysis_population", "target_label", "target_taxon", "assembly_arm",
    "profiler", "profile_id", "baseline_profile_id", "spike_fraction_total",
    "spike_fraction_target", "implanted_read_pairs_target", "native_unit",
    "baseline_native_abundance", "observed_native_abundance",
    "baseline_abundance_fraction", "observed_abundance_fraction",
    "baseline_retained_after_dilution", "read_proportional_reference",
    "recovered_spike_signal", "response_ratio", "signed_reference_error",
    "absolute_reference_error", "baseline_detected_native_nonzero",
    "observed_detected_native_nonzero", "source_baseline_profile",
    "source_profile", "source_design",
    # Profiler-scale reference (additive; read-proportional fields above are
    # retained unchanged as the labelled sensitivity representation).
    "reference_type", "target_genome_size_bp",
    "effective_community_genome_size_bp",
    "implanted_genome_equivalent_fraction_target",
    "implanted_genome_equivalent_fraction_total",
    "retained_baseline_profiler_scale", "expected_abundance_profiler_scale",
    "implanted_signal_profiler_scale", "recovered_spike_signal_profiler_scale",
    "response_ratio_profiler_scale", "response_residual_profiler_scale",
]

GENOME_EQUIVALENT_PROFILERS = {"metaphlan4"}
# sum_k f_ik must reproduce F_i for a community perturbation. A wider gap means
# the canonical table does not contain every implanted member, so Q_i would be
# understated and the correction silently wrong.
COMMUNITY_FRACTION_TOLERANCE = 0.01


def read_table(path: Path, required: list[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no data rows")
    missing = [field for field in required if field not in rows[0]]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
    return rows


def positive(value: str, path: Path, label: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{path}: {label} is not numeric") from error
    if not number > 0:
        raise ValueError(f"{path}: {label} must be > 0")
    return number


def load_target_genome_sizes(path: Path) -> dict[str, float]:
    rows = read_table(path, ["target_label", "genome_size_bp"])
    sizes: dict[str, float] = {}
    for entry in rows:
        label = entry["target_label"]
        if label in sizes:
            raise ValueError(f"{path}: duplicate target_label {label}")
        sizes[label] = positive(entry["genome_size_bp"], path, f"genome_size_bp for {label}")
    return sizes


def load_effective_genome_sizes(path: Path) -> dict[tuple[str, str, str], float]:
    rows = read_table(path, ["cohort", "sample_id", "effective_genome_size_bp"])
    sizes: dict[tuple[str, str, str], float] = {}
    for entry in rows:
        scope = (entry["cohort"], entry["sample_id"],
                 entry.get("analysis_population", "") or "")
        if scope in sizes:
            raise ValueError(f"{path}: duplicate scope {scope}")
        sizes[scope] = positive(
            entry["effective_genome_size_bp"], path,
            f"effective_genome_size_bp for {scope}")
    return sizes


def effective_size(sizes: dict[tuple[str, str, str], float],
                   row: dict[str, str]) -> float | None:
    """Prefer a population-specific entry, else a sample-wide one."""
    specific = (row["cohort"], row["sample_id"], row["analysis_population"])
    generic = (row["cohort"], row["sample_id"], "")
    if specific in sizes:
        return sizes[specific]
    return sizes.get(generic)


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field] for field in KEY_FIELDS)


def render(value: float) -> str:
    return format(value, ".17g")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--metaphlan-reference", choices=["genome_equivalent", "read_proportional"],
        help="Quantitative reference scale for MetaPhlAn rows. Required whenever "
             "the input contains MetaPhlAn rows; there is no default.")
    parser.add_argument(
        "--target-genome-sizes", type=Path,
        help="TSV of target_label, genome_size_bp measured from the exact spike FASTA.")
    parser.add_argument(
        "--effective-genome-sizes", type=Path,
        help="TSV of cohort, sample_id, [analysis_population], effective_genome_size_bp, "
             "derived independently of post-spike recovery.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    validator = Path(__file__).with_name("validate_canonical_input.py")
    validation_dir = args.outdir / "input_validation"
    subprocess.run([
        sys.executable, str(validator), "--input", str(args.input),
        "--outdir", str(validation_dir),
    ], check=True)

    with args.input.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    included = [row for row in rows if row["include"] == "1"]
    excluded = [row for row in rows if row["include"] == "0"]
    baselines = {
        key(row): row for row in included
        if float(row["spike_fraction_target"]) == 0
    }

    # --- profiler-scale reference gate -------------------------------------
    # Fail closed: if MetaPhlAn rows are present the operator must state which
    # reference scale applies. Silently defaulting would reintroduce the
    # read-fraction-versus-genome-equivalent estimand mismatch.
    ge_profilers = {r["profiler"] for r in included} & GENOME_EQUIVALENT_PROFILERS
    if ge_profilers and args.metaphlan_reference is None:
        raise SystemExit(
            "[FAIL] input contains MetaPhlAn rows; --metaphlan-reference "
            "{genome_equivalent|read_proportional} is required")
    use_ge = bool(ge_profilers) and args.metaphlan_reference == "genome_equivalent"

    target_sizes: dict[str, float] = {}
    eff_sizes: dict[tuple[str, str, str], float] = {}
    genome_inputs: list[Path] = []
    if use_ge:
        if not (args.target_genome_sizes and args.effective_genome_sizes):
            raise SystemExit(
                "[FAIL] --metaphlan-reference genome_equivalent requires both "
                "--target-genome-sizes and --effective-genome-sizes")
        target_sizes = load_target_genome_sizes(args.target_genome_sizes)
        eff_sizes = load_effective_genome_sizes(args.effective_genome_sizes)
        genome_inputs = [args.target_genome_sizes, args.effective_genome_sizes]

        ge_rows = [r for r in included if r["profiler"] in GENOME_EQUIVALENT_PROFILERS]
        missing_t = sorted({r["target_label"] for r in ge_rows
                            if r["target_label"] not in target_sizes})
        if missing_t:
            raise SystemExit(f"[FAIL] no genome size for target(s): {', '.join(missing_t)}")
        missing_e = sorted({(r["cohort"], r["sample_id"]) for r in ge_rows
                            if effective_size(eff_sizes, r) is None})
        if missing_e:
            raise SystemExit(
                "[FAIL] no effective community genome size for sample(s): "
                + ", ".join(f"{c}/{s}" for c, s in missing_e))

    # Q_i is a property of one perturbed profile, so group by profile_id: an
    # independent spike yields one member, a community spike yields all ten.
    q_total: dict[str, float] = {}
    if use_ge:
        members: dict[str, list[dict[str, str]]] = {}
        for row in included:
            if row["profiler"] not in GENOME_EQUIVALENT_PROFILERS:
                continue
            if float(row["spike_fraction_target"]) == 0:
                continue
            members.setdefault(row["profile_id"], []).append(row)
        for profile_id, rows_for_profile in members.items():
            first = rows_for_profile[0]
            g_eff = effective_size(eff_sizes, first)
            total = 0.0
            fraction_sum = 0.0
            for member in rows_for_profile:
                f_k = float(member["spike_fraction_target"])
                fraction_sum += f_k
                total += f_k * g_eff / target_sizes[member["target_label"]]
            declared = float(first["spike_fraction_total"])
            if first["analysis_population"] == "community" and declared > 0:
                if abs(fraction_sum - declared) / declared > COMMUNITY_FRACTION_TOLERANCE:
                    raise SystemExit(
                        f"[FAIL] profile {profile_id}: implanted member fractions sum to "
                        f"{fraction_sum:.6g} but spike_fraction_total is {declared:.6g}; "
                        "the canonical table is missing community members, so Q_i "
                        "would be understated")
            denominator = (1 - declared) + total
            if denominator <= 0:
                raise SystemExit(f"[FAIL] profile {profile_id}: non-positive renormalisation denominator")
            q_total[profile_id] = total

    endpoints: list[dict[str, str]] = []
    for row in included:
        target_fraction = float(row["spike_fraction_target"])
        if target_fraction == 0:
            continue
        baseline = baselines[key(row)]
        total_fraction = float(row["spike_fraction_total"])
        baseline_fraction = float(baseline["abundance_fraction"])
        observed_fraction = float(row["abundance_fraction"])
        retained = (1 - total_fraction) * baseline_fraction
        reference = retained + target_fraction
        recovered = observed_fraction - retained
        ratio = recovered / target_fraction
        signed_error = observed_fraction - reference

        endpoint = {
            field: row[field] for field in (
                "schema_version", "cohort", "study", "sample_id", "condition",
                "analysis_population", "target_label", "target_taxon",
                "assembly_arm", "profiler", "profile_id", "baseline_profile_id",
                "spike_fraction_total", "spike_fraction_target",
                "implanted_read_pairs_target", "native_unit", "source_profile",
                "source_design",
            )
        }
        endpoint.update({
            "baseline_native_abundance": baseline["native_abundance"],
            "observed_native_abundance": row["native_abundance"],
            "baseline_abundance_fraction": render(baseline_fraction),
            "observed_abundance_fraction": render(observed_fraction),
            "baseline_retained_after_dilution": render(retained),
            "read_proportional_reference": render(reference),
            "recovered_spike_signal": render(recovered),
            "response_ratio": render(ratio),
            "signed_reference_error": render(signed_error),
            "absolute_reference_error": render(abs(signed_error)),
            "baseline_detected_native_nonzero": baseline["detected_native_nonzero"],
            "observed_detected_native_nonzero": row["detected_native_nonzero"],
            "source_baseline_profile": baseline["source_profile"],
        })

        # Profiler-scale reference. Bracken reports a read fraction, so it is
        # already on the perturbation's own scale. MetaPhlAn reports a
        # genome-equivalent-like composition, so the implanted read fraction is
        # rescaled by G_eff/G_t and the whole composition is renormalised.
        if use_ge and row["profiler"] in GENOME_EQUIVALENT_PROFILERS:
            g_t = target_sizes[row["target_label"]]
            g_eff = effective_size(eff_sizes, row)
            q_it = target_fraction * g_eff / g_t
            Q_i = q_total[row["profile_id"]]
            denominator = (1 - total_fraction) + Q_i
            retained_ps = (1 - total_fraction) * baseline_fraction / denominator
            expected_ps = ((1 - total_fraction) * baseline_fraction + q_it) / denominator
            signal_ps = q_it / denominator
            recovered_ps = observed_fraction - retained_ps
            endpoint.update({
                "reference_type": "genome_equivalent",
                "target_genome_size_bp": render(g_t),
                "effective_community_genome_size_bp": render(g_eff),
                "implanted_genome_equivalent_fraction_target": render(q_it),
                "implanted_genome_equivalent_fraction_total": render(Q_i),
                "retained_baseline_profiler_scale": render(retained_ps),
                "expected_abundance_profiler_scale": render(expected_ps),
                "implanted_signal_profiler_scale": render(signal_ps),
                "recovered_spike_signal_profiler_scale": render(recovered_ps),
                "response_ratio_profiler_scale": render(recovered_ps / signal_ps),
                "response_residual_profiler_scale": render(observed_fraction - expected_ps),
            })
        else:
            endpoint.update({
                "reference_type": "read_proportional",
                "target_genome_size_bp": "",
                "effective_community_genome_size_bp": "",
                "implanted_genome_equivalent_fraction_target": "",
                "implanted_genome_equivalent_fraction_total": "",
                "retained_baseline_profiler_scale": render(retained),
                "expected_abundance_profiler_scale": render(reference),
                "implanted_signal_profiler_scale": render(target_fraction),
                "recovered_spike_signal_profiler_scale": render(recovered),
                "response_ratio_profiler_scale": render(ratio),
                "response_residual_profiler_scale": render(signed_error),
            })
        endpoints.append(endpoint)

    output = args.outdir / "paired_endpoints.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(endpoints)

    exclusions = args.outdir / "excluded_canonical_rows.tsv"
    with exclusions.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(excluded)

    summary = args.outdir / "endpoint_derivation_summary.tsv"
    summary.write_text(
        "metric\tvalue\n"
        f"canonical_rows\t{len(rows)}\n"
        f"included_baseline_rows\t{len(baselines)}\n"
        f"derived_positive_dose_rows\t{len(endpoints)}\n"
        f"excluded_rows\t{len(excluded)}\n"
        f"metaphlan_reference\t{args.metaphlan_reference or 'not_applicable'}\n"
        f"genome_equivalent_rows\t{sum(1 for e in endpoints if e['reference_type'] == 'genome_equivalent')}\n"
        "status\tPASS\n",
        encoding="utf-8",
    )
    checksum = args.outdir / "endpoint_derivation.sha256"
    paths = [args.input, output, exclusions, summary,
             validation_dir / "canonical_input_validation.tsv"] + genome_inputs
    checksum.write_text(
        "".join(f"{sha256(path)}  {path.resolve()}\n" for path in paths),
        encoding="utf-8",
    )
    (args.outdir / "SUCCESS").write_text(
        f"endpoint_rows\t{len(endpoints)}\nstatus\tPASS\n", encoding="utf-8"
    )
    print(f"[PASS] Derived {len(endpoints)} paired positive-dose endpoints")
    print(f"[INFO] {output}")


if __name__ == "__main__":
    main()
