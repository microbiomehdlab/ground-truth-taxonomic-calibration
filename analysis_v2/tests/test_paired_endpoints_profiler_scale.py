#!/usr/bin/env python3
"""Fixture tests for the profiler-scale (genome-equivalent) paired reference.

Bracken reports a read fraction, so its reference is unchanged. MetaPhlAn
reports a genome-equivalent-like composition, so the implanted read fraction is
rescaled by G_eff/G_t and the composition renormalised. These tests pin that
behaviour and the fail-closed gates around it.
"""

import csv
import math
import subprocess
import tempfile
from pathlib import Path

from test_canonical_input import HEADER, row

I = {name: idx for idx, name in enumerate(HEADER)}
SCRIPT = Path(__file__).resolve().parents[2] / "analysis_v2/scripts/derive_paired_endpoints.py"


def member(profiler, profile, target, f_target, f_total, native, population="community"):
    pairs = "0" if f_target == "0" else "1000"
    entry = row(profiler, profile, f_target, native, pairs)
    entry[I["target_label"]] = target
    entry[I["analysis_population"]] = population
    entry[I["spike_fraction_total"]] = f_total
    entry[I["spike_fraction_target"]] = f_target
    return entry


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def sizes(path, mapping):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_label", "genome_size_bp"])
        writer.writerows(mapping.items())


def effective(path, value, cohort="yachida", sample="S1"):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cohort", "sample_id", "effective_genome_size_bp", "method"])
        writer.writerow([cohort, sample, value, "fixture"])


def run(table, outdir, *extra):
    return subprocess.run(
        ["python3", str(SCRIPT), "--input", str(table), "--outdir", str(outdir), *extra],
        text=True, capture_output=True)


def read_rows(outdir):
    with (outdir / "paired_endpoints.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="profiler_scale.") as name:
        root = Path(name)

        # --- 1. MetaPhlAn rows without an explicit reference choice fail. ----
        table = root / "gate.tsv"
        write(table, [member("metaphlan4", "p1", "Fnuc", "0", "0", "1.0"),
                      member("metaphlan4", "p1", "Fnuc", "0.01", "0.01", "1.2")])
        result = run(table, root / "gate")
        assert result.returncode != 0, "missing --metaphlan-reference must fail closed"
        assert "metaphlan-reference" in result.stderr

        # --- 2. genome_equivalent without genome tables fails. ---------------
        result = run(table, root / "gate2", "--metaphlan-reference", "genome_equivalent")
        assert result.returncode != 0 and "requires both" in result.stderr

        # --- 3. Unmapped target fails closed. --------------------------------
        gt, ge = root / "gt.tsv", root / "ge.tsv"
        sizes(gt, {"Pmic": 1661863}); effective(ge, 3_240_000)
        result = run(table, root / "gate3", "--metaphlan-reference", "genome_equivalent",
                     "--target-genome-sizes", str(gt), "--effective-genome-sizes", str(ge))
        assert result.returncode != 0 and "no genome size" in result.stderr

        # --- 4. Independent spike, G_t == G_eff  =>  identical to read scale. -
        # q = f * G_eff/G_t = f, and Q = q, so renormalisation is the only
        # difference; with equal sizes the profiler scale must still be a valid
        # reference and the response ratio must be 1 for an ideal observation.
        gt_equal = root / "gt_equal.tsv"
        sizes(gt_equal, {"Fnuc": 3_240_000})
        f, o = 0.01, 0.01
        q = f  # equal genome sizes
        expected = ((1 - f) * o + q) / ((1 - f) + q)
        table = root / "equal.tsv"
        write(table, [member("metaphlan4", "p1", "Fnuc", "0", "0", str(o * 100), "independent"),
                      member("metaphlan4", "p1", "Fnuc", str(f), str(f),
                             str(expected * 100), "independent")])
        out = root / "equal_out"
        result = run(table, out, "--metaphlan-reference", "genome_equivalent",
                     "--target-genome-sizes", str(gt_equal), "--effective-genome-sizes", str(ge))
        assert result.returncode == 0, result.stderr
        got = read_rows(out)[0]
        assert got["reference_type"] == "genome_equivalent"
        assert math.isclose(float(got["expected_abundance_profiler_scale"]), expected, rel_tol=1e-9)
        assert math.isclose(float(got["response_ratio_profiler_scale"]), 1.0, rel_tol=1e-9)
        assert math.isclose(float(got["response_residual_profiler_scale"]), 0.0, abs_tol=1e-12)
        # Read-proportional fields are preserved unchanged alongside.
        assert math.isclose(float(got["read_proportional_reference"]), (1 - f) * o + f, rel_tol=1e-9)

        # --- 5. Unequal genome sizes: a small genome must raise the reference.
        # G_t < G_eff  =>  q > f  =>  expected target abundance exceeds the
        # read-proportional value. This is the whole point of the correction.
        gt_small = root / "gt_small.tsv"
        sizes(gt_small, {"Dpne": 1_247_407})
        table = root / "small.tsv"
        write(table, [member("metaphlan4", "p1", "Dpne", "0", "0", "0.5", "independent"),
                      member("metaphlan4", "p1", "Dpne", "0.01", "0.01", "2.0", "independent")])
        out = root / "small_out"
        assert run(table, out, "--metaphlan-reference", "genome_equivalent",
                   "--target-genome-sizes", str(gt_small),
                   "--effective-genome-sizes", str(ge)).returncode == 0
        got = read_rows(out)[0]
        q_it = 0.01 * 3_240_000 / 1_247_407
        assert math.isclose(float(got["implanted_genome_equivalent_fraction_target"]), q_it, rel_tol=1e-9)
        assert float(got["expected_abundance_profiler_scale"]) > float(got["read_proportional_reference"])

        # --- 6. Community: Q sums over all members of the same profile. ------
        gt_comm = root / "gt_comm.tsv"
        sizes(gt_comm, {"Dpne": 1_247_407, "Hhat": 5_697_783})
        table = root / "comm.tsv"
        write(table, [
            member("metaphlan4", "pB", "Dpne", "0", "0", "0.5"),
            member("metaphlan4", "pB", "Hhat", "0", "0", "0.5"),
            member("metaphlan4", "pC", "Dpne", "0.005", "0.01", "1.0"),
            member("metaphlan4", "pC", "Hhat", "0.005", "0.01", "0.4"),
        ])
        out = root / "comm_out"
        assert run(table, out, "--metaphlan-reference", "genome_equivalent",
                   "--target-genome-sizes", str(gt_comm),
                   "--effective-genome-sizes", str(ge)).returncode == 0
        rows = {r["target_label"]: r for r in read_rows(out)}
        q_expected = 0.005 * 3_240_000 / 1_247_407 + 0.005 * 3_240_000 / 5_697_783
        for entry in rows.values():
            assert math.isclose(float(entry["implanted_genome_equivalent_fraction_total"]),
                                q_expected, rel_tol=1e-9), "Q must sum over all members"
        assert float(rows["Dpne"]["implanted_genome_equivalent_fraction_target"]) > \
               float(rows["Hhat"]["implanted_genome_equivalent_fraction_target"])

        # --- 7. Incomplete community membership fails closed. ---------------
        table = root / "partial.tsv"
        write(table, [
            member("metaphlan4", "pB", "Dpne", "0", "0", "0.5"),
            member("metaphlan4", "pC", "Dpne", "0.005", "0.01", "1.0"),
        ])
        result = run(table, root / "partial_out", "--metaphlan-reference", "genome_equivalent",
                     "--target-genome-sizes", str(gt_comm), "--effective-genome-sizes", str(ge))
        assert result.returncode != 0 and "missing community members" in result.stderr

        # --- 8. Zero baseline is valid, not an error. ------------------------
        table = root / "zero.tsv"
        write(table, [member("metaphlan4", "p1", "Dpne", "0", "0", "0.0", "independent"),
                      member("metaphlan4", "p1", "Dpne", "0.01", "0.01", "2.0", "independent")])
        out = root / "zero_out"
        assert run(table, out, "--metaphlan-reference", "genome_equivalent",
                   "--target-genome-sizes", str(gt_small),
                   "--effective-genome-sizes", str(ge)).returncode == 0
        got = read_rows(out)[0]
        assert math.isclose(float(got["retained_baseline_profiler_scale"]), 0.0, abs_tol=1e-15)

        # --- 9. Bracken is untouched by the MetaPhlAn choice. ---------------
        table = root / "bracken.tsv"
        write(table, [member("kraken2_bracken", "k1", "Fnuc", "0", "0", "0.2", "independent"),
                      member("kraken2_bracken", "k1", "Fnuc", "0.01", "0.01", "0.19", "independent")])
        out = root / "bracken_out"
        assert run(table, out).returncode == 0, "Bracken-only input needs no reference flag"
        got = read_rows(out)[0]
        assert got["reference_type"] == "read_proportional"
        assert math.isclose(float(got["expected_abundance_profiler_scale"]),
                            float(got["read_proportional_reference"]), rel_tol=1e-12)
        assert math.isclose(float(got["response_ratio_profiler_scale"]),
                            float(got["response_ratio"]), rel_tol=1e-12)

        # --- 10. read_proportional opt-out keeps MetaPhlAn on the old scale. -
        table = root / "optout.tsv"
        write(table, [member("metaphlan4", "p1", "Dpne", "0", "0", "0.5", "independent"),
                      member("metaphlan4", "p1", "Dpne", "0.01", "0.01", "2.0", "independent")])
        out = root / "optout_out"
        assert run(table, out, "--metaphlan-reference", "read_proportional").returncode == 0
        got = read_rows(out)[0]
        assert got["reference_type"] == "read_proportional"
        assert got["target_genome_size_bp"] == ""

    print("[PASS] profiler-scale paired-endpoint fixture tests")


if __name__ == "__main__":
    main()
