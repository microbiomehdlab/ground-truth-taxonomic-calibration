#!/usr/bin/env python3
"""Tests for the runner-level MetaPhlAn reference helper.

Covers: no silent default, fail-closed on missing reference tables, separate
primary and sensitivity output directories, provenance capture, and Bracken-only
inputs remaining unaffected.
"""

import csv
import subprocess
import tempfile
from pathlib import Path

from test_canonical_input import HEADER, row

ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "analysis_v2/lib/metaphlan_reference.sh"
I = {name: idx for idx, name in enumerate(HEADER)}


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def pair(profiler, base, spike):
    return [row(profiler, f"spike-{profiler}", "0", base, "0"),
            row(profiler, f"spike-{profiler}", "0.01", spike, "1000")]


def run_helper(canonical, endpoint_root, env=None, sensitivity=True):
    script = (f'set -euo pipefail\n'
              f'ROOT="{ROOT}"\ncd "$ROOT"\n'
              f'source "{HELPER}"\n'
              f'derive_endpoints_with_references "{canonical}" "{endpoint_root}"\n')
    full = {"PATH": "/usr/bin:/bin:/usr/local/bin",
            "METAPHLAN_READ_REFERENCE_SENSITIVITY": "1" if sensitivity else "0"}
    full.update(env or {})
    return subprocess.run(["bash", "-c", script], text=True, capture_output=True, env=full)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="mpa_runner.") as name:
        root = Path(name)

        mpa = root / "mpa.tsv"; write(mpa, pair("metaphlan4", "1.0", "1.2"))
        brk = root / "brk.tsv"; write(brk, pair("kraken2_bracken", "0.2", "0.19"))

        targets = root / "targets.tsv"
        targets.write_text("target_label\tgenome_size_bp\nFnuc\t2180101\n", encoding="utf-8")
        geff = root / "geff.tsv"
        geff.write_text("cohort\tsample_id\teffective_genome_size_bp\nyachida\tS1\t3477000\n",
                        encoding="utf-8")

        # --- MetaPhlAn present but no reference tables -> fail closed --------
        result = run_helper(mpa, root / "o1")
        assert result.returncode != 0, "must not silently default"
        assert "TARGET_GENOME_SIZES" in result.stderr

        # Only one of the two supplied is still a failure.
        result = run_helper(mpa, root / "o2", {"TARGET_GENOME_SIZES": str(targets)})
        assert result.returncode != 0 and "EFFECTIVE_GENOME_SIZES" in result.stderr

        # An empty reference table is rejected.
        empty = root / "empty.tsv"; empty.write_text("", encoding="utf-8")
        result = run_helper(mpa, root / "o3", {"TARGET_GENOME_SIZES": str(targets),
                                               "EFFECTIVE_GENOME_SIZES": str(empty)})
        assert result.returncode != 0 and "missing or empty" in result.stderr

        # --- Happy path: primary + sensitivity in separate directories ------
        out = root / "endpoints"
        result = run_helper(mpa, out, {"TARGET_GENOME_SIZES": str(targets),
                                       "EFFECTIVE_GENOME_SIZES": str(geff)})
        assert result.returncode == 0, result.stderr
        sens = Path(str(out) + "_read_reference_sensitivity")
        assert (out / "SUCCESS").is_file(), "primary produced"
        assert (sens / "SUCCESS").is_file(), "sensitivity produced"
        assert out != sens and sens.is_dir(), "separate directories"
        assert (sens / "SENSITIVITY_ONLY.txt").is_file()

        with (out / "paired_endpoints.tsv").open(newline="") as handle:
            primary_rows = list(csv.DictReader(handle, delimiter="\t"))
        with (sens / "paired_endpoints.tsv").open(newline="") as handle:
            sens_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert primary_rows[0]["reference_type"] == "genome_equivalent"
        assert sens_rows[0]["reference_type"] == "read_proportional"
        assert primary_rows[0]["expected_abundance_profiler_scale"] != \
               sens_rows[0]["expected_abundance_profiler_scale"], \
            "the two references must give different expectations"
        # Read-proportional fields are preserved identically in both.
        assert primary_rows[0]["read_proportional_reference"] == \
               sens_rows[0]["read_proportional_reference"]

        prov = (out / "metaphlan_reference_provenance.tsv").read_text()
        assert "role\tPRIMARY" in prov and "genome_equivalent" in prov
        assert "target_genome_sizes_sha256" in prov
        assert "role\tSENSITIVITY" in (sens / "metaphlan_reference_provenance.tsv").read_text()

        # --- Re-running into an existing root refuses to overwrite ----------
        result = run_helper(mpa, out, {"TARGET_GENOME_SIZES": str(targets),
                                       "EFFECTIVE_GENOME_SIZES": str(geff)})
        assert result.returncode != 0 and "already exists" in result.stderr

        # --- Sensitivity can be disabled explicitly -------------------------
        out2 = root / "endpoints_nosens"
        result = run_helper(mpa, out2, {"TARGET_GENOME_SIZES": str(targets),
                                        "EFFECTIVE_GENOME_SIZES": str(geff)},
                            sensitivity=False)
        assert result.returncode == 0, result.stderr
        assert not Path(str(out2) + "_read_reference_sensitivity").exists()

        # --- Bracken-only input needs no reference tables and is unchanged --
        out3 = root / "endpoints_bracken"
        result = run_helper(brk, out3)
        assert result.returncode == 0, result.stderr
        with (out3 / "paired_endpoints.tsv").open(newline="") as handle:
            bracken_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert bracken_rows[0]["reference_type"] == "read_proportional"
        assert bracken_rows[0]["expected_abundance_profiler_scale"] == \
               bracken_rows[0]["read_proportional_reference"]

        # --- Hardening gates -------------------------------------------------
        # A canonical table with no profiler column is malformed, not Bracken-only.
        noprof = root / "noprof.tsv"
        rows = pair("metaphlan4", "1.0", "1.2")
        keep = [i for i, c in enumerate(HEADER) if c != "profiler"]
        with noprof.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([HEADER[i] for i in keep])
            writer.writerows([[r[i] for i in keep] for r in rows])
        result = run_helper(noprof, root / "o_noprof")
        assert result.returncode != 0 and "no 'profiler' column" in result.stderr

        # Reference-table schema is validated, not just non-emptiness.
        bad_schema = root / "bad_schema.tsv"
        bad_schema.write_text("label\tsize\nFnuc\t2180101\n", encoding="utf-8")
        result = run_helper(mpa, root / "o_schema",
                            {"TARGET_GENOME_SIZES": str(bad_schema),
                             "EFFECTIVE_GENOME_SIZES": str(geff)})
        assert result.returncode != 0 and "missing column" in result.stderr

        # A MetaPhlAn sample absent from the G_eff table fails before running.
        other = root / "geff_other.tsv"
        other.write_text("cohort\tsample_id\teffective_genome_size_bp\n"
                         "yachida\tOTHER\t3477000\n", encoding="utf-8")
        result = run_helper(mpa, root / "o_cover",
                            {"TARGET_GENOME_SIZES": str(targets),
                             "EFFECTIVE_GENOME_SIZES": str(other)})
        assert result.returncode != 0 and "absent from the G_eff table" in result.stderr

        # Disabling the sensitivity is recorded in the primary provenance.
        out4 = root / "endpoints_rec"
        assert run_helper(mpa, out4, {"TARGET_GENOME_SIZES": str(targets),
                                      "EFFECTIVE_GENOME_SIZES": str(geff)},
                          sensitivity=False).returncode == 0
        assert "sensitivity_disabled" in \
            (out4 / "metaphlan_reference_provenance.tsv").read_text()
        # Canonical checksum is recorded.
        assert "canonical_input_sha256" in \
            (out4 / "metaphlan_reference_provenance.tsv").read_text()

    print("[PASS] MetaPhlAn reference runner helper tests")


if __name__ == "__main__":
    main()
