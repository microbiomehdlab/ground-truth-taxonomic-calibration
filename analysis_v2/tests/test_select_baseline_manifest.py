#!/usr/bin/env python3
"""Fixture tests for header-aware baseline manifest selection.

Replaces positional `awk` field indexing, which breaks silently if the canonical
column order changes.
"""

import csv
import subprocess
import tempfile
from pathlib import Path

from test_canonical_input import HEADER, row

SCRIPT = Path(__file__).resolve().parents[2] / "analysis_v2/scripts/select_baseline_manifest.py"
I = {name: idx for idx, name in enumerate(HEADER)}


def run(canonical, outdir, *extra):
    return subprocess.run(
        ["python3", str(SCRIPT), "--canonical", str(canonical), "--outdir", str(outdir), *extra],
        text=True, capture_output=True)


def read(outdir, name="baseline_manifest.tsv"):
    with (outdir / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="selman.") as name:
        root = Path(name)

        # Canonical table: both profilers, baseline and positive dose.
        rows = []
        for profiler, base, spike in [("kraken2_bracken", "0.2", "0.19"),
                                      ("metaphlan4", "1.0", "1.2")]:
            rows.append(row(profiler, f"spike-{profiler}", "0", base, "0"))
            rows.append(row(profiler, f"spike-{profiler}", "0.01", spike, "1000"))
        table = root / "canonical.tsv"; write(table, rows)

        out = root / "o1"
        result = run(table, out)
        assert result.returncode == 0, result.stderr
        selected = read(out)
        assert len(selected) == 1, "only the MetaPhlAn baseline should survive"
        assert selected[0]["profiler"] == "metaphlan4"
        assert selected[0]["spike_fraction_total"] == "0"
        assert selected[0]["profile_id"] == selected[0]["baseline_profile_id"]
        audit = (out / "baseline_manifest_audit.tsv").read_text()
        assert "dropped_other_profiler\t2" in audit
        assert "dropped_positive_dose\t1" in audit
        assert (out / "SUCCESS").is_file()

        # Selection is by name, so reordered columns still work.
        order = [HEADER.index(c) for c in reversed(HEADER)]
        shuffled = root / "shuffled.tsv"
        with shuffled.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(list(reversed(HEADER)))
            writer.writerows([[r[i] for i in order] for r in rows])
        out2 = root / "o2"
        assert run(shuffled, out2).returncode == 0
        assert read(out2) == selected, "column order must not change the selection"

        # include=0 baselines are dropped and counted.
        excluded = row("metaphlan4", "spike-metaphlan4", "0", "1.0", "0")
        excluded[I["sample_id"]] = "S2"
        excluded[I["include"]] = "0"
        excluded[I["exclusion_reason"]] = "qc_fail"
        table3 = root / "with_excluded.tsv"; write(table3, rows + [excluded])
        out3 = root / "o3"
        assert run(table3, out3).returncode == 0
        assert len(read(out3)) == 1
        assert "dropped_include_0\t1" in (out3 / "baseline_manifest_audit.tsv").read_text()

        # A differing analysis_population is EXPECTED, not a conflict: the same
        # physical baseline is written for both populations.
        other_population = row("metaphlan4", "spike-metaphlan4", "0", "1.0", "0")
        other_population[I["analysis_population"]] = "community"  # helper default: independent
        table4 = root / "both_populations.tsv"
        write(table4, rows + [other_population])
        out4 = root / "o4"
        assert run(table4, out4).returncode == 0
        collapsed = read(out4)
        assert len(collapsed) == 1, "both populations collapse to one sample-wide row"
        assert collapsed[0]["analysis_population"] == ""
        assert "samples_in_both_populations\t1" in \
            (out4 / "baseline_manifest_audit.tsv").read_text()

        # A differing physical source path IS a conflict.
        conflicting = row("metaphlan4", "spike-metaphlan4", "0", "1.0", "0")
        conflicting[I["analysis_population"]] = "community"
        conflicting[I["source_profile"]] = "profiles/other.metaphlan.tsv"
        table4b = root / "conflict.tsv"; write(table4b, rows + [conflicting])
        result = run(table4b, root / "o4b")
        assert result.returncode != 0 and "disagree on source_profile" in result.stderr

        # --expected-profiles enforcement, both directions.
        out4c = root / "o4c"
        assert run(table, out4c, "--expected-profiles", "1").returncode == 0
        assert (out4c / "SUCCESS").is_file()
        out4d = root / "o4d"
        result = run(table, out4d, "--expected-profiles", "2")
        assert result.returncode != 0 and "expected 2" in result.stderr
        assert not (out4d / "SUCCESS").exists()

        # No MetaPhlAn rows at all is a hard failure, not an empty manifest.
        table5 = root / "bracken_only.tsv"
        write(table5, [r for r in rows if r[I["profiler"]] == "kraken2_bracken"])
        result = run(table5, root / "o5")
        assert result.returncode != 0 and "no metaphlan4 zero-dose" in result.stderr

        # Missing required column fails loudly.
        table6 = root / "missing.tsv"
        keep = [i for i, c in enumerate(HEADER) if c != "profiler"]
        with table6.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([HEADER[i] for i in keep])
            writer.writerows([[r[i] for i in keep] for r in rows])
        result = run(table6, root / "o6")
        assert result.returncode != 0 and "missing column" in result.stderr

    print("[PASS] baseline manifest selection fixture tests")


if __name__ == "__main__":
    main()
