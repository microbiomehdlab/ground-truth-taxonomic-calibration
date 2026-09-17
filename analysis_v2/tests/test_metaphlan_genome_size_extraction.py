#!/usr/bin/env python3
"""Fixture tests for MetaPhlAn database genome-size extraction.

The real database is only present on the cluster, so these tests use synthetic
pickles that mimic the schemas MetaPhlAn releases have used. They pin the
schema-tolerant lookup, the refusal to treat marker records as a size source,
and the fail-closed behaviour on an unrecognised schema.
"""

import bz2
import csv
import pickle
import subprocess
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis_v2/scripts/extract_metaphlan_genome_sizes.py"
# Real vJan25 taxonomy keys terminate in t__SGB..., confirmed by
# scripts/build_reference_representation_table.py::metaphlan_sgb_map.
A = "k__Bacteria|p__Firmicutes|g__G|s__Aaa|t__SGB1"
B = "k__Bacteria|p__Firmicutes|g__G|s__Bbb|t__SGB2"


def write_db(path: Path, obj, compress=False):
    if compress:
        with bz2.open(path, "wb") as handle:
            pickle.dump(obj, handle)
    else:
        path.write_bytes(pickle.dumps(obj))


def run(db, outdir, *extra):
    return subprocess.run(
        ["python3", str(SCRIPT), "--database", str(db), "--outdir", str(outdir), *extra],
        text=True, capture_output=True)


def read(outdir):
    with (outdir / "metaphlan_genome_sizes.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="mpa_db.") as name:
        root = Path(name)

        # --- 1. (taxid, genome_length) tuples, the MetaPhlAn 3/4 layout -----
        db = root / "v4.pkl"
        write_db(db, {"taxonomy": {A: ("2|1239|11|22", 2_000_000),
                                   B: ("2|1239|11|23", 6_000_000)},
                      "markers": {"m1": {"clade": "s__Aaa", "len": 1200}}})
        out = root / "o1"
        result = run(db, out)
        assert result.returncode == 0, result.stderr
        rows = read(out)
        assert [r["clade_lineage"] for r in rows] == [A, B], "output must be sorted"
        assert rows[0]["genome_size_bp"] == "2000000"
        assert rows[0]["ncbi_tax_id"] == "2|1239|11|22"
        assert rows[0]["rank"] == "t", "vJan25 keys are SGB-terminal"
        assert rows[0]["clade_leaf"] == "t__SGB1"
        # Marker length must never leak in as a genome size.
        assert all(r["genome_size_bp"] != "1200" for r in rows)
        assert (out / "SUCCESS").is_file()
        prov = (out / "genome_size_provenance.tsv").read_text()
        assert "database_sha256" in prov and "extraction_command" in prov

        # --- 2. Bare-length layout is also accepted -------------------------
        db2 = root / "bare.pkl"
        write_db(db2, {"taxonomy": {A: 2_000_000, B: 6_000_000}})
        out2 = root / "o2"
        assert run(db2, out2).returncode == 0
        assert read(out2)[0]["genome_size_bp"] == "2000000"

        # --- 3. bz2-compressed database is read transparently ---------------
        db3 = root / "comp.pkl.bz2"
        write_db(db3, {"taxonomy": {A: ("2", 2_000_000)}}, compress=True)
        out3 = root / "o3"
        assert run(db3, out3).returncode == 0
        assert len(read(out3)) == 1

        # --- 4. Implausible lengths are dropped, not rounded into range -----
        db4 = root / "implausible.pkl"
        write_db(db4, {"taxonomy": {A: ("2", 2_000_000), B: ("2", 900)}})
        out4 = root / "o4"
        assert run(db4, out4).returncode == 0
        rows4 = read(out4)
        assert [r["clade_lineage"] for r in rows4] == [A]
        assert "skipped_without_usable_length\t1" in \
            (out4 / "genome_size_provenance.tsv").read_text()

        # --- 5. Unrecognised schema fails closed and reports what it saw ----
        db5 = root / "unknown.pkl"
        write_db(db5, {"something_else": {"x": 1}})
        result = run(db5, root / "o5")
        assert result.returncode != 0
        assert "no taxonomy mapping" in result.stderr
        assert "observed schema" in result.stderr, "must print schema for a human to adapt"

        # --- 6. --inspect-only writes nothing -------------------------------
        out6 = root / "o6"
        result = run(db, out6, "--inspect-only")
        assert result.returncode == 0 and "top-level type" in result.stdout
        assert not (out6 / "metaphlan_genome_sizes.tsv").exists()

    print("[PASS] metaphlan genome-size extraction fixture tests")


if __name__ == "__main__":
    main()
