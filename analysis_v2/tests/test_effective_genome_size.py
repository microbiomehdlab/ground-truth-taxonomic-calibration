#!/usr/bin/env python3
"""Fixture tests for sample-specific effective community genome size G_eff,i.

Fixtures mimic the real vJan25 structure: database lineages terminate in
`t__SGB...`, and MetaPhlAn profiles carry the full hierarchy with both species
and terminal SGB rows.
"""

import csv
import math
import subprocess
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis_v2/scripts/compute_effective_genome_size.py"
MAN_FIELDS = ["cohort", "sample_id", "analysis_population", "profiler", "profile_id",
              "baseline_profile_id", "spike_fraction_total", "source_profile",
              "include", "exclusion_reason"]
K = "k__Bacteria"
P = f"{K}|p__Firmicutes"
G = f"{P}|c__C|o__O|f__F|g__G"


def sp(name):
    return f"{G}|s__{name}"


def sgb(name, ident):
    return f"{sp(name)}|t__SGB{ident}"


def profile(path, species_rows, sgb_rows, extra=()):
    """Hierarchical profile: parents, species rows, terminal SGB rows."""
    lines = ["#mpa_vJan25_CHOCOPhlAnSGB", "#clade_name\tNCBI_tax_id\trelative_abundance"]
    total = sum(v for _, v in species_rows)
    for clade in (K, P, G):
        lines.append(f"{clade}\t2\t{total}")
    for name, value in species_rows:
        lines.append(f"{sp(name)}\t2|1239\t{value}")
    for lineage, value in sgb_rows:
        lines.append(f"{lineage}\t2|1239\t{value}")
    for clade, value in extra:
        lines.append(f"{clade}\t-1\t{value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def manifest(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MAN_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def base_row(sample, source, **over):
    row = {"cohort": "yachida", "sample_id": sample, "analysis_population": "community",
           "profiler": "metaphlan4", "profile_id": f"base-{sample}",
           "baseline_profile_id": f"base-{sample}", "spike_fraction_total": "0",
           "source_profile": str(source), "include": "1", "exclusion_reason": ""}
    row.update(over)
    return row


def sizes_table(path, entries, column="clade_lineage"):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([column, "clade_leaf", "rank", "ncbi_tax_id", "genome_size_bp"])
        for key, size in entries:
            leaf = key.split("|")[-1]
            writer.writerow([key, leaf, leaf.split("__", 1)[0], "0", size])


def run(man, sizes, outdir, *extra):
    return subprocess.run(
        ["python3", str(SCRIPT), "--manifest", str(man), "--genome-sizes", str(sizes),
         "--outdir", str(outdir), *extra], text=True, capture_output=True)


def read(outdir, name="effective_genome_size.tsv"):
    with (outdir / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="geff.") as root_name:
        root = Path(root_name)
        SGB_A, SGB_B = sgb("Aaa", 1), sgb("Bbb", 2)
        sgb_sizes = root / "sgb_sizes.tsv"
        sizes_table(sgb_sizes, [(SGB_A, 2_000_000), (SGB_B, 6_000_000)])

        # --- 1. Weighted G_eff from terminal SGB rows; no double counting ----
        # Species rows carry the same 75/25 mass as their SGB children. If
        # species mass leaked in, G_eff would still be 3 Mb but eligible mass
        # would be 200, so the assertion on eligible_abundance is the real test.
        prof = root / "p1.tsv"
        profile(prof, [("Aaa", 75.0), ("Bbb", 25.0)],
                [(SGB_A, 75.0), (SGB_B, 25.0)])
        man = root / "m1.tsv"; manifest(man, [base_row("S1", prof)])
        out = root / "o1"
        result = run(man, sgb_sizes, out)
        assert result.returncode == 0, result.stderr
        got = read(out)[0]
        assert math.isclose(float(got["effective_genome_size_bp"]), 3_000_000, rel_tol=1e-9)
        assert math.isclose(float(got["eligible_abundance"]), 100.0, rel_tol=1e-9), \
            "species rows must not be added to their own SGB mass"
        assert got["n_mapped_features"] == "2" and got["profile_rank"] == "sgb"
        assert math.isclose(float(got["mapping_coverage"]), 1.0, rel_tol=1e-12)

        # --- 2. Fractions give the same G_eff as percentages ----------------
        prof2 = root / "p2.tsv"
        profile(prof2, [("Aaa", 0.75), ("Bbb", 0.25)],
                [(SGB_A, 0.75), (SGB_B, 0.25)])
        man2 = root / "m2.tsv"; manifest(man2, [base_row("S2", prof2)])
        out2 = root / "o2"
        assert run(man2, sgb_sizes, out2).returncode == 0
        got2 = read(out2)[0]
        assert got2["abundance_unit"] == "fraction"
        assert math.isclose(float(got2["effective_genome_size_bp"]), 3_000_000, rel_tol=1e-9)

        # --- 3. Rank incompatibility fails closed ---------------------------
        species_sizes = root / "species_sizes.tsv"
        sizes_table(species_sizes, [(sp("Aaa"), 2_000_000), (sp("Bbb"), 6_000_000)])
        result = run(man, species_sizes, root / "o3")
        assert result.returncode != 0, "species-only mapping must not serve sgb mode"
        assert "no identifiers of rank 't'" in result.stderr

        result = run(man, sgb_sizes, root / "o3b", "--profile-rank", "species")
        assert result.returncode != 0 and "no identifiers of rank 's'" in result.stderr

        # --- 4. species mode works with a species mapping -------------------
        out4 = root / "o4"
        assert run(man, species_sizes, out4, "--profile-rank", "species").returncode == 0
        got4 = read(out4)[0]
        assert math.isclose(float(got4["effective_genome_size_bp"]), 3_000_000, rel_tol=1e-9)
        assert got4["profile_rank"] == "species"

        # --- 5. Terminal species row without an SGB child stays eligible ----
        # It cannot match an SGB key, so it is unmapped and lowers coverage
        # rather than being silently dropped.
        prof5 = root / "p5.tsv"
        profile(prof5, [("Aaa", 70.0), ("Bbb", 20.0), ("Ccc", 10.0)],
                [(SGB_A, 70.0), (SGB_B, 20.0)])
        man5 = root / "m5.tsv"; manifest(man5, [base_row("S5", prof5)])
        out5 = root / "o5"
        assert run(man5, sgb_sizes, out5, "--min-coverage", "0.9").returncode == 0
        got5 = read(out5)[0]
        assert math.isclose(float(got5["mapping_coverage"]), 0.9, rel_tol=1e-9)
        assert math.isclose(float(got5["effective_genome_size_bp"]),
                            (70 * 2e6 + 20 * 6e6) / 90, rel_tol=1e-9)
        unmapped = read(out5, "unmapped_features.tsv")
        assert any(r["feature"].endswith("s__Ccc") for r in unmapped)

        # --- 6. Coverage below threshold excludes and fails ------------------
        result = run(man5, sgb_sizes, root / "o6", "--min-coverage", "0.95")
        assert result.returncode != 0
        assert read(root / "o6", "excluded_profiles.tsv")[0]["reason"] == "coverage_below_threshold"
        assert not (root / "o6" / "SUCCESS").exists()

        # --- 7. Exact SGB matching: a different SGB id does not match -------
        wrong = root / "wrong_sgb.tsv"
        sizes_table(wrong, [(sgb("Aaa", 99), 2_000_000), (SGB_B, 6_000_000)])
        result = run(man, wrong, root / "o7", "--min-coverage", "0.95")
        assert result.returncode != 0, "SGB1 must not match SGB99"

        # --- 8. Conflicting and invalid sizes -------------------------------
        dup_ok = root / "dup_ok.tsv"
        sizes_table(dup_ok, [(SGB_A, 2_000_000), (SGB_A, 2_000_000), (SGB_B, 6_000_000)])
        assert run(man, dup_ok, root / "o8").returncode == 0
        dup_bad = root / "dup_bad.tsv"
        sizes_table(dup_bad, [(SGB_A, 2_000_000), (SGB_A, 3_000_000), (SGB_B, 6_000_000)])
        result = run(man, dup_bad, root / "o8b")
        assert result.returncode != 0 and "conflicting sizes" in result.stderr
        bad = root / "bad.tsv"
        sizes_table(bad, [(SGB_A, 0), (SGB_B, 6_000_000)])
        result = run(man, bad, root / "o8c")
        assert result.returncode != 0 and "outside" in result.stderr

        # --- 9. UNCLASSIFIED and non-species terminal mass is excluded ------
        prof9 = root / "p9.tsv"
        profile(prof9, [("Aaa", 75.0), ("Bbb", 25.0)], [(SGB_A, 75.0), (SGB_B, 25.0)],
                extra=[("UNCLASSIFIED", 12.0), (f"{K}|p__Lonely", 5.0)])
        man9 = root / "m9.tsv"; manifest(man9, [base_row("S9", prof9)])
        out9 = root / "o9"
        assert run(man9, sgb_sizes, out9).returncode == 0
        got9 = read(out9)[0]
        assert math.isclose(float(got9["eligible_abundance"]), 100.0, rel_tol=1e-9)
        assert math.isclose(float(got9["excluded_abundance"]), 17.0, rel_tol=1e-9)

        # --- 10. Zero eligible abundance --------------------------------------
        prof10 = root / "p10.tsv"
        prof10.write_text("#h\nUNCLASSIFIED\t-1\t100.0\n", encoding="utf-8")
        man10 = root / "m10.tsv"; manifest(man10, [base_row("S10", prof10)])
        result = run(man10, sgb_sizes, root / "o10")
        assert result.returncode != 0
        assert read(root / "o10", "excluded_profiles.tsv")[0]["reason"] == "zero_eligible_abundance"

        # --- 11. Profiler enforcement: mixed and Bracken-only both fail -----
        man11 = root / "m11.tsv"
        manifest(man11, [base_row("S1", prof),
                         base_row("S1", prof, profiler="kraken2_bracken",
                                  profile_id="kbase-S1", baseline_profile_id="kbase-S1")])
        result = run(man11, sgb_sizes, root / "o11")
        assert result.returncode != 0 and "metaphlan4 rows only" in result.stderr
        man11b = root / "m11b.tsv"
        manifest(man11b, [base_row("S1", prof, profiler="kraken2_bracken")])
        result = run(man11b, sgb_sizes, root / "o11b")
        assert result.returncode != 0 and "metaphlan4 rows only" in result.stderr

        # --- 12. Inclusion policy --------------------------------------------
        man12 = root / "m12.tsv"
        manifest(man12, [base_row("S1", prof),
                         base_row("S12", prof, include="0", exclusion_reason="qc_fail")])
        result = run(man12, sgb_sizes, root / "o12")
        assert result.returncode != 0 and "require_included_only" in result.stderr
        out12 = root / "o12b"
        assert run(man12, sgb_sizes, out12, "--inclusion-policy", "filter_excluded",
                   "--max-excluded-samples", "1").returncode == 0
        assert any(r["reason"] == "canonical_include_0"
                   for r in read(out12, "excluded_profiles.tsv"))
        man12c = root / "m12c.tsv"
        manifest(man12c, [base_row("S1", prof, include="0", exclusion_reason="")])
        result = run(man12c, sgb_sizes, root / "o12c", "--inclusion-policy", "filter_excluded")
        assert result.returncode != 0 and "requires exclusion_reason" in result.stderr

        # --- Repeated baseline rows: repetition ok, conflicts fail ----------
        man13 = root / "m13.tsv"
        manifest(man13, [base_row("S1", prof), base_row("S1", prof)])
        out13 = root / "o13"
        assert run(man13, sgb_sizes, out13).returncode == 0
        assert len(read(out13)) == 1, "repeated baselines collapse to one row"

        # The same physical baseline repeated across populations is EXPECTED
        # and must collapse to one sample-wide row with blank population.
        man_pop = root / "m_pop.tsv"
        manifest(man_pop, [base_row("S1", prof, analysis_population="community"),
                           base_row("S1", prof, analysis_population="independent")])
        out_pop = root / "o_pop"
        assert run(man_pop, sgb_sizes, out_pop).returncode == 0
        pop_rows = read(out_pop)
        assert len(pop_rows) == 1, "one physical baseline yields one sample-wide G_eff"
        assert pop_rows[0]["analysis_population"] == "", \
            "blank population drives the generic downstream lookup"

        # Different sample ids are different samples, not a conflict.
        man_two = root / "m_two.tsv"
        manifest(man_two, [base_row("S1", prof), base_row("S2", prof2)])
        out_two = root / "o_two"
        assert run(man_two, sgb_sizes, out_two).returncode == 0
        assert len(read(out_two)) == 2

        # Conflicts on physical-baseline fields must fail.
        for field, value in [("source_profile", str(prof2)),
                             ("profiler", "kraken2_bracken"),
                             ("baseline_profile_id", "other-base")]:
            bad_man = root / f"m13_{field}.tsv"
            conflicting = base_row("S1", prof)
            conflicting[field] = value
            manifest(bad_man, [base_row("S1", prof), conflicting])
            result = run(bad_man, sgb_sizes, root / f"o13_{field}")
            # Some conflicts trip an earlier, more specific gate: a mismatched
            # baseline_profile_id breaks the zero-dose identity invariant, and a
            # foreign profiler is rejected before dedup.
            expected = {"profiler": "metaphlan4 rows only",
                        "baseline_profile_id": "post-spike"}.get(field, "disagree on")
            assert result.returncode != 0 and expected in result.stderr, field

        # Two cohorts may reuse a sample id; they stay distinct. Each cohort
        # has its own physical profile, as in production.
        feng_prof = root / "feng_S1.tsv"
        profile(feng_prof, [("Aaa", 75.0), ("Bbb", 25.0)], [(SGB_A, 75.0), (SGB_B, 25.0)])
        man_ch = root / "m_cohorts.tsv"
        manifest(man_ch, [base_row("S1", prof),
                          base_row("S1", feng_prof, cohort="feng")])
        out_ch = root / "o_cohorts"
        assert run(man_ch, sgb_sizes, out_ch).returncode == 0
        assert {r["cohort"] for r in read(out_ch)} == {"yachida", "feng"}

        # One physical file claimed by two different samples must fail.
        man_own = root / "m_owner.tsv"
        manifest(man_own, [base_row("S1", prof), base_row("S9", prof)])
        result = run(man_own, sgb_sizes, root / "o_owner")
        assert result.returncode != 0 and "claimed by both" in result.stderr

        # --- 14. Post-spike rejection ----------------------------------------
        man14 = root / "m14.tsv"
        manifest(man14, [base_row("S1", prof, spike_fraction_total="0.01",
                                  profile_id="spiked-S1")])
        result = run(man14, sgb_sizes, root / "o14")
        assert result.returncode != 0 and "post-spike" in result.stderr

        # --- 15. Post-spike abundances cannot change G_eff -------------------
        spiked = root / "p1_spiked.tsv"
        profile(spiked, [("Aaa", 1.0), ("Bbb", 99.0)], [(SGB_A, 1.0), (SGB_B, 99.0)])
        out15 = root / "o15"
        assert run(man, sgb_sizes, out15).returncode == 0
        assert read(out15)[0]["effective_genome_size_bp"] == got["effective_genome_size_bp"]

        # --- 16. Source-profile checksum ledger, deduplicated ----------------
        ledger = read(out13, "source_profile_checksums.tsv")
        assert len(ledger) == 1, "one physical file must be hashed once"
        assert ledger[0]["sha256"] and int(ledger[0]["bytes"]) > 0
        assert "source_profile_checksums.tsv" in (out13 / "effective_genome_size.sha256").read_text()

        # --- 17. Deterministic, order-independent output ---------------------
        rows17 = [base_row("S5", prof5), base_row("S1", prof)]
        man17 = root / "m17.tsv"
        manifest(man17, rows17)
        a = root / "o17a"
        assert run(man17, sgb_sizes, a, "--min-coverage", "0.9").returncode == 0
        manifest(man17, list(reversed(rows17)))
        b = root / "o17b"
        assert run(man17, sgb_sizes, b, "--min-coverage", "0.9").returncode == 0
        assert (a / "effective_genome_size.tsv").read_text() == \
               (b / "effective_genome_size.tsv").read_text()
        assert [r["sample_id"] for r in read(a)] == ["S1", "S5"]

    print("[PASS] effective-genome-size fixture tests")


if __name__ == "__main__":
    main()
