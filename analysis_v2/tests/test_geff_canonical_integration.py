#!/usr/bin/env python3
"""End-to-end integration: canonical input -> baseline manifest -> G_eff ->
paired endpoints on the MetaPhlAn genome-equivalent scale.

The fixture reproduces what `build_crc_cohort_canonical_input.py` actually
writes, so it cannot silently drift from the producer contract:

- two biological samples, both in the community population;
- one of them additionally in the independent population;
- community and independent zero-dose rows share `profile_id == sample_id`,
  `baseline_profile_id == sample_id` and the *same* native profile file,
  differing only in `analysis_population`;
- baselines repeat once per implanted target;
- both profilers present;
- positive-dose rows present for both populations.

The column list is imported from the canonical-input test helper rather than
restated, so a schema change breaks this test loudly.
"""

import csv
import math
import subprocess
import tempfile
from pathlib import Path

from test_canonical_input import HEADER

SCRIPTS = Path(__file__).resolve().parents[2] / "analysis_v2/scripts"
SELECT = SCRIPTS / "select_baseline_manifest.py"
GEFF = SCRIPTS / "compute_effective_genome_size.py"
ENDPOINTS = SCRIPTS / "derive_paired_endpoints.py"

TARGETS = [("Fnuc", "Fusobacterium nucleatum"), ("Dpne", "Dialister pneumosintes")]
G = "k__Bacteria|p__Firmicutes|c__C|o__O|f__F|g__G"
SGB_A = f"{G}|s__Aaa|t__SGB1"
SGB_B = f"{G}|s__Bbb|t__SGB2"
TARGET_GENOMES = {"Fnuc": 2_180_101, "Dpne": 1_247_407}


def canonical_row(**kw):
    """Build one canonical row with the producer's own defaults."""
    row = {
        "schema_version": "paired-dose-response-v2.1", "cohort": "yachida",
        "study": "YachidaS_2019", "sample_id": "S1", "condition": "CRC",
        "analysis_population": "community", "target_label": "Fnuc",
        "target_taxon": "Fusobacterium nucleatum", "assembly_arm": "original",
        "profiler": "metaphlan4", "profile_id": "S1", "baseline_profile_id": "S1",
        "spike_fraction_total": "0", "spike_fraction_target": "0",
        "implanted_read_pairs_target": "0", "native_abundance": "1.0",
        "native_unit": "relative_abundance_pct", "abundance_fraction": "0.01",
        "detected_native_nonzero": "1", "source_profile": "", "source_design": "BASELINE",
        "include": "1", "exclusion_reason": "",
    }
    row.update(kw)
    # Mirror make_row(): abundance_fraction is a pure unit conversion, which
    # validate_canonical_input.py enforces.
    native = float(row["native_abundance"])
    if row["profiler"] == "kraken2_bracken":
        row["native_unit"] = "fraction_total_reads"
        fraction = native
    else:
        row["native_unit"] = "relative_abundance_pct"
        fraction = native / 100.0
    row["abundance_fraction"] = format(fraction, ".17g")
    row["detected_native_nonzero"] = "1" if fraction > 0 else "0"
    return [row[c] for c in HEADER]


def write_profile(path, rows):
    lines = ["#mpa_vJan25_CHOCOPhlAnSGB", "#clade_name\tNCBI_tax_id\trelative_abundance"]
    total = sum(v for _, v in rows)
    for clade in ("k__Bacteria", "k__Bacteria|p__Firmicutes", G):
        lines.append(f"{clade}\t2\t{total}")
    for lineage, value in rows:
        species = lineage.rsplit("|t__", 1)[0]
        lines.append(f"{species}\t2\t{value}")      # parent species row
        lines.append(f"{lineage}\t2\t{value}")      # terminal SGB row
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_canonical(root, samples, independent_samples, baselines):
    """Mirror the producer: baseline row per target/profiler/population."""
    rows = []
    for sample in samples:
        for label, taxon in TARGETS:
            for profiler in ("kraken2_bracken", "metaphlan4"):
                suffix = ".bracken.S.tsv" if profiler == "kraken2_bracken" else ".metaphlan.tsv"
                base_file = baselines[(sample, profiler)]
                populations = ["community"]
                if sample in independent_samples:
                    populations.append("independent")
                for population in populations:
                    rows.append(canonical_row(
                        sample_id=sample, target_label=label, target_taxon=taxon,
                        profiler=profiler, profile_id=sample, baseline_profile_id=sample,
                        analysis_population=population, source_profile=str(base_file)))
                    # one positive dose per population
                    spiked = root / f"{sample}_{label}_{population}{suffix}"
                    spiked.write_text("x", encoding="utf-8")
                    rows.append(canonical_row(
                        sample_id=sample, target_label=label, target_taxon=taxon,
                        profiler=profiler, profile_id=f"{sample}_{label}_{population}",
                        baseline_profile_id=sample, analysis_population=population,
                        spike_fraction_total="0.01", spike_fraction_target="0.01",
                        implanted_read_pairs_target="1000", native_abundance="2.0",
                        source_profile=str(spiked), source_design="design.tsv"))
    return rows


def write_table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(script, *args):
    return subprocess.run(["python3", str(script), *map(str, args)],
                          text=True, capture_output=True)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="geff_integration.") as name:
        root = Path(name)
        samples = ["S1", "S2"]
        independent = {"S1"}

        baselines = {}
        for sample in samples:
            mpa = root / f"{sample}.metaphlan.tsv"
            write_profile(mpa, [(SGB_A, 75.0), (SGB_B, 25.0)])
            baselines[(sample, "metaphlan4")] = mpa
            brk = root / f"{sample}.bracken.S.tsv"
            brk.write_text("name\tfraction_total_reads\nAaa\t0.5\n", encoding="utf-8")
            baselines[(sample, "kraken2_bracken")] = brk

        canonical = root / "canonical_input.tsv"
        write_table(canonical, build_canonical(root, samples, independent, baselines))

        sizes = root / "sgb_sizes.tsv"
        with sizes.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["clade_lineage", "clade_leaf", "rank", "ncbi_tax_id", "genome_size_bp"])
            writer.writerow([SGB_A, "t__SGB1", "t", "0", 2_000_000])
            writer.writerow([SGB_B, "t__SGB2", "t", "0", 6_000_000])

        # --- Step 1: selection ------------------------------------------------
        sel = root / "sel"
        result = run(SELECT, "--canonical", canonical, "--cohort", "yachida",
                     "--profiler", "metaphlan4", "--expected-profiles", 2, "--outdir", sel)
        assert result.returncode == 0, result.stderr
        manifest_rows = read_tsv(sel / "baseline_manifest.tsv")
        assert len(manifest_rows) == 2, "one sample-wide baseline per biological sample"
        assert {r["sample_id"] for r in manifest_rows} == {"S1", "S2"}
        assert all(r["analysis_population"] == "" for r in manifest_rows), \
            "analysis_population must be blank to drive the generic lookup"
        assert all(r["profiler"] == "metaphlan4" for r in manifest_rows)
        audit = (sel / "baseline_manifest_audit.tsv").read_text()
        assert "samples_in_both_populations\t1" in audit, "S1 is in both populations"
        assert (sel / "SUCCESS").is_file()

        # Wrong expected count fails and withholds SUCCESS.
        bad = root / "sel_bad"
        result = run(SELECT, "--canonical", canonical, "--cohort", "yachida",
                     "--expected-profiles", 3, "--outdir", bad)
        assert result.returncode != 0 and "expected 3" in result.stderr
        assert not (bad / "SUCCESS").exists()

        # --- Step 2: G_eff ----------------------------------------------------
        geff_dir = root / "geff"
        result = run(GEFF, "--manifest", sel / "baseline_manifest.tsv",
                     "--genome-sizes", sizes, "--outdir", geff_dir)
        assert result.returncode == 0, result.stderr
        geff_rows = read_tsv(geff_dir / "effective_genome_size.tsv")
        assert len(geff_rows) == 2
        assert all(r["analysis_population"] == "" for r in geff_rows)
        for row in geff_rows:
            assert math.isclose(float(row["effective_genome_size_bp"]), 3_000_000, rel_tol=1e-9)
            assert math.isclose(float(row["eligible_abundance"]), 100.0, rel_tol=1e-9), \
                "species rows must not be counted alongside their SGB children"
        ledger = read_tsv(geff_dir / "source_profile_checksums.tsv")
        assert len(ledger) == 2, "each physical baseline hashed exactly once"

        # The uncollapsed canonical table itself must also work, collapsing the
        # community/independent repetitions to one row per sample.
        direct = root / "geff_direct"
        result = run(GEFF, "--manifest", canonical, "--genome-sizes", sizes,
                     "--outdir", direct)
        assert result.returncode != 0, "canonical table still carries Bracken rows"
        assert "metaphlan4 rows only" in result.stderr

        mpa_only = root / "canonical_mpa.tsv"
        write_table(mpa_only, [r for r in build_canonical(root, samples, independent, baselines)
                               if r[HEADER.index("profiler")] == "metaphlan4"
                               and r[HEADER.index("spike_fraction_total")] == "0"])
        direct2 = root / "geff_direct2"
        result = run(GEFF, "--manifest", mpa_only, "--genome-sizes", sizes, "--outdir", direct2)
        assert result.returncode == 0, result.stderr
        assert len(read_tsv(direct2 / "effective_genome_size.tsv")) == 2, \
            "community and independent repetitions collapse to one row per sample"

        # --- Step 3: paired endpoints consume the generic G_eff ---------------
        targets = root / "target_sizes.tsv"
        with targets.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["target_label", "genome_size_bp"])
            for label, size in TARGET_GENOMES.items():
                writer.writerow([label, size])

        mpa_canonical = root / "canonical_mpa_full.tsv"
        write_table(mpa_canonical, [r for r in build_canonical(root, samples, independent, baselines)
                                    if r[HEADER.index("profiler")] == "metaphlan4"])
        ep = root / "endpoints"
        result = run(ENDPOINTS, "--input", mpa_canonical, "--outdir", ep,
                     "--metaphlan-reference", "genome_equivalent",
                     "--target-genome-sizes", targets,
                     "--effective-genome-sizes", geff_dir / "effective_genome_size.tsv")
        assert result.returncode == 0, result.stderr
        endpoint_rows = read_tsv(ep / "paired_endpoints.tsv")
        assert endpoint_rows, "endpoints were produced"

        by_population = {}
        for row in endpoint_rows:
            assert row["reference_type"] == "genome_equivalent"
            by_population.setdefault(row["analysis_population"], set()).add(
                row["effective_community_genome_size_bp"])
        assert set(by_population) == {"community", "independent"}
        community = by_population["community"]
        independent_vals = by_population["independent"]
        assert community == independent_vals == {"3000000"}, \
            "both populations must resolve the same sample-wide G_eff"

        # Same target, same dose, two populations -> identical q_it.
        q_by_pop = {}
        for row in endpoint_rows:
            if row["sample_id"] == "S1" and row["target_label"] == "Fnuc":
                q_by_pop[row["analysis_population"]] = \
                    row["implanted_genome_equivalent_fraction_target"]
        assert q_by_pop["community"] == q_by_pop["independent"]

        # --- Step 4: conflicting physical baseline fails ----------------------
        rows = build_canonical(root, samples, independent, baselines)
        idx_pop = HEADER.index("analysis_population")
        idx_src = HEADER.index("source_profile")
        idx_dose = HEADER.index("spike_fraction_total")
        for row in rows:
            if row[idx_pop] == "independent" and row[idx_dose] == "0":
                row[idx_src] = str(root / "S2.metaphlan.tsv")   # wrong physical file
        conflict = root / "conflict.tsv"; write_table(conflict, rows)
        result = run(SELECT, "--canonical", conflict, "--cohort", "yachida", "--outdir", root / "sel_conflict")
        assert result.returncode != 0 and "disagree on source_profile" in result.stderr

        # --- Step 5: same profile_id in two cohorts stays distinct ------------
        two_cohort = build_canonical(root, ["S1"], set(), baselines)
        idx_cohort = HEADER.index("cohort")
        for row in build_canonical(root, ["S1"], set(), baselines):
            clone = list(row); clone[idx_cohort] = "feng"
            two_cohort.append(clone)
        table = root / "two_cohort.tsv"; write_table(table, two_cohort)
        out = root / "sel_two_cohort"
        result = run(SELECT, "--canonical", table, "--profiler", "metaphlan4",
                     "--expected-profiles", 2, "--outdir", out)
        assert result.returncode == 0, result.stderr
        cohorts = {r["cohort"] for r in read_tsv(out / "baseline_manifest.tsv")}
        assert cohorts == {"yachida", "feng"}, "profile_id S1 must not collide across cohorts"

        # --- Step 6: include=0 without a reason fails -------------------------
        rows = build_canonical(root, ["S1"], set(), baselines)
        idx_inc = HEADER.index("include")
        rows[0][idx_inc] = "0"
        table = root / "no_reason.tsv"; write_table(table, rows)
        result = run(SELECT, "--canonical", table, "--outdir", root / "sel_noreason")
        assert result.returncode != 0 and "exclusion_reason" in result.stderr

        # include=0 with a reason is dropped and appears in the ledger.
        rows[0][HEADER.index("exclusion_reason")] = "qc_fail"
        table = root / "with_reason.tsv"; write_table(table, rows)
        out = root / "sel_reason"
        assert run(SELECT, "--canonical", table, "--outdir", out).returncode == 0
        ledger_rows = read_tsv(out / "excluded_canonical_rows.tsv")
        assert ledger_rows and ledger_rows[0]["reason"] == "qc_fail"
        assert ledger_rows[0]["sample_id"] == "S1" and ledger_rows[0]["line"]

    print("[PASS] G_eff canonical integration tests")


if __name__ == "__main__":
    main()
