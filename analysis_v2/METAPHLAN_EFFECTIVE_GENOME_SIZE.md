# Effective community genome size for the MetaPhlAn profiler-scale reference

**Status: PROVISIONAL.** The estimator, fail-closed policy, and tests are
implemented and fixture-tested. The production MetaPhlAn database was **not**
available in the environment where this was written, so the authoritative
genome-size mapping has not been extracted and the sealed Yachida coverage
audit has not been run. The coverage thresholds below are therefore a stated
starting policy, not a frozen decision. This becomes `DECIDED` only after the
cluster steps in "Required cluster execution" succeed.

## Estimand

For each **unspiked baseline** MetaPhlAn profile `i`:

```
G_eff,i = sum_j(c_ij * G_j) / sum_j(c_ij)
```

- `c_ij` — baseline MetaPhlAn relative abundance of eligible feature `j`;
- `G_j` — the database representative genome length for feature `j`;
- the sum runs over eligible features that map to a genome size.

`G_eff,i` is the abundance-weighted mean representative genome length of the
community that MetaPhlAn reports for that sample. It is a community property of
the unperturbed sample, not a property of any spike.

It supports the MetaPhlAn profiler-scale reference in
`analysis_v2/ENDPOINTS.md` and `CLAUDE.md`:

```
q_it = f_it * G_eff,i / G_t
Q_i  = sum_k(f_ik * G_eff,i / G_k)
```

## Scope: one physical baseline, one sample-wide G_eff

`build_crc_cohort_canonical_input.py` writes the **same physical** zero-dose
MetaPhlAn profile for both analytical populations. The community and the
independent baseline row are both produced by
`make_row(meta, population, target, profiler, sample, sample, 0.0, 0.0, 0, ...)`
against the same `unique_file(baseline_root, sample + ".metaphlan.tsv")`, so they
share `profile_id == sample_id`, `baseline_profile_id == sample_id` and one
`source_profile`, differing only in `analysis_population`.

`G_eff,i` is therefore a property of that single physical baseline, not a
per-population estimate:

- identity is `(cohort, sample_id)` — never bare `profile_id`, which equals the
  sample id and would collide across cohorts;
- exactly one G_eff row per biological sample;
- `analysis_population` is written **empty** in the G_eff table. This is
  deliberate: `derive_paired_endpoints.py` looks up
  `(cohort, sample_id, analysis_population)` first and falls back to
  `(cohort, sample_id, "")`, so the blank value is what lets one sample-wide
  value serve both community and independent endpoints;
- the physical baseline is hashed exactly once;
- neither population label is arbitrarily retained.

Repetition across targets and across populations is expected and collapses
silently. Disagreement on any physical field — cohort, sample_id, profiler,
profile_id, baseline_profile_id, source_profile, either spike fraction, or
inclusion state — fails and names the field. A resolved source path claimed by
two different `(cohort, sample_id)` identities also fails rather than being
absorbed by `setdefault`.

`tests/test_geff_canonical_integration.py` runs the whole chain — canonical
fixture, selector, G_eff, paired endpoints — and asserts that the community and
independent endpoint rows resolve the same `G_eff` and the same `q_it`.

## Rank and identifier — resolved from repository evidence

The production MetaPhlAn vJan25 taxonomy is **SGB-based**: database lineages
terminate in `t__SGB...`. Evidence:

- `scripts/build_reference_representation_table.py::metaphlan_sgb_map()` reads
  `database["taxonomy"]` and requires **both** `s__` and `t__SGB\d+` in the same
  key, i.e. keys run `...|s__Species|t__SGB123`.
- `workflows/metaphlan4/postprocess_local.sh` exists to strip `t__` rows out of
  the raw profile ("Extract species-only … WITHOUT any `|t__`"), which it would
  not need to do unless raw profiles contained them.
- `source_profile` in the canonical table points at the **raw** `<sample>.metaphlan.tsv`
  (`build_crc_cohort_canonical_input.py`), not the species-only derivative, so
  the full hierarchy including SGB rows is available.

**Default `--profile-rank sgb`** therefore weights **terminal** rows at or below
species rank, keyed on the full lineage, which matches database keys exactly.

**Why this does not double-count.** MetaPhlAn repeats community mass at every
rank. Eligibility is restricted to *terminal* rows — rows with no descendant row
in the same profile — so a species row that has an SGB child is non-terminal and
never contributes alongside it. Terminal rows at or below species rank form a
non-overlapping partition of species-level mass. A test asserts
`eligible_abundance == 100` for a profile whose species and SGB rows each carry
the same 100 units; if species mass leaked in it would read 200.

A terminal `s__` row with **no** SGB child remains eligible. It cannot match an
SGB key, so it is reported as unmapped and lowers `mapping_coverage` rather than
being silently dropped.

`--profile-rank species` is the documented alternative: it weights every `s__`
row, matching the canonical input's species-level estimand, and requires a
species-rank mapping table. It would need an explicit SGB-to-species aggregation
policy on the mapping side, which is **not** implemented; no aggregation is
invented silently.

**Rank incompatibility fails closed.** The estimator inspects the mapping
table's ranks and refuses to run if it holds no identifier of the rank the
selected mode needs — an SGB-only mapping cannot serve the species estimator,
and a species-only mapping cannot serve the SGB estimator. Both directions are
tested. This is what prevents a silent zero-coverage run against the real
database.

## Authoritative inputs

| Quantity | Source | Notes |
|---|---|---|
| `G_t` (implanted targets) | exact spike FASTA lengths + SHA-256 in `CLAUDE.md` | measured from the assemblies used to simulate reads |
| `G_j` (community features) | the exact MetaPhlAn vJan25 database used for production | extracted by `scripts/extract_metaphlan_genome_sizes.py` |
| `c_ij` | raw `<sample>.metaphlan.tsv` baseline profiles | zero-dose, MetaPhlAn-only rows |

UHGG and generic NCBI lengths are **not** primary; they may be a labelled
sensitivity. Marker length is never a size source — the extractor reads the
database taxonomy mapping, and a test asserts a marker record cannot enter it.

The extractor records database path, filename, SHA-256, byte size, extraction
command, entry counts, and per-rank clade counts.

## Input enforcement

**Selection.** `scripts/select_baseline_manifest.py` keeps rows matching the
requested cohort with `profiler == metaphlan4`, `spike_fraction_total == 0`,
`spike_fraction_target == 0`, `profile_id == baseline_profile_id` and
`include == 1`; collapses them on `(cohort, sample_id)`; blanks
`analysis_population`; writes a per-row exclusion ledger
(`excluded_canonical_rows.tsv`) rather than only a count; requires a non-empty
`exclusion_reason` on every `include == 0` row even though that row is dropped;
and with `--expected-profiles N` withholds `SUCCESS` and exits non-zero unless
exactly N sample-wide baselines are selected.

**MetaPhlAn only.** `profiler` is a required manifest column. Every included row
must be `profiler == metaphlan4`; mixed-profiler and Bracken-only manifests fail
with a clear error. The check lives **in the script**, not only in the selection
command, so a hand-built manifest cannot bypass it. This prevents a Bracken
table being parsed as a MetaPhlAn profile.

**Inclusion policy.** `include` and `exclusion_reason` are required. The default
`--inclusion-policy require_included_only` **fails** if any `include == 0` row is
present, so excluded evidence cannot enter silently. `filter_excluded` is the
documented alternative: it drops them while writing each one to the exclusion
audit with its canonical reason. `include == 0` without a reason always fails.

**Repeated baseline rows.** A baseline legitimately repeats once per implanted
target. Rows sharing a `profile_id` are collapsed **only after** verifying they
agree on cohort, sample_id, analysis_population, profiler, profile_id,
baseline_profile_id, source_profile, spike fraction, and inclusion state. Any
disagreement fails and names the offending field. Conflicting sample IDs, source
paths, and population metadata are each tested.

**Identifier mapping** is exact: no fuzzy matching, no prefix matching, no name
normalisation, no taxon-name collapsing. Unmapped features are never imputed.

## Numerator, denominator, and coverage

- numerator `sum_j(c_ij * G_j)` over **mapped** eligible features;
- denominator `sum_j(c_ij)` over the same mapped features;
- `eligible_abundance` — all eligible mass, mapped or not;
- `mapping_coverage = mapped_abundance / eligible_abundance`.

Unmapped eligible mass leaves both numerator and denominator rather than being
imputed; coverage states how much of the community the estimate rests on.

Abundance units are detected per profile (`--abundance-unit auto`): an eligible
total above 10 is a percentage, otherwise a fraction. The estimator is a ratio
and so is invariant to this scaling; the detected unit is recorded per profile.

## Fail-closed behaviour

No override flag exists for any of these:

- **post-spike input** — the manifest must declare `spike_fraction_total == 0`
  and `profile_id == baseline_profile_id`, using the canonical contract's own
  invariants rather than a filename heuristic;
- non-MetaPhlAn rows;
- `include == 0` rows under the default policy, or without a reason;
- repeated rows that disagree on an invariant field;
- a mapping table with no identifier of the required rank;
- conflicting genome sizes for one identifier (agreeing duplicates are fine);
- sizes outside [100 kb, 50 Mb], or non-numeric;
- missing/malformed profile, non-numeric or negative abundance;
- zero eligible abundance, or no mapped features;
- coverage below threshold;
- more excluded profiles than `--max-excluded-samples` (default 0).

## Provenance

`source_profile_checksums.tsv` records cohort, sample_id, profile_id, resolved
path, byte size, and SHA-256 for every native profile that contributed to a
`G_eff` value. Each physical file is hashed **once**, after repeated-row
validation, so a baseline appearing once per target is not hashed ten times. The
ledger is included in `effective_genome_size.sha256` alongside the manifest,
mapping table, and every output.

## Diagnostics

`effective_genome_size.tsv` carries, per profile: `G_eff`, coverage, eligible /
mapped / unmapped / excluded abundance, mapped and unmapped feature counts, the
detected unit, and the source profile path.

`effective_genome_size_audit.tsv` carries counts, the coverage min/median/max,
the `G_eff` min/median/max, and the distinct unmapped-feature count.
`unmapped_features.tsv` ranks every unmapped identifier by how many profiles it
affects and how much abundance it carries — this is the table that shows whether
a coverage failure is one systematic gap or a long tail.

All outputs are sorted and order-independent; a test builds the same manifest in
two orders and asserts byte-identical output.

## Prespecified sensitivity alternatives

1. 90% coverage threshold instead of 95%.
2. Cohort-level median `G_eff` instead of sample-specific.
3. A single global `G_eff` constant.
4. UHGG v2.0.2 representative lengths as `G_j`.
5. `--match-on leaf` instead of lineage.

These are **labelled sensitivities**. The primary estimator is sample-specific,
and none of these may be selected by comparing recovery performance.

## Independence from spike outcomes

- Only zero-dose rows are readable; post-spike manifests are rejected outright.
- Nothing in either script reads observed post-spike abundance, recovery,
  response ratios, or any endpoint output.
- Test 11 of `tests/test_effective_genome_size.py` writes a spiked profile with
  a very different composition on disk and asserts `G_eff,i` is byte-identical
  to the baseline-only result.
- The provisional 3.10 Mb value fitted from development recovery data is **not**
  used anywhere. It was derived from post-spike outcomes and is therefore
  precisely what this policy forbids; it remains a diagnostic that the
  genome-size mechanism is real, nothing more.

## Limitations

- MetaPhlAn normalises by **marker** coverage within an SGB, not by whole-genome
  length. Whole-genome scaling is a defensible first-order approximation of that
  behaviour, not an exact reproduction of the estimator.
- Residual mechanisms remain: mappability, marker selection and copy number,
  strain divergence from the SGB representative, and database representation.
- Unmapped community mass biases `G_eff,i` toward the genome sizes of mapped
  taxa. The coverage threshold bounds, but does not eliminate, this.
- The database representative genome length is itself an estimate for
  uncultivated SGBs.
- `G_eff,i` is estimated from a compositional profile, so it inherits any
  compositional bias already present in the baseline.

## Required cluster execution

These require the production database and sealed profiles. Until they succeed
this policy stays **PROVISIONAL** and no definitive MetaPhlAn quantitative claim
is licensed. All outputs are `DEVELOPMENT_ONLY`.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
source config/yachida.strict-production.env

: "${METAPHLAN_DB:?METAPHLAN_DB is not set}"
: "${METAPHLAN_INDEX:?METAPHLAN_INDEX is not set}"
DB_PKL="$METAPHLAN_DB/$METAPHLAN_INDEX.pkl"
test -s "$DB_PKL" || { echo "missing MetaPhlAn database pickle: $DB_PKL" >&2; exit 1; }

STAMP=$(date +%Y%m%d_%H%M%S)
MAP_DIR="work/metaphlan_genome_sizes_${STAMP}"
SEL_DIR="work/yachida_baseline_manifest_${STAMP}"
CANONICAL=work/<yachida_definitive_run>/canonical_input.tsv   # set explicitly

# 1. Inspect the real schema before trusting the extractor's defaults.
#    Confirm keys terminate in t__SGB... and that lengths are whole-genome.
python3 analysis_v2/scripts/extract_metaphlan_genome_sizes.py \
  --database "$DB_PKL" --outdir /tmp/mpa_inspect --inspect-only

# 2. Extract the authoritative mapping (exact output path, no wildcard).
python3 analysis_v2/scripts/extract_metaphlan_genome_sizes.py \
  --database "$DB_PKL" --outdir "$MAP_DIR"
test -s "$MAP_DIR/SUCCESS" || { echo "extraction failed"; exit 1; }
MAPPING="$MAP_DIR/metaphlan_genome_sizes.tsv"

# 3. Select sample-wide zero-dose MetaPhlAn baselines by column name.
python3 analysis_v2/scripts/select_baseline_manifest.py \
  --canonical "$CANONICAL" \
  --cohort yachida \
  --profiler metaphlan4 \
  --expected-profiles 201 \
  --outdir "$SEL_DIR"

# Verify the count independently; never rely on documentation saying it
# "should" be 201.
test -s "$SEL_DIR/SUCCESS" || { echo "selection failed"; exit 1; }
awk -F '\t' 'NR > 1 {n++} END {exit !(n == 201)}' \
  "$SEL_DIR/baseline_manifest.tsv" || { echo "not 201 baselines"; exit 1; }
MANIFEST="$SEL_DIR/baseline_manifest.tsv"

# Confirm terminal SGB rows really exist in a raw production profile before
# trusting the rank policy. Expect a non-zero count.
grep -c 't__SGB' "$(awk -F '\t' 'NR == 2 {print $8}' "$MANIFEST")"

# 4. Primary audit: 95% coverage, zero tolerated exclusions.
PRIMARY="work/yachida_geff_dev_${STAMP}"
python3 analysis_v2/scripts/compute_effective_genome_size.py \
  --manifest "$MANIFEST" --genome-sizes "$MAPPING" \
  --outdir "$PRIMARY" --profile-rank sgb --min-coverage 0.95
echo DEVELOPMENT_ONLY > "$PRIMARY/DEVELOPMENT_ONLY.txt"

# 5. Prespecified 90% sensitivity, separate directory. The raised exclusion
#    limit only observes the distribution; the primary run must exclude none.
SENS="work/yachida_geff_cov90_dev_${STAMP}"
python3 analysis_v2/scripts/compute_effective_genome_size.py \
  --manifest "$MANIFEST" --genome-sizes "$MAPPING" \
  --outdir "$SENS" --profile-rank sgb --min-coverage 0.90 \
  --max-excluded-samples 201
echo DEVELOPMENT_ONLY > "$SENS/DEVELOPMENT_ONLY.txt"
```

Expected sealed Yachida state: 201 verified samples, 201 retained-output
receipts, 3,408 MetaPhlAn entries, state directory as recorded in `CLAUDE.md`.
Step 3 should select 201 baseline profiles.

Before freezing, inspect the coverage distribution in
`effective_genome_size_audit.tsv` and the ranked gaps in
`unmapped_features.tsv`. If real coverage makes 95% inappropriate, record the
evidence and choose a defensible alternative — but never select a threshold by
comparing recovery performance. The policy is **not** `DECIDED` until these
commands succeed against the real database and sealed baselines.
