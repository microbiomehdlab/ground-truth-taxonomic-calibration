# Bracken all-input recovery: additive comparison

## What changes and what does not

This standalone module produces target recovery tables beside the existing
native-fraction comparison. It does not edit upstream profiles, seals, canonical
tables, main analysis policy or DA inputs. It is not yet wired into the full
downstream runner. MetaPhlAn endpoints are unchanged.

The denominator audit completed on the cluster as job3097806, source commit
a1580c8, output work/bracken_denominator_audit_20260930T145109Z. User-reported
counts: Yachida3408/Feng3032/Zeller3048 profiles. A subsequent read-only terminal
check matched all8977 spike profiles to R+N and their baselines to R, with no
unused design rows. This implementation repeats that reconciliation rather than
trusting the terminal report. Do not overwrite that audit's NOT_YET_VERIFIED flag.

## Equations

Let R be baseline input pairs, N the total inserted pairs, n_t target inserted
pairs, B_t baseline Bracken new_est_reads and O_t spiked new_est_reads.

- Baseline all-input abundance: B_t/R.
- Spiked all-input abundance: O_t/(R+N).
- Retained baseline contribution: B_t/(R+N).
- Expected abundance under ideal read-proportional recovery: (B_t+n_t)/(R+N).
- Recovered added signal: (O_t-B_t)/(R+N).
- Recovery ratio: (O_t-B_t)/n_t.
- Signed error: [O_t-(B_t+n_t)]/(R+N).

Negative signals are retained. These are estimated classified pair counts,
not directly observed per-target origin assignments or absolute cell abundances.
The expectation assumes perfect recovery of the addition and baseline retention;
changes in reassignment can affect it and are part of measured response.

The legacy comparison is reproduced from the original printed native fractions
and canonical achieved dose fields, including their existing eight-decimal
rounding. Exact count-based fractions are separate columns. Thus old results
are not silently relabelled as all-input values. Neither fraction is renormalized.
Count-positive but printed-fraction-zero targets are flagged separately; the
native detection definition is not silently changed by this quantitation update.

## Community evidence limitation

Retained community design TSVs record total inserted pairs, not each member's
observed count. The existing canonical builder reconstructs member counts from
the frozen ORDERED panel weights and the deterministic integer allocator.
This module uses that same allocator and cross-checks the resulting counts
against canonical target counts. It never assumes an exact total/10 share.

This is provenance-based reconstruction, not independent per-member read-origin
verification. Output says frozen_panel_integer_reconstruction explicitly; the
required CLI option acknowledges it. Historical execution must correspond to
the frozen panel/order/weights. No claim that the comparison independently
proves historical member identity or weights should be made. Input panel and
alias hashes are checked against taxon_identity_freeze.sha256.

## Validation and outputs

Verify audit checksums, consumed Kraken/Bracken/design hashes, exact profile
coverage, all design-row coverage, every target canonical row, aliases/native
values, total/mate counts and achieved fraction rounding (half a unit at eight
decimal places). Missing/corrupt inputs fail before any output directory is
created. Unknown targets, duplicate keys and overlapping/existing output paths
fail. Inputs are rehashed again before writing outputs.

- bracken_recovery_comparison.tsv: exact all-input and legacy endpoints,
  counts/fractions/identities and allocation basis.
- pair_reconciliation.tsv: one row per spiked profile.
- input_hashes.tsv: consumed input and relevant code provenance.
- summary.json: PASS_ENDPOINT_CONSTRUCTION, not proof of good recovery.
- SHA256SUMS: new output integrity manifest.

Expected endpoint counts with the full original design: Yachida15870,
Feng12580, Zeller12720 (41170 total). Community profiles contribute ten targets;
individual profiles contribute only their implanted target. Baseline rows support
the calculation but are not positive-dose recovery endpoints.

## Run after publication of these code files to the branch

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
python3 -B analysis_v2/tests/test_bracken_denominators.py
python3 -B analysis_v2/tests/test_bracken_all_input_recovery.py
export PROJECT="$PWD"
export AUDIT_ROOT="$PWD/work/bracken_denominator_audit_20260930T145109Z"
export PREFLIGHT_STAMP=20260929T155304Z
export RECOVERY_ROOT="$PWD/work/bracken_all_input_recovery_$(date -u +%Y%m%dT%H%M%SZ)"
RECOVERY_JOB="$(sbatch --parsable --export=ALL analysis_v2/run_bracken_all_input_recovery.sbatch)"
printf 'Job: %s\nOutput: %s\n' "$RECOVERY_JOB" "$RECOVERY_ROOT"
```

Only submit if tests pass. After completion inspect accounting and stdout/stderr;
verify all three output checksum files and summaries before interpreting values.
No source commit/push or cluster run is implied by creating this document.
Use the common CLI with explicit canonical paths if the preflight timestamp differs.

## Next scientific step

Inspect paired old/new distributions by cohort, target, population, dose and
condition, including nondetections and negative values. Quantify denominator
effects before ranking profilers or changing manuscript conclusions. Native
Bracken values remain the planned DA input. Any all-input DA sensitivity must
be separately specified; this module does not authorize changing the primary DA.
