# Bracken denominator audit

Read-only source inspection; only a new output directory is written. Standard
library Python, compatible with Python 3.8. No containers or raw FASTQs needed.

## Scope

Select samples from each checksummed production_seal_v2 manifest and expected
profile counts from its sample_flow. Select profiles only from sample receipts,
require expected sample-root/category placement and verify size/SHA256 for each
Kraken/Bracken file. Quarantine directories are not scanned. Verify all three
categories: baseline, community and independent. Never overwrite an output.

Counts are Kraken observations, not automatically individual reads. The production
workflow uses --paired; independent reconciliation with achieved spike-pair design
counts is NOT implemented here. This audit explicitly marks that fact. It does
not yet compute per-target recovery or alter canonical abundance/DA inputs.

Output profile_denominators.tsv contains total=classified+unclassified, summed
Bracken species estimates, unclassified fraction, represented fraction and native
fraction sum. Classified-minus-estimated counts are separate from unclassified.
Baseline deltas use the same biological sample, not an arbitrary first profile.
Native fractions may sum slightly below one due to printing; no renormalization
is applied. A decrease in unclassified fraction alone does not prove target recovery.

input_hashes.tsv records selected files and receipts. Design TSVs found in the
receipt are inventoried with headers for the next pair/dose reconciliation step.
summary.json explicitly says recovery_validated=false. SHA256SUMS seals audit
outputs, not the entire upstream dataset anew. Only consumed profile/design
files are rehashed, not every unrelated retained output.

## Run after these files have been committed/pushed and pulled onto the cluster

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
python3 -B analysis_v2/tests/test_bracken_denominators.py
export PROJECT="$PWD"
export AUDIT_ROOT="$PWD/work/bracken_denominator_audit_$(date -u +%Y%m%dT%H%M%SZ)"
BRACKEN_AUDIT_JOB="$(sbatch --parsable --export=ALL analysis_v2/run_bracken_denominator_audit.sbatch)"
printf 'Job: %s\nOutput: %s\n' "$BRACKEN_AUDIT_JOB" "$AUDIT_ROOT"
```

After completion, inspect accounting, stdout/stderr, each summary and checksums.
Expected production totals: Yachida3408, Feng3032, Zeller3048 profiles. The script
also checks counts per sample against the seal; aggregate counts alone are not
sufficient. Three per-cohort PASS_DENOMINATOR_AUDIT_ONLY results mean this stage
passed, not that read-pair/dose reconciliation or biological recovery was validated.

If one cohort fails, earlier successful cohort outputs remain; use a fresh output
root for a complete rerun or the explicit Python CLI with a new cohort outdir.
Never delete seals, edit receipts or ignore count mismatches to force a pass.

## Next stage

Inspect inventoried design schemas; map each profile to achieved baseline/spike
pair counts and target fractions. Test that mapping before deriving all-input
per-target abundance and baseline-dilution-adjusted recovery. Preserve native
Bracken abundance and the old recovery view as separately labelled outputs.
