# HC3 comparison and resumable DA3 batches

No automatic authorization. Do not launch the whole DA3 grid solely because
the engineering tests pass. Partial-null checks and the clinical comparison
retain their documented assumptions and uncertainty.

DA1 comparison keeps identical OLS coefficients, age/sex design, transformation
and frozen family, and replaces standard errors with HC3 leverage-adjusted
sandwich estimates. P values use a residual-df t approximation, NOT exact
finite-sample inference. Source: https://sandwich.r-forge.r-project.org/reference/vcovHC.html
Local tests compare an independent covariance-matrix formula and lm
coefficients; sandwich package cross-check is conditional on installation.
HC3 is a candidate comparison, not a claim that heteroskedasticity or sparse
microbiome inference is solved. Actual MaAsLin2 outputs remain unchanged.

Submit only `submit_robust_clinical_checks.sh` now: 12 actual native comparisons
and 36 paired-method synthetic checks on the actual designs. Same synthetic
draws are used for ordinary and HC3 inference; report both 200-repetition
Monte Carlo intervals. This does not rerun the completed DA3 validation.
The script uses committed source snapshots; parallel code work is safe.

DA3 worker now locks each batch, writes into a separate attempt, publishes
only a complete checksummed result, and verifies plan/code/source identity
before skipping an already completed batch. Failed attempts are retained;
do not manually promote them. A changed source or policy rejects reuse.
Older canary outputs without resume identities are not migrated or reused.

Collector streams compact targets and summaries, checks every planned context
and ten unique targets per context, handles exact/MC schemas explicitly and
reports missing/corrupt batches. No incomplete run is labelled complete.
`run_da3_batches.sbatch` supports an offset for bounded arrays: split 4,800
batches into six 800-task arrays to avoid high task indices, limiting combined
concurrency explicitly. Final submission/resource/retry runbook remains to
be reviewed before launch. No production jobs are submitted by these changes.
