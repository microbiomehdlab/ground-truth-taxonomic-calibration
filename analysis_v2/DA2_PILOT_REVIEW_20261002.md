# DA2 paired pilot and next validation — 2026-10-02

User-reported cluster evidence: job 3097925, all 12 contexts passed in
21–29 seconds; all context checksums verified. Source commit 73a2ebe;
results `work/da2_paired_difference_pilot_20261002T002553Z`.

All ten corresponding target taxa appear among positive community discoveries
in all six cohort/profiler contexts. Extra discoveries require separate taxon
identity and compositional interpretation; do not call them false positives.
Individual Fnuc: five of six contexts q<=0.05. Yachida–Bracken has positive
effect, p=6.88581805282181e-5, q=0.22606140667414 over 3,283 features.
This is failed discovery at the chosen threshold, not absent recovery.
Large transformed effects involving zeros are pseudocount-dependent.
These are artificial within-person increases, not clinical validation.

## Conditional paired-null diagnostic

Run 1,000 seeded person-level sign reversals per context, a single sign shared
by all species for each person. Use the saved differences, original numerical
tolerances and unchanged full BH family. Save every sign vector and every draw.
Report two-sided any-discovery rates and conditional Monte Carlo binomial CIs;
also positive discovery counts. Do not substitute randomized p values for
primary p values or select method settings using favorable results.

This diagnostic assumes joint sign symmetry/exchangeability of differences.
Spiked differences are not empirical biological-null samples. The conditional
distribution retains their magnitudes and feature dependence but cannot prove
t-test validity under skewness, establish biological FDR, or replace DA3's
baseline-only randomized-label null. n=10 has only 1,024 possible sign vectors;
sampling is with replacement, with duplicates retained and unique counts shown.
Intervals quantify Monte Carlo uncertainty, not generalization to new studies.
No automatic production PASS threshold is imposed by this diagnostic.

Remaining gates: assess diagnostics, add independent simulated-null scenarios
(including sparse/asymmetric stress cases) with enough replicates, evaluate
pseudocount sensitivity and DA3 baseline-null calibration, then review and
freeze production policy. A problematic sign-flip result requires investigation,
not outcome-selected switching to a more significant test.
