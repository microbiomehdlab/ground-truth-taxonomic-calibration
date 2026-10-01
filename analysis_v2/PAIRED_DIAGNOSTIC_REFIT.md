# Paired diagnostic refits (not primary analysis)

Uses corrected audit tables and immutable failed pilot inputs. No species is
selected by significance, and no result replaces a primary result. All eligible
contexts select the lexically first neutral feature ID per diagnostic category:
optimizer failure, code-zero convergence message, clean comparison, deterministic
within-pair difference, plus F. nucleatum in every context with a saved model.
Overlapping selections are combined with all selection reasons recorded.
Maximum five species per context (60 across the twelve failed contexts).

Five records per selected species: saved original, raw-response nloptwrap refit,
centred/SD-scaled response with nloptwrap, bobyqa, and Nelder_Mead. Same binary
group/random person intercept, same REML setting as the saved model; verify
saved response and person/group identities against reconstructed native input.
Transformation remains log2(1+a/1e-8). Scaling is across all observations of that
species, not separately by group/pair. Return beta, standard errors and variance
SDs to original response units. Capture all warnings/messages and optimizer
diagnostics; save every successfully constructed model, including warned models.
No optimizer chosen as winner and no threshold declaring warnings harmless.

Paired-difference one-sample t tests are separately labelled diagnostics, not
primary DA2 or replacements. Constant differences produce an explicit error and
NA p, not an invented extreme p. No q values on this selected diagnostic subset.
Agreement of fits alone is not statistical calibration or validity proof.
Errors of individual refits are retained in comparison rows; bundle PASS means
the diagnostic run completed, not that every model converged or was acceptable.

Outputs: selection.tsv, refit_comparisons.tsv, paired_difference_checks.tsv,
models/, sessionInfo.txt, source commit/image checksum, SUCCESS and SHA256SUMS.
Input audit/source checksums reverified before and after; separate fresh output.
Use refit_paired_diagnostics.sbatch with PROJECT, ANALYSIS_SIF, DA_PILOT_RESULTS,
PAIRED_AUDIT_ROOT and PAIRED_REFIT_ROOT. Production remains paused.

Tests: deterministic diagnostic selection/input-order invariance locally. Actual
mixed-model/rescaling fixture requires lmerTest; enforce REQUIRE_REFIT_BACKEND=1
in the pinned cluster image. Local backend absence is reported, not hidden.

Reference: https://lme4.github.io/lme4/reference/convergence.html describes
optimizer comparison as a convergence diagnostic. This limited set is not the
complete allFit suite and does not automatically adjudicate all warnings.
