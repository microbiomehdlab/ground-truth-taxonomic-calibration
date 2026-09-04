# Assembly-choice sensitivity model

**Status:** prespecified and fixture-tested; final-data execution pending the
canonical original-versus-clean input gate.

This analysis asks whether replacing the original Pana and Pint assemblies with
cleaner alternatives changes profiler response. It is an assembly-choice and
assembly-quality sensitivity analysis, not a causal contamination experiment.
The original arm remains part of the historical experimental record and must
not be overwritten.

The input is the paired-endpoint table for the same 30 Yachida samples, two
targets, six positive doses, two profilers, and both assembly arms. The primary
Gaussian GAM is:

`recovered_spike_signal ~ profiler * assembly_arm * dose_pp + condition +`
`s(sample_id, bs="re") + s(target_label, bs="re")`

where `dose_pp = 100 * spike_fraction_target`. The reported response slopes are
scaled so that one represents read-proportional response. Primary sensitivity
contrasts are clean minus original response slope within each profiler; the
profiler difference-in-differences is secondary. No rows are filtered based on
their observed response, and negative recovered signals remain in the model.

The script records convergence and Hessian checks, cell counts, all four
profiler-by-arm slopes, contrasts with 95% intervals, fitted residuals, session
information, checksums, and a `SUCCESS` marker. A final-data result is not valid
unless the canonical builder and endpoint derivation both pass first.
