# Assembly-choice sensitivity model

**Status:** prespecified and fixture-tested; final-data execution pending the
canonical original-versus-clean input gate.

This analysis asks whether replacing the original Pana and Pint assemblies with
cleaner alternatives changes profiler response. It is an assembly-choice and
assembly-quality sensitivity analysis, not a causal contamination experiment.
The original arm remains part of the historical experimental record and must
not be overwritten.

The input is the paired-endpoint table for the same 30 Yachida samples, two
targets, six positive doses, two profilers, and both assembly arms.

## Primary biological-sample analysis

The biological sample is the inferential unit. For every sample, target,
profiler, and assembly arm, an ordinary least-squares line with an estimated
intercept is fitted across the six positive implanted fractions. The clean-arm
slope minus the original-arm slope is then calculated within each sample. A
slope of one represents read-proportional response.

The four target-by-profiler mean paired differences are reported with
biological-sample bootstrap percentile intervals and two-sided sign-flip tests.
Sign flips are enumerated for at most 16 samples and otherwise use a
deterministic Monte Carlo procedure. BH adjustment is performed across these
four primary tests. Medians, standard deviations, every sample-level slope, and
every paired difference are retained so that heterogeneity remains visible.

Equal-weight target averages and within-sample profiler
difference-in-differences are secondary. The profiler comparison asks whether
assembly choice changes the two profiler outputs differently; it does not imply
that their native abundance definitions are identical.

## Secondary trajectory model

The secondary Gaussian GAM is:

`recovered_spike_signal ~ profiler * assembly_arm * dose_pp * target_label +`
`condition + sample random intercept + sample-target-profiler random intercept +`
`shared random dose slope + clean-arm random dose-slope deviation`

where `dose_pp = 100 * spike_fraction_target`. Its random terms represent
repeated-measures dependence and heterogeneous trajectories. The GAM records
convergence, Hessian checks, fitted residuals, and model-based contrasts. It is
a trajectory diagnostic, not the source of primary assembly inference. This
ordering was frozen because repeated-row models can produce very narrow
model-based intervals when biological replication is modest, even with random
slopes.

No rows are filtered based on observed response, and negative recovered signals
remain in both analyses. Exact zero p-values caused by floating-point underflow
are prohibited in the secondary model.

The primary script records per-sample fits, paired differences, bootstrap
intervals, sign-flip metadata, adjusted values, diagnostics, session
information, checksums, and a `SUCCESS` marker. A final-data result is invalid
unless canonical construction, native-profile auditing, and endpoint derivation
all pass first.
