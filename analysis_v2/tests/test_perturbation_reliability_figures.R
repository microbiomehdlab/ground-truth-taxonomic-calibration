#!/usr/bin/env Rscript
script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(script_arg), "../.."))
tmp <- tempfile("reliability_figures."); dir.create(tmp)
input <- file.path(tmp, "scores.tsv"); out <- file.path(tmp, "out")
d <- expand.grid(cohort = c("feng", "zeller"), analysis_population = c("community", "independent"),
                 profiler = c("kraken2_bracken", "metaphlan4"), feature = paste0("taxon", 1:4),
                 stringsAsFactors = FALSE)
d$contrast <- "CRC_vs_Control"; d$baseline_effect <- 1
d$eligible_stress_tests <- 12; d$eligible_implanted_targets <- 2
d$perturbation_persistence <- rep(c(.4, .7, .9, 1), length.out = nrow(d))
d$direction_stability <- rep(c(.75, 1), length.out = nrow(d))
d$effect_fidelity <- rep(c(.45, .7, .85, .95), length.out = nrow(d))
d$artificial_artifact_score <- .1; d$artifact_resistance <- .9; d$score_components <- 4
d$perturbation_reliability_score <- rowMeans(d[c("perturbation_persistence", "direction_stability", "effect_fidelity", "artifact_resistance")])
d$external_cohorts_evaluable <- 1
d$external_replication_status <- rep(c("DIRECTIONALLY_REPLICATED", "NOT_SIGNIFICANT_ELSEWHERE"), length.out = nrow(d))
d$evidence_tier <- "fixture"
write.table(d, input, sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(file.path(root, "analysis_v2/scripts/plot_perturbation_reliability_scores.R"),
  "--scores", input, "--outdir", out))
stopifnot(status == 0, file.exists(file.path(out, "SUCCESS")),
          file.exists(file.path(out, "figures/reliability_landscape.pdf")),
          file.exists(file.path(out, "figures/reliability_by_replication.pdf")),
          file.exists(file.path(out, "reliability_figures.sha256")))
message("[PASS] perturbation-reliability figure fixture")
