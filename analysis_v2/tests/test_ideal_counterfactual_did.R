#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)),
                                "../.."))
root <- tempfile("ideal_counterfactual_did.")
dir.create(root, recursive = TRUE)
manifest_path <- file.path(root, "manifest.tsv")
abundance_path <- file.path(root, "abundance.tsv")
endpoints_path <- file.path(root, "endpoints.tsv")
outdir <- file.path(root, "model")

samples <- paste0("S", seq_len(30))
condition <- rep(c("Control", "Adenoma", "CRC"), each = 10)
within_group <- rep(seq_len(10), 3)
age <- 44 + within_group
sex <- rep(rep(c("Female", "Male"), 5), 3)
doses <- data.frame(
  dose_level = c("baseline", "dose_01", "dose_02"),
  F = c(0, 0.10, 0.20), f = c(0, 0.01, 0.02),
  stringsAsFactors = FALSE
)
manifest <- merge(data.frame(
  sample_id = samples, condition = condition, age = age, sex = sex,
  stringsAsFactors = FALSE
), doses, by = NULL)
manifest <- manifest[order(manifest$dose_level, manifest$sample_id), ]
manifest$cohort <- "synthetic"
manifest$study <- "Study1"
manifest$analysis_population <- "community"
manifest$target_label <- "Fnuc"
manifest$assembly_arm <- "original"
manifest$profiler <- "metaphlan4"
manifest$profile_id <- paste(manifest$sample_id, manifest$dose_level,
                             sep = "_")
manifest$baseline_profile_id <- paste0(manifest$sample_id, "_baseline")
manifest$spike_fraction_target <- manifest$f
manifest$source_profile <- file.path(root, paste0(manifest$profile_id, ".tsv"))
manifest$target_feature <- "Fusobacterium nucleatum"
manifest$bmi <- 23
manifest$include <- 1
manifest$exclusion_reason <- ""
manifest$F <- NULL
manifest$f <- NULL
write.table(manifest, manifest_path, sep = "\t", quote = FALSE,
            row.names = FALSE)

features <- c("Fusobacterium nucleatum", "Stable bystander",
              "Distorted bystander", "Rare")
abundance_rows <- list()
counter <- 1L
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  sample_index <- match(row$sample_id, samples)
  dose_index <- match(row$dose_level, doses$dose_level)
  F <- doses$F[dose_index]
  f <- doses$f[dose_index]
  baseline <- c(
    0.010 + sample_index * 1e-5,
    0.020 + ifelse(row$condition == "CRC", 0.004, 0),
    0.015 + sample_index * 2e-5,
    ifelse(sample_index == 1, 0.001, 0)
  )
  if (row$dose_level == "baseline") {
    observed <- baseline
  } else {
    expected <- (1 - F) * baseline
    expected[1] <- expected[1] + f
    noise <- ((sample_index - 1L) %% 5L - 2L) * 0.01
    residual <- c(0, 0,
                  ifelse(row$condition == "CRC", 0.70, 0) + noise, 0)
    observed <- (expected + 1e-8) * 2^residual - 1e-8
  }
  abundance_rows[[counter]] <- data.frame(
    profiler = row$profiler, source_profile = row$source_profile,
    feature = features, abundance_fraction = observed,
    stringsAsFactors = FALSE
  )
  counter <- counter + 1L
}
abundance <- do.call(rbind, abundance_rows)
write.table(abundance, abundance_path, sep = "\t", quote = FALSE,
            row.names = FALSE)

positive <- manifest[manifest$dose_level != "baseline", ]
endpoints <- positive[c(
  "cohort", "study", "sample_id", "condition", "analysis_population",
  "target_label", "assembly_arm", "profiler", "profile_id",
  "baseline_profile_id", "source_profile"
)]
dose_index <- match(positive$dose_level, doses$dose_level)
endpoints$spike_fraction_total <- doses$F[dose_index]
endpoints$spike_fraction_target <- doses$f[dose_index]
endpoints$source_baseline_profile <- file.path(
  root, paste0(positive$sample_id, "_baseline.tsv")
)
write.table(endpoints, endpoints_path, sep = "\t", quote = FALSE,
            row.names = FALSE)

script <- file.path(repo,
                    "analysis_v2/scripts/fit_ideal_counterfactual_did.R")
result <- system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path, "--outdir", outdir,
  "--reference-scale", "read_proportional",
  "--write-sample-residuals"
), stdout = TRUE, stderr = TRUE)
status <- attr(result, "status")
if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))

stopifnot(file.exists(file.path(outdir, "SUCCESS")))
sample_residuals <- read.delim(file.path(
  outdir, "ideal_counterfactual_sample_residuals.tsv"
))
models <- read.delim(file.path(
  outdir, "ideal_counterfactual_did_results.tsv"
))
stopifnot(nrow(sample_residuals) == 180L, nrow(models) == 6L)
stopifnot(setequal(unique(models$feature_role),
                   c("implanted_target", "bystander")))
target_row <- sample_residuals[
  sample_residuals$sample_id == "S1" &
    sample_residuals$dose_level == "dose_01" &
    sample_residuals$feature == "Fusobacterium nucleatum", ][1, ]
expected <- (1 - target_row$spike_fraction_total) *
  target_row$baseline_abundance_fraction + target_row$spike_fraction_target
stopifnot(abs(target_row$ideal_expected_abundance_fraction - expected) < 1e-12)
stable <- models[models$feature == "Stable bystander", ]
distorted <- models[models$feature == "Distorted bystander", ]
target <- models[models$feature == "Fusobacterium nucleatum", ]
stopifnot(all(abs(stable$effect) < 1e-10), all(stable$p_value == 1),
          all(abs(target$effect) < 1e-10), all(target$p_value == 1),
          all(distorted$effect > 0.65), all(distorted$q_value < 0.05),
          !any(models$feature == "Rare"),
          all(models$contrast ==
                "CRC_vs_Control_in_counterfactual_residual"))

# Full runs do not materialize the very large sample-by-feature table unless
# explicitly requested.
compact_outdir <- file.path(root, "compact_model")
compact <- system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path, "--outdir", compact_outdir,
  "--reference-scale", "read_proportional",
  "--require-complete-panel"
), stdout = TRUE, stderr = TRUE)
compact_status <- attr(compact, "status")
if (is.null(compact_status)) compact_status <- 0L
stopifnot(compact_status == 0L,
          !file.exists(file.path(compact_outdir,
                                 "ideal_counterfactual_sample_residuals.tsv")),
          file.exists(file.path(compact_outdir, "SUCCESS")))

# The fraction contract is fail-closed.
bad_endpoints <- endpoints
bad_endpoints$spike_fraction_total[1] <-
  bad_endpoints$spike_fraction_target[1] / 2
bad_path <- file.path(root, "bad_endpoints.tsv")
write.table(bad_endpoints, bad_path, sep = "\t", quote = FALSE,
            row.names = FALSE)
bad <- suppressWarnings(system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", bad_path,
  "--reference-scale", "read_proportional",
  "--outdir", file.path(root, "bad_model")
), stdout = TRUE, stderr = TRUE))
bad_status <- attr(bad, "status")
if (is.null(bad_status)) bad_status <- 0L
stopifnot(bad_status != 0L)

# A development intersection is permitted, but a definitive complete-panel
# request rejects the same missing context.
incomplete_manifest <- manifest[-which(
  manifest$sample_id == "S1" & manifest$dose_level == "dose_02"
)[1], ]
incomplete_path <- file.path(root, "incomplete_manifest.tsv")
write.table(incomplete_manifest, incomplete_path, sep = "\t", quote = FALSE,
            row.names = FALSE)
incomplete <- suppressWarnings(system2("Rscript", c(
  script, "--profile-manifest", incomplete_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path,
  "--outdir", file.path(root, "incomplete_model"),
  "--require-complete-panel"
), stdout = TRUE, stderr = TRUE))
incomplete_status <- attr(incomplete, "status")
if (is.null(incomplete_status)) incomplete_status <- 0L
stopifnot(incomplete_status != 0L)

cat("[PASS] ideal-counterfactual difference-in-differences fixture\n")

# No silent default: omitting --reference-scale must fail.
no_scale <- suppressWarnings(system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path,
  "--outdir", file.path(root, "no_scale")
), stdout = TRUE, stderr = TRUE))
no_scale_status <- attr(no_scale, "status")
stopifnot(!is.null(no_scale_status), no_scale_status != 0L,
          any(grepl("--reference-scale", no_scale, fixed = TRUE)))
cat("[PASS] reference-scale gate\n")

# --- profiler_scale path -----------------------------------------------------
# The expected abundance is taken verbatim from the response table, so the
# residual must equal the manual log2 ratio computed from that same column.
resp_path <- file.path(root, "responses.tsv")
resp_rows <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  if (row$dose_level == "baseline") next
  sample_index <- match(row$sample_id, samples)
  for (j in seq_along(features)) {
    # A deliberately distinct expected value: not (1-F)o, so a regression to the
    # read-proportional formula would change the residual and fail below.
    expected_ps <- 0.001 * j + 0.0001 * sample_index
    observed_ps <- expected_ps * 1.25
    resp_rows[[length(resp_rows) + 1L]] <- data.frame(
      cohort = "synthetic", study = row$study,
      analysis_population = "community", sample_id = row$sample_id,
      condition = row$condition, profiler = "metaphlan4",
      assembly_arm = "original", dose_level = row$dose_level,
      feature = features[j],
      baseline_abundance_fraction = 0.01,
      observed_abundance_fraction = observed_ps,
      expected_abundance_profiler_scale = expected_ps,
      reference_type = "genome_equivalent",
      stringsAsFactors = FALSE)
  }
}
resp <- do.call(rbind, resp_rows)
write.table(resp, resp_path, sep = "\t", quote = FALSE, row.names = FALSE)

ps_outdir <- file.path(root, "model_profiler_scale")
ps <- suppressWarnings(system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path,
  "--reference-scale", "profiler_scale",
  "--response-table", resp_path,
  "--outdir", ps_outdir, "--write-sample-residuals"
), stdout = TRUE, stderr = TRUE))
ps_status <- attr(ps, "status")
if (is.null(ps_status)) ps_status <- 0L
if (ps_status != 0L) stop(paste(ps, collapse = "\n"))

ps_residuals <- read.delim(
  file.path(ps_outdir, "ideal_counterfactual_sample_residuals.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(ps_residuals) > 0)

pseudocount <- 1e-8
stopifnot(all(c("sample_id", "dose_level", "feature",
                "counterfactual_log2_residual") %in% names(ps_residuals)))
merged <- merge(
  ps_residuals,
  resp[, c("sample_id", "dose_level", "feature",
           "expected_abundance_profiler_scale")],
  by = c("sample_id", "dose_level", "feature"), all.x = TRUE)
# The merge must actually match; an empty result would make the check vacuous.
stopifnot(nrow(merged) == nrow(ps_residuals), nrow(merged) > 0,
          !anyNA(merged$expected_abundance_profiler_scale))
manual <- log2(merged$observed_abundance_fraction + pseudocount) -
          log2(merged$expected_abundance_profiler_scale + pseudocount)
stopifnot(length(manual) == nrow(merged), all(is.finite(manual)))
stopifnot(max(abs(merged$counterfactual_log2_residual - manual)) < 1e-9)
cat("[PASS] profiler-scale DiD residual matches manual calculation\n")

# profiler_scale without --response-table must fail.
no_table <- suppressWarnings(system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path,
  "--reference-scale", "profiler_scale",
  "--outdir", file.path(root, "no_table")
), stdout = TRUE, stderr = TRUE))
no_table_status <- attr(no_table, "status")
stopifnot(!is.null(no_table_status), no_table_status != 0L,
          any(grepl("--response-table", no_table, fixed = TRUE)))

# A wrong within-profiler reference type must fail.
bad_resp <- resp; bad_resp$reference_type <- "read_proportional"
bad_resp_path <- file.path(root, "bad_responses.tsv")
write.table(bad_resp, bad_resp_path, sep = "\t", quote = FALSE, row.names = FALSE)
bad_ref <- suppressWarnings(system2("Rscript", c(
  script, "--profile-manifest", manifest_path,
  "--abundance-long", abundance_path,
  "--paired-endpoints", endpoints_path,
  "--reference-scale", "profiler_scale",
  "--response-table", bad_resp_path,
  "--outdir", file.path(root, "bad_ref")
), stdout = TRUE, stderr = TRUE))
bad_ref_status <- attr(bad_ref, "status")
stopifnot(!is.null(bad_ref_status), bad_ref_status != 0L,
          any(grepl("requires genome_equivalent", bad_ref, fixed = TRUE)))
cat("[PASS] profiler-scale DiD validation gates\n")
