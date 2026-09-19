#!/usr/bin/env Rscript

# Fit phenotype-associated departures from an ideal read-proportional spike
# counterfactual. This is deliberately separate from call-set propagation:
# loss of statistical significance is not evidence that a baseline call was a
# false positive.

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}

# --- quantitative reference scale --------------------------------------------
# Bracken keeps the read-proportional reference, which this script reconstructs
# as (1-F)o with +f on the focal target. MetaPhlAn's primary scale is
# genome-equivalent and is NOT reconstructed here: under profiler_scale the
# expected abundance is read verbatim from the migrated all-feature response
# table, which already carries the complete-mixture Q_i. The scale is always
# stated explicitly; there is no default.
reference_scale <- value("--reference-scale", NULL)
if (is.null(reference_scale))
  stop("--reference-scale {read_proportional|profiler_scale} is required; ",
       "there is no default reference scale.")
if (!reference_scale %in% c("read_proportional", "profiler_scale"))
  stop("Unknown --reference-scale: ", reference_scale)
profiler_scale <- identical(reference_scale, "profiler_scale")
# profiler_scale takes the expected abundance from the already-migrated
# all-feature response table, which carries the correctly reconstructed
# complete-mixture Q_i. Q_i is never recalculated here.
response_table <- value("--response-table", NULL)
if (profiler_scale && is.null(response_table))
  stop("--reference-scale profiler_scale requires --response-table ",
       "(paired_feature_responses.parquet or a lossless TSV export).")
if (profiler_scale && !file.exists(response_table))
  stop("Response table not found: ", response_table)

profiler_scale_expected <- NULL
if (profiler_scale) {
  if (grepl("\\.parquet$", response_table, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("Reading a Parquet response table requires the 'arrow' package; ",
           "export a lossless TSV instead.")
    resp <- as.data.frame(arrow::read_parquet(response_table))
  } else {
    resp <- read.delim(response_table, check.names = FALSE,
                       stringsAsFactors = FALSE)
  }
  needed <- c("cohort", "study", "analysis_population", "sample_id", "condition",
              "profiler", "assembly_arm", "dose_level", "feature",
              "baseline_abundance_fraction", "observed_abundance_fraction",
              "expected_abundance_profiler_scale", "reference_type")
  absent <- setdiff(needed, names(resp))
  if (length(absent))
    stop("Response table lacks column(s): ", paste(absent, collapse = ", "))
  # Validate the profiler/reference mapping within profiler.
  expected_reference <- c(kraken2_bracken = "read_proportional",
                          metaphlan4 = "genome_equivalent")
  for (prof in unique(resp$profiler)) {
    if (!prof %in% names(expected_reference))
      stop("Unknown profiler in response table: ", prof)
    types <- unique(resp$reference_type[resp$profiler == prof])
    if (length(types) != 1L)
      stop("Profiler ", prof, " carries multiple reference types: ",
           paste(types, collapse = ", "))
    if (!identical(types, unname(expected_reference[[prof]])))
      stop("Profiler ", prof, " has reference_type ", types,
           "; the primary scale requires ", expected_reference[[prof]])
  }
  profiler_scale_expected <- resp
}

manifest_path <- value("--profile-manifest")
abundance_path <- value("--abundance-long")
endpoints_path <- value("--paired-endpoints")
outdir <- value("--outdir")
min_prevalence <- as.numeric(value("--min-prevalence", "0.10"))
pseudocount <- as.numeric(value("--pseudocount", "1e-8"))
write_sample_residuals <- "--write-sample-residuals" %in% args
require_complete_panel <- "--require-complete-panel" %in% args

if (is.null(manifest_path) || is.null(abundance_path) ||
    is.null(endpoints_path) || is.null(outdir)) {
  stop(paste(
    "Required: --profile-manifest FILE --abundance-long FILE",
    "--paired-endpoints FILE --outdir DIR"
  ))
}
if (!is.finite(min_prevalence) || min_prevalence < 0 || min_prevalence > 1 ||
    !is.finite(pseudocount) || pseudocount <= 0) {
  stop("Invalid prevalence threshold or pseudocount.")
}
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE,
                                            no.. = TRUE))) {
  stop("Output directory exists and is not empty: ", outdir)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_path, check.names = FALSE,
                       stringsAsFactors = FALSE)
abundance <- read.delim(abundance_path, check.names = FALSE,
                        stringsAsFactors = FALSE)
endpoints <- read.delim(endpoints_path, check.names = FALSE,
                        stringsAsFactors = FALSE)

manifest_required <- c(
  "cohort", "study", "analysis_population", "sample_id", "condition",
  "target_label", "assembly_arm", "profiler", "profile_id",
  "baseline_profile_id", "spike_fraction_target", "dose_level",
  "source_profile", "target_feature", "age", "sex", "bmi", "include",
  "exclusion_reason"
)
abundance_required <- c(
  "profiler", "source_profile", "feature", "abundance_fraction"
)
endpoint_key <- c(
  "cohort", "study", "sample_id", "condition", "analysis_population",
  "target_label", "assembly_arm", "profiler", "profile_id",
  "baseline_profile_id"
)
endpoint_required <- c(
  endpoint_key, "spike_fraction_total", "spike_fraction_target",
  "source_baseline_profile", "source_profile"
)
if (length(setdiff(manifest_required, names(manifest))))
  stop("Profile manifest columns missing: ",
       paste(setdiff(manifest_required, names(manifest)), collapse = ", "))
if (length(setdiff(abundance_required, names(abundance))))
  stop("Abundance columns missing: ",
       paste(setdiff(abundance_required, names(abundance)), collapse = ", "))
if (length(setdiff(endpoint_required, names(endpoints))))
  stop("Paired-endpoint columns missing: ",
       paste(setdiff(endpoint_required, names(endpoints)), collapse = ", "))

manifest$include <- as.character(manifest$include)
manifest$exclusion_reason <- trimws(as.character(manifest$exclusion_reason))
manifest$exclusion_reason[is.na(manifest$exclusion_reason)] <- ""
if (any(!manifest$include %in% c("0", "1")) ||
    any((manifest$include == "0") != nzchar(manifest$exclusion_reason))) {
  stop("Invalid include/exclusion encoding in profile manifest.")
}
manifest <- manifest[manifest$include == "1", , drop = FALSE]
if (!nrow(manifest)) stop("No included profile contexts.")

manifest$spike_fraction_target <- suppressWarnings(
  as.numeric(manifest$spike_fraction_target)
)
manifest$age_numeric <- suppressWarnings(as.numeric(manifest$age))
manifest$sex <- trimws(as.character(manifest$sex))
abundance$abundance_fraction <- suppressWarnings(
  as.numeric(abundance$abundance_fraction)
)
endpoints$spike_fraction_total <- suppressWarnings(
  as.numeric(endpoints$spike_fraction_total)
)
endpoints$spike_fraction_target <- suppressWarnings(
  as.numeric(endpoints$spike_fraction_target)
)
if (anyNA(manifest$spike_fraction_target) ||
    any(manifest$spike_fraction_target < 0) || anyNA(manifest$age_numeric) ||
    any(!manifest$sex %in% c("Female", "Male"))) {
  stop("Invalid manifest fractions or primary age/sex covariates.")
}
if (anyNA(abundance$abundance_fraction) ||
    any(!is.finite(abundance$abundance_fraction)) ||
    any(abundance$abundance_fraction < 0) ||
    any(abundance$abundance_fraction > 1.00001)) {
  stop("Abundance values must be finite fractions on [0,1].")
}
if (anyNA(endpoints$spike_fraction_total) ||
    anyNA(endpoints$spike_fraction_target) ||
    any(!is.finite(endpoints$spike_fraction_total)) ||
    any(!is.finite(endpoints$spike_fraction_target))) {
  stop("Invalid endpoint fractions.")
}

manifest$source_profile <- normalizePath(manifest$source_profile,
                                         mustWork = FALSE)
abundance$source_profile <- normalizePath(abundance$source_profile,
                                          mustWork = FALSE)
endpoints$source_profile <- normalizePath(endpoints$source_profile,
                                         mustWork = FALSE)
endpoints$source_baseline_profile <- normalizePath(
  endpoints$source_baseline_profile, mustWork = FALSE
)

# Index the long abundance table once. Repeatedly scanning the complete table
# for every sample and dose is prohibitively slow in the full three-cohort run.
profile_token <- function(profiler, source_profile) {
  paste(profiler, source_profile, sep = "\r")
}
abundance_by_profile <- split(
  seq_len(nrow(abundance)),
  profile_token(abundance$profiler, abundance$source_profile)
)

key_string <- function(data, columns) {
  do.call(paste, c(data[columns], sep = "\r"))
}
manifest_identity <- c(endpoint_key, "dose_level")
if (anyDuplicated(key_string(manifest, manifest_identity)))
  stop("Duplicate included profile context.")
if (anyDuplicated(key_string(abundance, c("profiler", "source_profile",
                                          "feature"))))
  stop("Duplicate abundance feature within a native profile.")
if (anyDuplicated(key_string(endpoints, endpoint_key)))
  stop("Duplicate paired endpoint context.")

positive <- manifest$spike_fraction_target > 0
if (any(manifest$dose_level == "baseline" & positive) ||
    any(manifest$dose_level != "baseline" & !positive)) {
  stop("dose_level and spike_fraction_target disagree.")
}
manifest_endpoint_key <- key_string(manifest[positive, , drop = FALSE],
                                    endpoint_key)
endpoint_keys <- key_string(endpoints, endpoint_key)
endpoint_index <- match(manifest_endpoint_key, endpoint_keys)
if (anyNA(endpoint_index))
  stop("A positive profile context lacks an exact paired endpoint.")
matched_endpoints <- endpoints[endpoint_index, , drop = FALSE]
if (any(abs(manifest$spike_fraction_target[positive] -
            matched_endpoints$spike_fraction_target) > 1e-12)) {
  stop("Manifest and endpoint target fractions disagree.")
}
F <- matched_endpoints$spike_fraction_total
f <- matched_endpoints$spike_fraction_target
if (any(f <= 0) || any(F < f - 1e-12) || any(F >= 1))
  stop("Positive endpoints require 0 < f <= F < 1.")
if (any(manifest$source_profile[positive] !=
        matched_endpoints$source_profile)) {
  stop("Manifest and endpoint perturbed profile paths disagree.")
}
manifest$spike_fraction_total <- 0
manifest$endpoint_fraction_target <- 0
manifest$endpoint_source_baseline <- ""
manifest$spike_fraction_total[positive] <- F
manifest$endpoint_fraction_target[positive] <- f
manifest$endpoint_source_baseline[positive] <-
  matched_endpoints$source_baseline_profile

# Match the disease-model policy for incomplete development inputs: freeze one
# complete biological-sample intersection per cohort/study/population. Fully
# complete definitive inputs retain every eligible sample.
panel_columns <- c("target_label", "assembly_arm", "profiler", "dose_level")
panel_group_columns <- c("cohort", "study", "analysis_population")
panel_groups <- split(
  manifest,
  interaction(manifest[panel_group_columns], drop = TRUE, lex.order = TRUE)
)
retained_panel_groups <- list()
sample_panel_audit <- list()
for (panel_index in seq_along(panel_groups)) {
  panel_group <- panel_groups[[panel_index]]
  first_panel <- panel_group[1, , drop = FALSE]
  expected_cells <- nrow(unique(panel_group[panel_columns]))
  observed <- unique(panel_group[c("sample_id", panel_columns)])
  observed_counts <- table(observed$sample_id)
  retained_samples <- names(observed_counts)[observed_counts == expected_cells]
  if (length(retained_samples) < 2L)
    stop("Common ideal-counterfactual sample intersection has fewer than two samples.")
  retained_panel_groups[[panel_index]] <- panel_group[
    panel_group$sample_id %in% retained_samples, , drop = FALSE
  ]
  sample_panel_audit[[panel_index]] <- data.frame(
    cohort = first_panel$cohort, study = first_panel$study,
    analysis_population = first_panel$analysis_population,
    input_samples = length(unique(panel_group$sample_id)),
    retained_samples = length(retained_samples),
    excluded_samples = length(unique(panel_group$sample_id)) -
      length(retained_samples),
    expected_cells_per_sample = expected_cells,
    policy = "complete_common_target_arm_profiler_dose_intersection",
    stringsAsFactors = FALSE
  )
}
manifest <- do.call(rbind, retained_panel_groups)
sample_panel_audit <- do.call(rbind, sample_panel_audit)
if (require_complete_panel && any(sample_panel_audit$excluded_samples != 0L)) {
  stop("A definitive ideal-reference analysis requires a complete sample panel.")
}

make_matrix <- function(records, species) {
  result <- matrix(0, nrow = nrow(records), ncol = length(species),
                   dimnames = list(records$sample_id, species))
  for (i in seq_len(nrow(records))) {
    indices <- abundance_by_profile[[profile_token(
      records$profiler[i], records$source_profile[i]
    )]]
    if (is.null(indices)) next
    hits <- abundance[indices, , drop = FALSE]
    hits <- hits[hits$feature %in% species, , drop = FALSE]
    if (nrow(hits)) {
      result[i, match(hits$feature, species)] <- hits$abundance_fraction
    }
  }
  result
}

fit_crc_hc3 <- function(y, metadata) {
  used <- character()
  omitted <- character()
  metadata$condition <- factor(metadata$condition,
                               levels = c("Control", "Adenoma", "CRC"))
  if (length(unique(metadata$age_numeric)) > 1L) {
    metadata$age_scaled <- as.numeric(scale(metadata$age_numeric))
    used <- c(used, "age")
  } else {
    omitted <- c(omitted, "age_invariant")
  }
  if (length(unique(metadata$sex)) > 1L) {
    metadata$sex <- factor(metadata$sex)
    used <- c(used, "sex")
  } else {
    omitted <- c(omitted, "sex_invariant")
  }
  terms <- c("condition", if ("age" %in% used) "age_scaled",
             if ("sex" %in% used) "sex")
  X <- model.matrix(reformulate(terms), metadata)
  if (!"conditionCRC" %in% colnames(X))
    stop("CRC versus Control is not estimable.")
  if (sd(y) < 1e-12) {
    return(c(
      effect = 0, standard_error = 0, lower_95 = 0, upper_95 = 0,
      p_value = 1, df = nrow(metadata) - ncol(X),
      covariates_used = paste(used, collapse = ","),
      covariates_omitted = paste(omitted, collapse = ",")
    ))
  }
  fit <- lm.fit(X, y)
  if (fit$rank != ncol(X) || fit$df.residual < 2L)
    stop("Rank-deficient ideal-counterfactual model.")
  inverse <- solve(crossprod(X))
  hat <- rowSums((X %*% inverse) * X)
  adjusted <- fit$residuals / pmax(1 - hat, 1e-8)
  weighted <- sweep(X, 1L, adjusted, `*`)
  covariance <- inverse %*% crossprod(weighted) %*% inverse
  index <- match("conditionCRC", colnames(X))
  se <- sqrt(max(covariance[index, index], 0))
  estimate <- unname(fit$coefficients[index])
  if (!is.finite(se) || se <= .Machine$double.eps) {
    p <- if (abs(estimate) < 1e-12) 1 else .Machine$double.xmin
  } else {
    p <- pmax(2 * pt(abs(estimate / se), fit$df.residual,
                     lower.tail = FALSE), .Machine$double.xmin)
  }
  c(
    effect = estimate, standard_error = se,
    lower_95 = estimate + qt(.025, fit$df.residual) * se,
    upper_95 = estimate + qt(.975, fit$df.residual) * se,
    p_value = p, df = fit$df.residual,
    covariates_used = paste(used, collapse = ","),
    covariates_omitted = paste(omitted, collapse = ",")
  )
}

global_columns <- c("cohort", "study", "analysis_population", "profiler")
family_columns <- c(global_columns, "target_label", "assembly_arm")
global_groups <- split(
  manifest,
  interaction(manifest[global_columns], drop = TRUE, lex.order = TRUE)
)
global_configs <- list()
for (global in global_groups) {
  first <- global[1, , drop = FALSE]
  baseline_candidates <- global[global$dose_level == "baseline", , drop = FALSE]
  identity <- c("sample_id", "condition", "source_profile", "age_numeric", "sex")
  per_sample <- split(baseline_candidates[identity],
                      baseline_candidates$sample_id)
  if (any(vapply(per_sample, function(x) nrow(unique(x)) != 1L, logical(1))))
    stop("Inconsistent duplicated baseline profiles or metadata.")
  baseline <- baseline_candidates[!duplicated(baseline_candidates$sample_id),
                                  , drop = FALSE]
  baseline <- baseline[order(baseline$sample_id), , drop = FALSE]
  if (!all(c("Control", "Adenoma", "CRC") %in% baseline$condition))
    stop("Control, Adenoma, and CRC are required in every analysis panel.")
  observed_species <- sort(unique(abundance$feature[
    abundance$profiler == first$profiler &
      abundance$source_profile %in% baseline$source_profile
  ]))
  if (!length(observed_species)) stop("Baseline profiles contain no species.")
  observed <- make_matrix(baseline, observed_species)
  prevalence <- colMeans(observed > 0)
  targets <- sort(unique(global$target_feature))
  universe <- sort(union(observed_species[prevalence >= min_prevalence],
                         targets))
  baseline_matrix <- make_matrix(baseline, universe)
  global_configs[[key_string(first, global_columns)]] <- list(
    metadata = baseline, universe = universe,
    baseline_matrix = baseline_matrix,
    excluded_species = setdiff(observed_species, universe)
  )
}

sample_path <- file.path(outdir, "ideal_counterfactual_sample_residuals.tsv")
result_path <- file.path(outdir, "ideal_counterfactual_did_results.tsv")
sample_written <- FALSE
result_written <- FALSE
sample_rows_total <- 0L
result_rows_total <- 0L
families <- split(
  manifest,
  interaction(manifest[family_columns], drop = TRUE, lex.order = TRUE)
)
panel_audit <- list()

for (family_index in seq_along(families)) {
  family <- families[[family_index]]
  first <- family[1, , drop = FALSE]
  config <- global_configs[[key_string(first, global_columns)]]
  metadata <- config$metadata
  universe <- config$universe
  target_feature <- unique(family$target_feature)
  if (length(target_feature) != 1L || !target_feature %in% universe)
    stop("A family must have exactly one target feature in the frozen universe.")
  baseline_family <- family[family$dose_level == "baseline", , drop = FALSE]
  if (anyDuplicated(baseline_family$sample_id) ||
      !identical(sort(baseline_family$sample_id), sort(metadata$sample_id)))
    stop("A family lacks the frozen unique baseline sample panel.")
  baseline_family <- baseline_family[match(metadata$sample_id,
                                           baseline_family$sample_id), , drop = FALSE]
  positive_levels <- sort(unique(family$dose_level[family$dose_level != "baseline"]))
  if (!length(positive_levels)) stop("A family has no positive dose.")
  panel_audit[[family_index]] <- data.frame(
    first[1, family_columns, drop = FALSE],
    samples = nrow(metadata), positive_doses = length(positive_levels),
    frozen_features = length(universe), target_feature = target_feature,
    stringsAsFactors = FALSE
  )
  for (level in positive_levels) {
    records <- family[family$dose_level == level, , drop = FALSE]
    if (anyDuplicated(records$sample_id) ||
        !identical(sort(records$sample_id), sort(metadata$sample_id)))
      stop("A positive dose lacks the frozen unique baseline sample panel.")
    records <- records[match(metadata$sample_id, records$sample_id), , drop = FALSE]
    if (any(records$endpoint_source_baseline != baseline_family$source_profile))
      stop("Endpoint baseline path does not match the family baseline profile.")
    F <- records$spike_fraction_total
    f <- records$endpoint_fraction_target
    observed <- make_matrix(records, universe)
    if (profiler_scale) {
      # Authoritative expected abundance, already renormalised with the
      # complete-mixture Q_i. Nothing is reconstructed from F and the focal f.
      lookup <- profiler_scale_expected
      key_resp <- paste(lookup$cohort, lookup$analysis_population,
                        lookup$profiler, lookup$assembly_arm,
                        lookup$dose_level, lookup$sample_id, lookup$feature,
                        sep = "\r")
      n_rows <- nrow(records); k_cols <- length(universe)
      want <- paste(rep(first$cohort, n_rows * k_cols),
                    rep(first$analysis_population, n_rows * k_cols),
                    rep(first$profiler, n_rows * k_cols),
                    rep(first$assembly_arm, n_rows * k_cols),
                    rep(level, n_rows * k_cols),
                    rep(records$sample_id, times = k_cols),
                    rep(universe, each = n_rows), sep = "\r")
      hit <- match(want, key_resp)
      if (anyNA(hit))
        stop("Response table lacks ", sum(is.na(hit)),
             " feature/sample cell(s) required by this family.")
      expected <- matrix(
        as.numeric(lookup$expected_abundance_profiler_scale[hit]),
        nrow = n_rows, ncol = k_cols,
        dimnames = dimnames(config$baseline_matrix))
      target_index <- match(target_feature, universe)
    } else {
      expected <- sweep(config$baseline_matrix, 1L, 1 - F, `*`)
      target_index <- match(target_feature, universe)
      expected[, target_index] <- expected[, target_index] + f
    }
    if (any(!is.finite(expected)) || any(expected < 0) ||
        any(expected > 1.00001))
      stop("Ideal counterfactual abundance is outside the fraction scale.")
    residual <- log2(observed + pseudocount) - log2(expected + pseudocount)
    if (any(!is.finite(residual))) stop("Non-finite counterfactual residual.")

    n <- nrow(records)
    k <- length(universe)
    sample_rows_total <- sample_rows_total + n * k
    if (write_sample_residuals) {
      sample_rows <- data.frame(
        cohort = rep(first$cohort, n * k),
        study = rep(first$study, n * k),
        analysis_population = rep(first$analysis_population, n * k),
        sample_id = rep(records$sample_id, times = k),
        condition = rep(records$condition, times = k),
        target_label = rep(first$target_label, n * k),
        assembly_arm = rep(first$assembly_arm, n * k),
        profiler = rep(first$profiler, n * k),
        dose_level = rep(level, n * k),
        spike_fraction_total = rep(F, times = k),
        spike_fraction_target = rep(f, times = k),
        target_feature = rep(target_feature, n * k),
        feature = rep(universe, each = n),
        feature_role = rep(ifelse(universe == target_feature,
                                  "implanted_target", "bystander"), each = n),
        baseline_abundance_fraction = as.vector(config$baseline_matrix),
        observed_abundance_fraction = as.vector(observed),
        ideal_expected_abundance_fraction = as.vector(expected),
        counterfactual_log2_residual = as.vector(residual),
        age = rep(records$age_numeric, times = k),
        sex = rep(records$sex, times = k),
        stringsAsFactors = FALSE
      )
      write.table(sample_rows, sample_path, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = !sample_written,
                  append = sample_written)
      sample_written <- TRUE
    }

    fits <- lapply(seq_along(universe), function(j) {
      fit_crc_hc3(residual[, j], records)
    })
    p_values <- vapply(fits, function(x) as.numeric(x[["p_value"]]),
                       numeric(1))
    q_values <- p.adjust(p_values, method = "BH")
    model_rows <- do.call(rbind, lapply(seq_along(universe), function(j) {
      fit <- fits[[j]]
      data.frame(
        first[1, family_columns, drop = FALSE], dose_level = level,
        median_spike_fraction_total = median(F),
        median_spike_fraction_target = median(f),
        contrast = "CRC_vs_Control_in_counterfactual_residual",
        feature = universe[j],
        feature_role = ifelse(universe[j] == target_feature,
                              "implanted_target", "bystander"),
        effect = as.numeric(fit[["effect"]]),
        standard_error = as.numeric(fit[["standard_error"]]),
        lower_95 = as.numeric(fit[["lower_95"]]),
        upper_95 = as.numeric(fit[["upper_95"]]),
        p_value = p_values[j], q_value = q_values[j],
        model_spec = "ideal_counterfactual_residual_age_sex_hc3",
        n_samples = n,
        covariates_used = as.character(fit[["covariates_used"]]),
        covariates_omitted = as.character(fit[["covariates_omitted"]]),
        stringsAsFactors = FALSE
      )
    }))
    write.table(model_rows, result_path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = !result_written,
                append = result_written)
    result_written <- TRUE
    result_rows_total <- result_rows_total + nrow(model_rows)
  }
}

panel_audit <- do.call(rbind, panel_audit)
write.table(panel_audit,
            file.path(outdir, "ideal_counterfactual_panel_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_panel_audit,
            file.path(outdir, "ideal_counterfactual_sample_panel_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
settings <- data.frame(
  setting = c(
    "estimand", "expected_bystander", "expected_target", "residual",
    "primary_formula", "pseudocount_fraction", "minimum_prevalence",
    "standard_errors", "multiplicity_family", "complete_panel_required",
    "interpretation"
  ),
  value = c(
    "CRC-minus-Control difference in counterfactual log2 residual",
    "e=(1-F)*o", "e=(1-F)*o+f", "u=log2(a+1e-8)-log2(e+1e-8)",
    "u ~ condition + varying(scaled_age, sex)", pseudocount,
    min_prevalence, "HC3",
    "BH within cohort/study/population/target/arm/profiler/dose",
    require_complete_panel,
    "phenotype-dependent departure from an ideal read-proportional reference; not false-positive adjudication"
  ), stringsAsFactors = FALSE
)
write.table(settings, file.path(outdir, "ideal_counterfactual_settings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  metric = c("families", "sample_residual_rows_available",
             "sample_residual_rows_written", "model_rows", "status"),
  value = c(length(families), sample_rows_total,
            ifelse(write_sample_residuals, sample_rows_total, 0),
            result_rows_total, "PASS"),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(outdir, "ideal_counterfactual_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(),
               file = file.path(outdir, "ideal_counterfactual_session_info.txt"))

source_files <- normalizePath(c(manifest_path, abundance_path, endpoints_path))
output_files <- c(
  if (write_sample_residuals) sample_path else character(), result_path,
  file.path(outdir, "ideal_counterfactual_panel_audit.tsv"),
  file.path(outdir, "ideal_counterfactual_sample_panel_audit.tsv"),
  file.path(outdir, "ideal_counterfactual_settings.tsv"),
  file.path(outdir, "ideal_counterfactual_summary.tsv"),
  file.path(outdir, "ideal_counterfactual_session_info.txt")
)
checksum_path <- file.path(outdir, "ideal_counterfactual.sha256")
status <- system2("sha256sum", c(source_files, output_files),
                  stdout = checksum_path)
if (!identical(status, 0L)) stop("Could not seal ideal-counterfactual model.")
writeLines(c(
  "analysis\tideal_counterfactual_difference_in_differences",
  paste0("sample_residual_rows_available\t", sample_rows_total),
  paste0("sample_residual_rows_written\t",
         ifelse(write_sample_residuals, sample_rows_total, 0)),
  paste0("model_rows\t", result_rows_total), "status\tPASS"
), file.path(outdir, "SUCCESS"))
message("[PASS] Ideal-counterfactual difference-in-differences completed: ",
        outdir)
