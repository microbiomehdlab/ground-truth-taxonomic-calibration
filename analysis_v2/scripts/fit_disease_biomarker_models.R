#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args); if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
manifest_path <- value("--profile-manifest")
abundance_path <- value("--abundance-long")
outdir <- value("--outdir")
min_prevalence <- as.numeric(value("--min-prevalence", "0.10"))
pseudocount <- as.numeric(value("--pseudocount", "1e-8"))
if (is.null(manifest_path) || is.null(abundance_path) || is.null(outdir))
  stop("Required: --profile-manifest FILE --abundance-long FILE --outdir DIR")
if (!is.finite(min_prevalence) || min_prevalence < 0 || min_prevalence > 1 ||
    !is.finite(pseudocount) || pseudocount <= 0) stop("Invalid filter or pseudocount.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
abundance <- read.delim(abundance_path, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("cohort", "study", "analysis_population", "sample_id", "condition",
            "target_label", "assembly_arm", "profiler", "dose_level",
            "spike_fraction_target", "source_profile", "target_feature",
            "age", "sex", "bmi", "include")
if (length(setdiff(needed, names(manifest)))) stop("Profile manifest columns missing.")
if (length(setdiff(c("profiler", "source_profile", "feature", "abundance_fraction"), names(abundance))))
  stop("Abundance columns missing.")
manifest <- manifest[manifest$include == 1, , drop = FALSE]
manifest$dose <- as.numeric(manifest$spike_fraction_target)
manifest$age_numeric <- suppressWarnings(as.numeric(manifest$age))
manifest$bmi_numeric <- suppressWarnings(as.numeric(manifest$bmi))
manifest$sex <- trimws(manifest$sex)
abundance$value <- as.numeric(abundance$abundance_fraction)
if (!nrow(manifest) || anyNA(manifest$dose) || any(manifest$dose < 0) ||
    anyNA(abundance$value) || any(abundance$value < 0)) stop("Invalid numeric input.")
if (anyNA(manifest$age_numeric) || any(!manifest$sex %in% c("Female", "Male")))
  stop("Primary age/sex covariates must be complete and coded Female/Male.")
manifest$source_profile <- normalizePath(manifest$source_profile, mustWork = FALSE)
abundance$source_profile <- normalizePath(abundance$source_profile, mustWork = FALSE)

family_columns <- c("cohort", "study", "analysis_population", "target_label",
                    "assembly_arm", "profiler")
families <- split(manifest, interaction(manifest[family_columns], drop = TRUE, lex.order = TRUE))
if (!length(families)) stop("No disease-analysis families.")

fit_hc3 <- function(y, metadata, include_bmi = FALSE) {
  used <- character(); omitted <- character()
  if (length(unique(metadata$age_numeric)) > 1L) {
    metadata$age_scaled <- as.numeric(scale(metadata$age_numeric)); used <- c(used, "age")
  } else omitted <- c(omitted, "age_invariant")
  if (length(unique(metadata$sex)) > 1L) used <- c(used, "sex")
  else omitted <- c(omitted, "sex_invariant")
  if (include_bmi) {
    if (length(unique(metadata$bmi_numeric)) > 1L) {
      metadata$bmi_scaled <- as.numeric(scale(metadata$bmi_numeric)); used <- c(used, "bmi")
    } else omitted <- c(omitted, "bmi_invariant")
  }
  formula_terms <- c("condition", if ("age" %in% used) "age_scaled",
                     if ("sex" %in% used) "sex", if ("bmi" %in% used) "bmi_scaled")
  metadata$condition <- factor(metadata$condition, levels = c("Control", "Adenoma", "CRC"))
  if ("sex" %in% used) metadata$sex <- factor(metadata$sex)
  formula <- reformulate(formula_terms)
  X <- model.matrix(formula, metadata)
  if (sd(y) < .Machine$double.eps) {
    null <- c(effect = 0, standard_error = 0, lower_95 = 0, upper_95 = 0,
              p_value = 1, df = nrow(metadata) - ncol(X))
    return(list(CRC_vs_Control = null, Adenoma_vs_Control = null,
                covariates_used = paste(used, collapse = ","),
                covariates_omitted = paste(omitted, collapse = ",")))
  }
  fit <- lm.fit(X, y)
  if (fit$rank != ncol(X) || fit$df.residual < 2L) return(NULL)
  inverse <- solve(crossprod(X))
  hat <- rowSums((X %*% inverse) * X)
  adjusted <- fit$residuals / pmax(1 - hat, 1e-8)
  weighted <- sweep(X, 1L, adjusted, `*`)
  covariance <- inverse %*% crossprod(weighted) %*% inverse
  extract <- function(name) {
    index <- match(name, colnames(X))
    se <- unname(sqrt(covariance[index, index]))
    estimate <- unname(fit$coefficients[index])
    statistic <- estimate / se
    c(effect = estimate, standard_error = se,
      lower_95 = estimate + qt(.025, fit$df.residual) * se,
      upper_95 = estimate + qt(.975, fit$df.residual) * se,
      p_value = pmax(2 * pt(abs(statistic), fit$df.residual, lower.tail = FALSE),
                     .Machine$double.xmin), df = fit$df.residual)
  }
  list(CRC_vs_Control = extract("conditionCRC"),
       Adenoma_vs_Control = extract("conditionAdenoma"),
       covariates_used = paste(used, collapse = ","),
       covariates_omitted = paste(omitted, collapse = ","))
}

make_matrix <- function(records, species) {
  result <- matrix(0, nrow = nrow(records), ncol = length(species),
                   dimnames = list(records$sample_id, species))
  for (i in seq_len(nrow(records))) {
    hits <- abundance[abundance$profiler == records$profiler[i] &
                        abundance$source_profile == records$source_profile[i] &
                        abundance$feature %in% species, , drop = FALSE]
    result[records$sample_id[i], hits$feature] <- hits$value
  }
  result
}

key_string <- function(row, columns) paste(vapply(columns, function(x) as.character(row[[x]][1]), character(1)), collapse = "\r")
global_columns <- c("cohort", "study", "analysis_population", "profiler")

fit_collection <- function(matrix, metadata, context) {
  primary <- list(); sensitivity <- list(); bmi_complete <- is.finite(metadata$bmi_numeric)
  for (feature in colnames(matrix)) {
    y <- log2(matrix[, feature] + pseudocount)
    primary_fit <- fit_hc3(y, metadata, FALSE)
    if (is.null(primary_fit)) stop("Rank-deficient primary disease model: ", context)
    primary[[feature]] <- primary_fit
    if (sum(bmi_complete) >= 12L && all(c("Control", "Adenoma", "CRC") %in% unique(metadata$condition[bmi_complete]))) {
      sensitivity_fit <- fit_hc3(y[bmi_complete], metadata[bmi_complete, ], TRUE)
      if (is.null(sensitivity_fit)) stop("Rank-deficient BMI sensitivity model: ", context)
      sensitivity[[feature]] <- sensitivity_fit
    }
  }
  list(primary = primary, sensitivity = sensitivity, bmi_complete = bmi_complete)
}

# Freeze one baseline-derived universe and one set of baseline fits per cohort/profiler.
global_groups <- split(manifest, interaction(manifest[global_columns], drop = TRUE, lex.order = TRUE))
global_configs <- list()
for (global in global_groups) {
  first <- global[1, , drop = FALSE]
  candidates <- global[global$dose == 0, , drop = FALSE]
  identity_columns <- c("sample_id", "condition", "source_profile", "age_numeric", "sex", "bmi_numeric")
  per_sample <- split(candidates[identity_columns], candidates$sample_id)
  if (any(vapply(per_sample, function(x) nrow(unique(x)) != 1L, logical(1))))
    stop("Inconsistent duplicated baseline metadata or profile paths.")
  baseline <- candidates[!duplicated(candidates$sample_id), , drop = FALSE]
  baseline <- baseline[order(baseline$sample_id), , drop = FALSE]
  if (!all(c("Control", "Adenoma", "CRC") %in% unique(baseline$condition)))
    stop("Control, Adenoma, and CRC samples are required in every cohort/profiler baseline.")
  observed_species <- sort(unique(abundance$feature[abundance$profiler == first$profiler &
                                                     abundance$source_profile %in% baseline$source_profile]))
  observed_matrix <- make_matrix(baseline, observed_species)
  prevalence <- colMeans(observed_matrix > 0)
  implanted_targets <- sort(unique(global$target_feature))
  universe <- sort(union(observed_species[prevalence >= min_prevalence], implanted_targets))
  baseline_matrix <- make_matrix(baseline, universe)[baseline$sample_id, , drop = FALSE]
  context <- paste(unlist(first[1, global_columns, drop = FALSE]), collapse = "/")
  fits <- fit_collection(baseline_matrix, baseline, paste0(context, "/baseline"))
  global_configs[[key_string(first, global_columns)]] <- list(
    metadata = baseline, universe = universe, baseline_fits = fits,
    excluded_baseline_species = setdiff(observed_species, universe))
}

render_fits <- function(fits, model_spec, target_list, common, level, dose, metadata, bmi_complete) {
  for (contrast in c("CRC_vs_Control", "Adenoma_vs_Control")) {
    valid <- names(fits)
    pvalues <- vapply(valid, function(feature) fits[[feature]][[contrast]][["p_value"]], numeric(1))
    qvalues <- p.adjust(pvalues, "BH")
    for (feature in valid) {
      stats <- fits[[feature]][[contrast]]
      target_list[[length(target_list) + 1L]] <- c(common, list(
        dose_level = level, spike_fraction_target = dose, contrast = contrast,
        feature = feature, effect = stats[["effect"]], standard_error = stats[["standard_error"]],
        lower_95 = stats[["lower_95"]], upper_95 = stats[["upper_95"]],
        p_value = stats[["p_value"]], q_value = qvalues[feature], model_spec = model_spec,
        n_samples = if (model_spec == "primary_age_sex") nrow(metadata) else sum(bmi_complete),
        covariates_used = fits[[feature]]$covariates_used,
        covariates_omitted = fits[[feature]]$covariates_omitted,
        include = 1, exclusion_reason = "", baseline_reference_kind = "observed_calls"))
    }
  }
  target_list
}

primary_rows <- list(); sensitivity_rows <- list(); exclusions <- list(); ei <- 1L
for (family in families) {
  first <- family[1, , drop = FALSE]
  config <- global_configs[[key_string(first, global_columns)]]
  metadata <- config$metadata; universe <- config$universe
  levels <- c("baseline", sort(unique(family$dose_level[family$dose > 0])))
  level_records <- lapply(levels, function(level) family[family$dose_level == level, , drop = FALSE])
  sample_sets <- lapply(level_records, function(x) sort(x$sample_id))
  if (any(vapply(level_records, function(x) anyDuplicated(x$sample_id) > 0, logical(1))) ||
      !all(vapply(sample_sets, identical, logical(1), sort(metadata$sample_id))))
    stop("Dose levels do not contain the frozen unique baseline samples.")
  common <- as.list(first[1, family_columns, drop = FALSE])
  family_species <- sort(unique(abundance$feature[abundance$profiler == first$profiler &
                                                   abundance$source_profile %in% family$source_profile]))
  excluded_species <- setdiff(family_species, universe)
  for (level_index in seq_along(levels)) {
    records <- level_records[[level_index]][match(metadata$sample_id, level_records[[level_index]]$sample_id), , drop = FALSE]
    dose <- if (levels[level_index] == "baseline") 0 else median(records$dose)
    for (feature in excluded_species) {
      exclusions[[ei]] <- c(common, list(dose_level = levels[level_index], spike_fraction_target = dose,
        contrast = "ALL", feature = feature, model_spec = "not_tested", include = 0,
        exclusion_reason = "outside_frozen_baseline_universe")); ei <- ei + 1L
    }
    fits <- if (levels[level_index] == "baseline") config$baseline_fits else {
      matrix <- make_matrix(records, universe)[metadata$sample_id, , drop = FALSE]
      fit_collection(matrix, metadata, paste(unlist(common), levels[level_index], sep = "/"))
    }
    primary_rows <- render_fits(fits$primary, "primary_age_sex", primary_rows, common,
                                levels[level_index], dose, metadata, fits$bmi_complete)
    if (length(fits$sensitivity))
      sensitivity_rows <- render_fits(fits$sensitivity, "sensitivity_age_sex_bmi_complete_case",
                                      sensitivity_rows, common, levels[level_index], dose,
                                      metadata, fits$bmi_complete)
  }
}

bind <- function(x) if (length(x)) do.call(rbind.data.frame, c(x, stringsAsFactors = FALSE)) else data.frame()
primary <- bind(primary_rows); sensitivity <- bind(sensitivity_rows); excluded <- bind(exclusions)
numeric_fields <- c("spike_fraction_target", "effect", "standard_error", "lower_95",
                    "upper_95", "p_value", "q_value", "n_samples")
for (field in intersect(numeric_fields, names(primary))) primary[[field]] <- as.numeric(primary[[field]])
for (field in intersect(numeric_fields, names(sensitivity))) sensitivity[[field]] <- as.numeric(sensitivity[[field]])
write.table(primary, file.path(outdir, "primary_disease_da_results.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(sensitivity, file.path(outdir, "sensitivity_bmi_disease_da_results.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(excluded, file.path(outdir, "disease_da_exclusions.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
settings <- data.frame(setting = c("primary_formula", "sensitivity_formula", "transformation",
                                   "pseudocount_fraction", "minimum_prevalence", "standard_errors",
                                   "feature_universe", "multiplicity_family", "primary_contrast", "secondary_contrast"),
                       value = c("condition + varying(scaled_age, sex)",
                                 "condition + varying(scaled_age, sex, scaled_bmi)",
                                 "log2(abundance_fraction + fixed pseudocount)", pseudocount,
                                 min_prevalence, "HC3",
                                 "baseline prevalence within cohort/profiler plus all implanted targets",
                                 "BH within cohort/profiler/target/arm/dose/contrast/model over an identical universe",
                                 "CRC_vs_Control", "Adenoma_vs_Control"))
write.table(settings, file.path(outdir, "disease_da_settings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(metric = c("families", "primary_tests", "bmi_sensitivity_tests",
                                 "excluded_rows", "status"),
                      value = c(length(families), nrow(primary), nrow(sensitivity), nrow(excluded), "PASS"))
write.table(summary, file.path(outdir, "disease_da_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
files <- c(normalizePath(manifest_path), normalizePath(abundance_path), list.files(outdir, full.names = TRUE))
status <- system2("sha256sum", files, stdout = file.path(outdir, "disease_da.sha256"))
if (!identical(status, 0L)) stop("Could not seal disease model.")
writeLines(c("analysis\tnative_disease_biomarker_models", "primary_covariates\tage,sex",
             "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Disease biomarker models completed: ", outdir)
