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
manifest_required <- c("cohort", "study", "analysis_population", "sample_id", "condition",
                       "target_label", "assembly_arm", "profiler", "profile_id",
                       "baseline_profile_id", "spike_fraction_target", "dose_level", "source_profile",
                       "target_feature", "include", "exclusion_reason")
abundance_required <- c("profiler", "source_profile", "feature", "abundance_fraction")
if (length(setdiff(manifest_required, names(manifest)))) stop("Profile manifest columns missing.")
if (length(setdiff(abundance_required, names(abundance)))) stop("Abundance columns missing.")
manifest <- manifest[manifest$include == 1, , drop = FALSE]
manifest$dose <- as.numeric(manifest$spike_fraction_target)
abundance$value <- as.numeric(abundance$abundance_fraction)
if (!nrow(manifest) || anyNA(manifest$dose) || any(manifest$dose < 0) ||
    anyNA(abundance$value) || any(abundance$value < 0)) stop("Invalid numeric inputs.")
manifest$source_profile <- normalizePath(manifest$source_profile, mustWork = FALSE)
abundance$source_profile <- normalizePath(abundance$source_profile, mustWork = FALSE)
if (anyDuplicated(abundance[c("profiler", "source_profile", "feature")]))
  stop("Duplicate native-profile feature rows.")

family_columns <- c("cohort", "study", "analysis_population", "condition",
                    "target_label", "assembly_arm", "profiler")
families <- split(manifest, interaction(manifest[family_columns], drop = TRUE, lex.order = TRUE))
if (!length(families)) stop("No analysis families.")

feature_rows <- list(); exclusion_rows <- list(); paired_rows <- list()
result_index <- exclusion_index <- pair_index <- 1L
safe_test <- function(delta) {
  if (all(abs(delta) < .Machine$double.eps)) return(c(effect = 0, se = 0, p = 1))
  if (sd(delta) < .Machine$double.eps) return(c(effect = mean(delta), se = 0,
                                                p = .Machine$double.xmin))
  test <- t.test(delta, mu = 0)
  c(effect = mean(delta), se = sd(delta) / sqrt(length(delta)), p = test$p.value)
}

for (family in families) {
  first <- family[1L, , drop = FALSE]
  baseline <- family[family$dose_level == "baseline", , drop = FALSE]
  positive_levels <- sort(unique(family$dose_level[family$dose > 0]))
  if (!length(positive_levels) || nrow(baseline) < 2L || anyDuplicated(baseline$sample_id))
    stop("Each family needs unique baselines and positive dose levels.")
  contexts <- lapply(positive_levels, function(level) family[family$dose_level == level, , drop = FALSE])
  if (any(vapply(contexts, function(context)
    nrow(context) != nrow(baseline) || !setequal(context$sample_id, baseline$sample_id), logical(1))))
    stop("Every dose level requires one matched profile per baseline sample.")
  species <- sort(unique(abundance$feature[
    abundance$profiler == first$profiler & abundance$source_profile %in%
      family$source_profile]))
  species <- union(species, first$target_feature)
  matrix_for <- function(records) {
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
  base <- matrix_for(baseline)
  dose_matrices <- lapply(contexts, matrix_for)
  prevalence <- colMeans(do.call(rbind, c(list(base), dose_matrices)) > 0)
  keep <- prevalence >= min_prevalence | species == first$target_feature
  common <- as.list(first[1L, family_columns, drop = FALSE])
  contrast <- paste0("spiked_vs_matched_baseline__background_", first$condition)
  baseline_rows <- list()
  for (feature in species) {
    if (!keep[feature]) {
      exclusion_rows[[exclusion_index]] <- c(common, list(
        spike_fraction_target = 0, contrast = contrast, feature = feature,
        effect = 0, standard_error = NA, p_value = 1, q_value = 1,
        include = 0, exclusion_reason = "prevalence_below_threshold"))
      exclusion_index <- exclusion_index + 1L
    } else {
      baseline_rows[[feature]] <- c(common, list(
        spike_fraction_target = 0, contrast = contrast, feature = feature,
        effect = 0, standard_error = 0, p_value = 1, q_value = 1,
        include = 1, exclusion_reason = ""))
    }
  }
  for (row in baseline_rows) {
    feature_rows[[result_index]] <- row; result_index <- result_index + 1L
  }
  for (level_index in seq_along(contexts)) {
    context <- contexts[[level_index]]
    spiked <- dose_matrices[[level_index]]
    base_aligned <- base[rownames(spiked), , drop = FALSE]
    transformed_delta <- log2(spiked + pseudocount) - log2(base_aligned + pseudocount)
    achieved_dose <- median(context$dose)
    tested <- list()
    for (feature in species) {
      if (!keep[feature]) {
        exclusion_rows[[exclusion_index]] <- c(common, list(
          spike_fraction_target = achieved_dose, contrast = contrast, feature = feature,
          effect = 0, standard_error = NA, p_value = 1, q_value = 1,
          include = 0, exclusion_reason = "prevalence_below_threshold"))
        exclusion_index <- exclusion_index + 1L
        next
      }
      delta <- transformed_delta[, feature]
      stats <- safe_test(delta); tested[[feature]] <- stats
      for (sample in names(delta)) {
        paired_rows[[pair_index]] <- c(common, list(
          dose_level = positive_levels[level_index], spike_fraction_target = achieved_dose,
          contrast = contrast, sample_id = sample, feature = feature,
          paired_log2_change = unname(delta[sample])))
        pair_index <- pair_index + 1L
      }
    }
    pvalues <- vapply(tested, function(x) x[["p"]], numeric(1))
    qvalues <- p.adjust(pvalues, method = "BH")
    for (feature in names(tested)) {
      stats <- tested[[feature]]
      feature_rows[[result_index]] <- c(common, list(
        spike_fraction_target = achieved_dose, contrast = contrast, feature = feature,
        effect = stats[["effect"]], standard_error = stats[["se"]],
        p_value = stats[["p"]], q_value = qvalues[feature], include = 1,
        exclusion_reason = ""))
      result_index <- result_index + 1L
    }
  }
}

included <- do.call(rbind.data.frame, c(feature_rows, stringsAsFactors = FALSE))
excluded <- if (length(exclusion_rows))
  do.call(rbind.data.frame, c(exclusion_rows, stringsAsFactors = FALSE)) else included[0, ]
calls <- rbind(included, excluded)
pairs <- do.call(rbind.data.frame, c(paired_rows, stringsAsFactors = FALSE))
for (name in c("spike_fraction_target", "effect", "standard_error", "p_value", "q_value"))
  calls[[name]] <- as.numeric(calls[[name]])
calls$include <- as.integer(calls$include)
write.table(calls, file.path(outdir, "paired_da_results.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
write.table(pairs, file.path(outdir, "sample_feature_log2_changes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
settings <- data.frame(setting = c("estimand", "transformation", "pseudocount_fraction",
                                   "minimum_prevalence", "target_filter_exception",
                                   "test", "multiplicity_family"),
                       value = c("mean paired log2 abundance change",
                                 "log2(abundance_fraction + fixed pseudocount)", pseudocount,
                                 min_prevalence, "always retain intended target",
                                 "two-sided one-sample t-test of paired changes",
                                 "BH across included species within each context"))
write.table(settings, file.path(outdir, "paired_da_settings.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
summary <- data.frame(metric = c("contexts", "included_feature_tests", "excluded_feature_rows",
                                 "paired_sample_feature_rows", "status"),
                      value = c(sum(vapply(families, function(x) length(unique(x$dose_level[x$dose > 0])), integer(1))),
                                nrow(included), nrow(excluded), nrow(pairs), "PASS"))
write.table(summary, file.path(outdir, "paired_da_summary.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
files <- c(normalizePath(manifest_path), normalizePath(abundance_path),
           list.files(outdir, full.names = TRUE))
status <- system2("sha256sum", files, stdout = file.path(outdir, "paired_da.sha256"))
if (!identical(status, 0L)) stop("Could not create checksum seal.")
writeLines(c("analysis\tpaired_perturbation_feature_models", "inference_unit\tbiological_sample",
             "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Paired biomarker models completed: ", outdir)
