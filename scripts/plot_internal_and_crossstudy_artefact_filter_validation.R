#!/usr/bin/env Rscript

# Cross-study validation of spike-derived abundance-artifact filtering
#
# This script learns non-target artifact lists in one cohort from sample-level
# community spike traces, then applies those frozen lists to the ORIGINAL
# significant DA calls in the other cohort. It does not rerun MaAsLin2.
# Therefore, it evaluates cross-cohort transfer of a post hoc exclusion rule.
# A later workflow can reuse the same learned lists to regenerate feature tables
# and rerun MaAsLin2 if desired.

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

as_num <- function(x) suppressWarnings(as.numeric(x))
as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  z %in% c("true", "t", "1", "yes", "y")
}

norm_taxon <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("^s__", "", x)
  x <- gsub("[|;]", " ", x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x)
}

f1_score <- function(tp, fp, fn) {
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  c(precision = precision, recall = recall, f1 = f1)
}

aggregate_rows <- function(df, by_cols, fun) {
  keys <- interaction(df[by_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, keys)
  out <- lapply(pieces, fun)
  do.call(rbind, out)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

trace_file <- args[["trace"]] %||% "RUNS/spike_metrics/community/species_trace_with_condition.csv"
da_file <- args[["significant-file"]] %||% "RUNS/maaslin_spike/maaslin_significant_features_ALLFILTERS.csv"
outdir <- args[["outdir"]] %||% "RUNS/manuscript_figures/crossstudy_artefact_filter_validation"
q_threshold <- as_num(args[["q-threshold"]] %||% 0.10)
main_threshold <- as_num(args[["threshold"]] %||% 0.05)
thresholds <- as_num(strsplit(args[["thresholds"]] %||% "0.025,0.05,0.10,0.20", ",", fixed = TRUE)[[1]])
min_samples <- as.integer(args[["min-samples"]] %||% 20)
min_median_expected <- as_num(args[["min-median-expected"]] %||% 0)
min_fraction_recurrence <- as_num(args[["min-fraction-recurrence"]] %||% 0)
community_label <- args[["community-label"]] %||% "CRCpanel"
filter_mode <- args[["filter-mode"]] %||% "original"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

message("[INFO] Reading trace: ", trace_file)
trace <- safe_read_csv(trace_file)
message("[INFO] Trace dimensions: ", nrow(trace), " x ", ncol(trace))
message("[INFO] Reading DA results: ", da_file)
da <- safe_read_csv(da_file)
message("[INFO] DA dimensions: ", nrow(da), " x ", ncol(da))

required_trace <- c("base_id", "spike_mode", "spike_label", "spike_fraction_total",
                    "tool", "taxon", "target_flag", "expected", "relative_error",
                    "Target_Condition", "Study")
missing_trace <- setdiff(required_trace, names(trace))
if (length(missing_trace)) stop("Trace is missing columns: ", paste(missing_trace, collapse = ", "))

required_da <- c("feature", "qval", "coef", "background_condition", "background_study",
                 "tool", "spike_mode", "spike_label", "spike_fraction", "filter_mode")
missing_da <- setdiff(required_da, names(da))
if (length(missing_da)) stop("DA file is missing columns: ", paste(missing_da, collapse = ", "))

trace$spike_fraction_total <- as_num(trace$spike_fraction_total)
trace$expected <- as_num(trace$expected)
trace$relative_error <- as_num(trace$relative_error)
trace$target_flag <- as_logical_safe(trace$target_flag)
trace$taxon_key <- norm_taxon(trace$taxon)
trace$tool <- as.character(trace$tool)
trace$Study <- as.character(trace$Study)
trace$Target_Condition <- as.character(trace$Target_Condition)

trace <- trace[
  trace$spike_mode == "community" &
    trace$spike_label == community_label &
    !trace$target_flag &
    is.finite(trace$spike_fraction_total),
  , drop = FALSE
]
if (!nrow(trace)) stop("No community non-target rows remained after filtering.")

# The DA table already contains the important derived flags in most workflows,
# but recreate them robustly if needed.
da$qval <- as_num(da$qval)
da$coef <- as_num(da$coef)
da$spike_fraction <- as_num(da$spike_fraction)
da$is_significant2 <- if ("is_significant" %in% names(da)) {
  as_logical_safe(da$is_significant)
} else {
  is.finite(da$qval) & da$qval <= q_threshold
}
da$is_positive2 <- if ("is_positive" %in% names(da)) {
  as_logical_safe(da$is_positive)
} else {
  is.finite(da$coef) & da$coef > 0
}
da$is_target2 <- if ("is_target" %in% names(da)) {
  as_logical_safe(da$is_target)
} else {
  rep(FALSE, nrow(da))
}
feature_source <- if ("feature_norm" %in% names(da)) da$feature_norm else da$feature
da$taxon_key <- norm_taxon(feature_source)

da <- da[
  da$spike_mode == "community" &
    da$spike_label == community_label &
    da$filter_mode == filter_mode &
    da$is_significant2 &
    da$is_positive2 &
    is.finite(da$spike_fraction),
  , drop = FALSE
]
if (!nrow(da)) stop("No positive significant community DA rows remained after filtering.")

study_levels <- sort(intersect(unique(trace$Study), unique(da$background_study)))
if (length(study_levels) < 2L) {
  stop("At least two shared studies are required. Shared studies: ", paste(study_levels, collapse = ", "))
}
# Use the first two shared studies, matching the manuscript cohorts.
study_levels <- study_levels[1:2]
message("[INFO] Studies: ", paste(study_levels, collapse = ", "))

# Compute per-study, per-profiler, per-fraction, per-taxon artifact statistics.
stat_fun <- function(z) {
  rr <- z$relative_error
  data.frame(
    Study = z$Study[1],
    tool = z$tool[1],
    spike_fraction_total = z$spike_fraction_total[1],
    taxon = z$taxon[1],
    taxon_key = z$taxon_key[1],
    n_samples = sum(is.finite(rr)),
    mean_abs_relative_error = safe_mean(abs(rr)),
    sd_relative_error = safe_sd(rr),
    median_expected = safe_median(z$expected),
    stringsAsFactors = FALSE
  )
}
artifact_stats <- aggregate_rows(
  trace,
  c("Study", "tool", "spike_fraction_total", "taxon_key"),
  stat_fun
)
rownames(artifact_stats) <- NULL
artifact_stats$eligible <- artifact_stats$n_samples >= min_samples &
  is.finite(artifact_stats$median_expected) &
  artifact_stats$median_expected >= min_median_expected

write.csv(artifact_stats, file.path(outdir, "artifact_statistics_by_study_tool_fraction_taxon.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# Internal (within-benchmark) filtering
# -----------------------------------------------------------------------------
# Internal artifact statistics are pooled across both cohorts, matching the
# original proof-of-concept filtering analysis. These lists are applied to DA
# contexts from the same benchmark and therefore represent internal, not
# externally transferred, performance.

internal_stat_fun <- function(z) {
  rr <- z$relative_error
  data.frame(
    tool = z$tool[1],
    spike_fraction_total = z$spike_fraction_total[1],
    taxon = z$taxon[1],
    taxon_key = z$taxon_key[1],
    n_samples = sum(is.finite(rr)),
    mean_abs_relative_error = safe_mean(abs(rr)),
    sd_relative_error = safe_sd(rr),
    median_expected = safe_median(z$expected),
    stringsAsFactors = FALSE
  )
}

internal_stats <- aggregate_rows(
  trace,
  c("tool", "spike_fraction_total", "taxon_key"),
  internal_stat_fun
)
rownames(internal_stats) <- NULL
internal_stats$eligible <- internal_stats$n_samples >= min_samples &
  is.finite(internal_stats$median_expected) &
  internal_stats$median_expected >= min_median_expected
write.csv(internal_stats,
          file.path(outdir, "artifact_statistics_internal_tool_fraction_taxon.csv"),
          row.names = FALSE)

build_internal_artifact_list <- function(threshold) {
  x <- internal_stats[internal_stats$eligible, , drop = FALSE]
  x$artifact <- (is.finite(x$mean_abs_relative_error) & x$mean_abs_relative_error > threshold) |
    (is.finite(x$sd_relative_error) & x$sd_relative_error > threshold)

  if (min_fraction_recurrence > 0) {
    rec_fun <- function(z) {
      data.frame(
        tool = z$tool[1],
        taxon_key = z$taxon_key[1],
        taxon = z$taxon[1],
        n_fractions = nrow(z),
        artifact_fraction = mean(z$artifact),
        artifact = mean(z$artifact) >= min_fraction_recurrence,
        spike_fraction_total = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    return(aggregate_rows(x, c("tool", "taxon_key"), rec_fun))
  }

  x[x$artifact, c("tool", "spike_fraction_total", "taxon", "taxon_key",
                  "n_samples", "mean_abs_relative_error", "sd_relative_error",
                  "median_expected", "artifact"), drop = FALSE]
}

apply_internal_filter <- function(threshold) {
  art <- build_internal_artifact_list(threshold)
  context_cols <- c("tool", "background_study", "background_condition", "spike_fraction")
  context_key <- interaction(da[context_cols], drop = TRUE, lex.order = TRUE)
  contexts <- split(da, context_key)

  out <- lapply(contexts, function(z) {
    tool_i <- z$tool[1]
    frac_i <- z$spike_fraction[1]
    study_i <- z$background_study[1]
    cond_i <- z$background_condition[1]

    if (min_fraction_recurrence > 0) {
      art_keys <- art$taxon_key[art$tool == tool_i & art$artifact]
    } else {
      art_keys <- art$taxon_key[
        art$tool == tool_i &
          abs(art$spike_fraction_total - frac_i) < 1e-12 &
          art$artifact
      ]
    }

    target_keys <- unique(z$taxon_key[z$is_target2])
    offtarget_keys <- unique(z$taxon_key[!z$is_target2])
    retained_offtarget <- setdiff(offtarget_keys, art_keys)

    tp_before <- length(target_keys)
    fp_before <- length(offtarget_keys)
    fn_before <- max(10L - tp_before, 0L)
    m_before <- f1_score(tp_before, fp_before, fn_before)

    tp_after <- tp_before
    fp_after <- length(retained_offtarget)
    fn_after <- fn_before
    m_after <- f1_score(tp_after, fp_after, fn_after)

    universe <- unique(trace$taxon_key[
      trace$tool == tool_i &
        trace$Study == study_i &
        trace$Target_Condition == cond_i &
        abs(trace$spike_fraction_total - frac_i) < 1e-12
    ])
    universe <- universe[nzchar(universe)]

    data.frame(
      threshold = threshold,
      tool = tool_i,
      background_study = study_i,
      background_condition = cond_i,
      spike_fraction = frac_i,
      n_artifact_taxa = length(unique(art_keys)),
      n_test_nontarget_universe = length(universe),
      n_offtarget_before = fp_before,
      n_offtarget_after = fp_after,
      n_offtarget_removed = fp_before - fp_after,
      fraction_test_universe_flagged = if (length(universe)) mean(universe %in% art_keys) else NA_real_,
      fraction_offtarget_DA_flagged = if (length(offtarget_keys)) mean(offtarget_keys %in% art_keys) else NA_real_,
      target_recall_before = unname(m_before["recall"]),
      target_recall_after = unname(m_after["recall"]),
      precision_before = unname(m_before["precision"]),
      precision_after = unname(m_after["precision"]),
      f1_before = unname(m_before["f1"]),
      f1_after = unname(m_after["f1"]),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  attr(result, "artifact_list") <- art
  result
}

internal_results_list <- lapply(thresholds, apply_internal_filter)
internal_results <- do.call(rbind, internal_results_list)
internal_lists <- do.call(rbind, lapply(seq_along(thresholds), function(i) {
  x <- attr(internal_results_list[[i]], "artifact_list")
  x$threshold <- thresholds[i]
  x
}))
rownames(internal_results) <- NULL
rownames(internal_lists) <- NULL
write.csv(internal_results,
          file.path(outdir, "internal_filter_context_metrics.csv"),
          row.names = FALSE)
write.csv(internal_lists,
          file.path(outdir, "internal_training_artefact_lists.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Cross-study transfer filtering
# -----------------------------------------------------------------------------
build_crossstudy_artifact_list <- function(train_study, threshold) {
  x <- artifact_stats[artifact_stats$Study == train_study & artifact_stats$eligible, , drop = FALSE]
  x$artifact <- (is.finite(x$mean_abs_relative_error) & x$mean_abs_relative_error > threshold) |
    (is.finite(x$sd_relative_error) & x$sd_relative_error > threshold)

  if (min_fraction_recurrence > 0) {
    rec_fun <- function(z) {
      data.frame(
        Study = train_study,
        tool = z$tool[1],
        taxon_key = z$taxon_key[1],
        taxon = z$taxon[1],
        n_fractions = nrow(z),
        artifact_fraction = mean(z$artifact),
        artifact = mean(z$artifact) >= min_fraction_recurrence,
        spike_fraction_total = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    return(aggregate_rows(x, c("tool", "taxon_key"), rec_fun))
  }

  x[x$artifact, c("Study", "tool", "spike_fraction_total", "taxon", "taxon_key",
                  "n_samples", "mean_abs_relative_error", "sd_relative_error",
                  "median_expected", "artifact"), drop = FALSE]
}

apply_crossstudy_filter <- function(train_study, test_study, threshold) {
  art <- build_crossstudy_artifact_list(train_study, threshold)
  test_da <- da[da$background_study == test_study, , drop = FALSE]
  if (!nrow(test_da)) return(NULL)

  context_cols <- c("tool", "background_study", "background_condition", "spike_fraction")
  context_key <- interaction(test_da[context_cols], drop = TRUE, lex.order = TRUE)
  contexts <- split(test_da, context_key)

  out <- lapply(contexts, function(z) {
    tool_i <- z$tool[1]
    frac_i <- z$spike_fraction[1]
    cond_i <- z$background_condition[1]

    if (min_fraction_recurrence > 0) {
      art_keys <- art$taxon_key[art$tool == tool_i & art$artifact]
    } else {
      art_keys <- art$taxon_key[
        art$tool == tool_i &
          abs(art$spike_fraction_total - frac_i) < 1e-12 &
          art$artifact
      ]
    }

    target_keys <- unique(z$taxon_key[z$is_target2])
    offtarget_keys <- unique(z$taxon_key[!z$is_target2])
    retained_offtarget <- setdiff(offtarget_keys, art_keys)

    tp_before <- length(target_keys)
    fp_before <- length(offtarget_keys)
    fn_before <- max(10L - tp_before, 0L)
    m_before <- f1_score(tp_before, fp_before, fn_before)

    tp_after <- tp_before
    fp_after <- length(retained_offtarget)
    fn_after <- fn_before
    m_after <- f1_score(tp_after, fp_after, fn_after)

    data.frame(
      train_study = train_study,
      test_study = test_study,
      transfer = paste(train_study, "→", test_study),
      threshold = threshold,
      tool = tool_i,
      background_condition = cond_i,
      spike_fraction = frac_i,
      n_artifact_taxa_train = length(unique(art_keys)),
      n_offtarget_before = fp_before,
      n_offtarget_after = fp_after,
      n_offtarget_removed = fp_before - fp_after,
      target_recall_before = unname(m_before["recall"]),
      target_recall_after = unname(m_after["recall"]),
      precision_before = unname(m_before["precision"]),
      precision_after = unname(m_after["precision"]),
      f1_before = unname(m_before["f1"]),
      f1_after = unname(m_after["f1"]),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  attr(result, "artifact_list") <- art
  result
}

cross_results <- list()
cross_lists <- list()
k <- 1L
for (thr in thresholds) {
  for (i in seq_along(study_levels)) {
    train <- study_levels[i]
    test <- study_levels[setdiff(seq_along(study_levels), i)][1]
    res <- apply_crossstudy_filter(train, test, thr)
    if (!is.null(res)) {
      cross_results[[k]] <- res
      art <- attr(res, "artifact_list")
      art$train_study <- train
      art$test_study <- test
      art$threshold <- thr
      cross_lists[[k]] <- art
      k <- k + 1L
    }
  }
}
crossstudy_results <- do.call(rbind, cross_results)
crossstudy_lists <- do.call(rbind, cross_lists)
rownames(crossstudy_results) <- NULL
rownames(crossstudy_lists) <- NULL
write.csv(crossstudy_results,
          file.path(outdir, "crossstudy_filter_transfer_context_metrics.csv"),
          row.names = FALSE)
write.csv(crossstudy_lists,
          file.path(outdir, "crossstudy_training_artefact_lists.csv"),
          row.names = FALSE)

internal_main <- internal_results[abs(internal_results$threshold - main_threshold) < 1e-12, , drop = FALSE]
cross_main <- crossstudy_results[abs(crossstudy_results$threshold - main_threshold) < 1e-12, , drop = FALSE]
if (!nrow(internal_main) || !nrow(cross_main)) stop("Main threshold was not found among --thresholds.")

# Friendly labels
study_label <- function(x) gsub("_", " ", x)
tool_label_map <- c(
  "kraken2_bracken" = "Kraken2 + Bracken",
  "metaphlan4" = "MetaPhlAn 4"
)
add_labels <- function(x) {
  x$tool_label <- unname(tool_label_map[x$tool])
  x$tool_label[is.na(x$tool_label)] <- x$tool[is.na(x$tool_label)]
  x
}
internal_main <- add_labels(internal_main)
internal_results <- add_labels(internal_results)
cross_main <- add_labels(cross_main)
crossstudy_results <- add_labels(crossstudy_results)
cross_main$transfer_label <- paste(study_label(cross_main$train_study), "→", study_label(cross_main$test_study))
crossstudy_results$transfer_label <- paste(study_label(crossstudy_results$train_study), "→", study_label(crossstudy_results$test_study))

profiler_cols <- c("Kraken2 + Bracken" = "#009E73", "MetaPhlAn 4" = "#6F5BD3")

library(ggplot2)
theme_manuscript <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = "grey55"),
    strip.text = element_text(size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# Panel A: internal mechanistic enrichment.
enrich_long <- rbind(
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             group = "All non-target taxa",
             fraction_flagged = internal_main$fraction_test_universe_flagged),
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             group = "Off-target enriched DA taxa",
             fraction_flagged = internal_main$fraction_offtarget_DA_flagged)
)
enrich_sum <- aggregate(fraction_flagged ~ tool_label + group,
                        data = enrich_long,
                        FUN = function(x) median(x, na.rm = TRUE))

pA <- ggplot(enrich_sum, aes(x = group, y = fraction_flagged,
                             colour = tool_label, group = tool_label)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.75) +
  geom_text(
    aes(label = paste0(round(100 * fraction_flagged, 1), "%")),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = profiler_cols, name = "Profiler") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "A. Off-target DA calls are enriched\namong artefact-prone taxa",
    subtitle = "Internal community benchmark;\nartefacts defined from abundance deviation",
    x = NULL,
    y = "Median fraction flagged"
  ) + theme_manuscript

# Panel B: internal before/after off-target burden by context.
b_long <- rbind(
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             state = "Unfiltered", burden = internal_main$n_offtarget_before),
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             state = "Filtered", burden = internal_main$n_offtarget_after)
)
b_long$context_id <- interaction(b_long$tool_label, b_long$background_study,
                                 b_long$background_condition, b_long$spike_fraction, drop = TRUE)
b_long$state <- factor(b_long$state, levels = c("Unfiltered", "Filtered"))

pB <- ggplot(b_long, aes(x = state, y = burden, group = context_id, colour = tool_label)) +
  geom_line(alpha = 0.40, linewidth = 0.55) +
  geom_point(size = 1.7) +
  facet_wrap(~ tool_label, scales = "free_y") +
  scale_colour_manual(values = profiler_cols, name = "Profiler", guide = "none") +
  labs(
    title = "B. Within-benchmark filtering reduces\noff-target DA burden",
    subtitle = paste0("Pooled artefact lists; threshold = ", 100 * main_threshold, "%;\nprofiler facets use separate y-scales"),
    x = NULL,
    y = "Off-target enriched taxa"
  ) + theme_manuscript

# Panel C: internal precision and F1 before/after.
# Implanted targets are protected from exclusion by design, so target recall is
# not plotted as an empirical filtering outcome.
c_long <- rbind(
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             metric = "Precision", state = "Unfiltered", value = internal_main$precision_before),
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             metric = "Precision", state = "Filtered", value = internal_main$precision_after),
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             metric = "F1", state = "Unfiltered", value = internal_main$f1_before),
  data.frame(internal_main[, c("tool_label", "background_study", "background_condition", "spike_fraction")],
             metric = "F1", state = "Filtered", value = internal_main$f1_after)
)
c_sum <- aggregate(value ~ tool_label + metric + state,
                   data = c_long,
                   FUN = function(x) median(x, na.rm = TRUE))
c_sum$state <- factor(c_sum$state, levels = c("Unfiltered", "Filtered"))
c_sum$metric <- factor(c_sum$metric, levels = c("F1", "Precision"))

pC <- ggplot(c_sum, aes(x = state, y = value, colour = tool_label, group = tool_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  facet_grid(metric ~ tool_label) +
  scale_colour_manual(values = profiler_cols, name = "Profiler", guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  labs(
    title = "C. Internal filtering improves precision and F1",
    subtitle = "Implanted targets and aliases were protected from exclusion by design",
    x = NULL,
    y = "Median performance"
  ) + theme_manuscript

# Panel D: compact cross-study transfer validation using precision and F1.
d_long <- rbind(
  data.frame(cross_main[, c("transfer_label", "tool_label", "background_condition", "spike_fraction")],
             metric = "Precision", state = "Unfiltered", value = cross_main$precision_before),
  data.frame(cross_main[, c("transfer_label", "tool_label", "background_condition", "spike_fraction")],
             metric = "Precision", state = "Filtered", value = cross_main$precision_after),
  data.frame(cross_main[, c("transfer_label", "tool_label", "background_condition", "spike_fraction")],
             metric = "F1", state = "Unfiltered", value = cross_main$f1_before),
  data.frame(cross_main[, c("transfer_label", "tool_label", "background_condition", "spike_fraction")],
             metric = "F1", state = "Filtered", value = cross_main$f1_after)
)
d_sum <- aggregate(value ~ transfer_label + tool_label + metric + state,
                   data = d_long,
                   FUN = function(x) median(x, na.rm = TRUE))
d_sum$state <- factor(d_sum$state, levels = c("Unfiltered", "Filtered"))
d_sum$metric <- factor(d_sum$metric, levels = c("F1", "Precision"))

pD <- ggplot(d_sum, aes(x = state, y = value, colour = tool_label, group = tool_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  facet_grid(metric ~ transfer_label) +
  scale_colour_manual(values = profiler_cols, name = "Profiler") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  labs(
    title = "D. Part of the filtering benefit transfers between cohorts",
    subtitle = "Frozen artefact lists learned in one cohort and applied to the other",
    x = NULL,
    y = "Median performance"
  ) + theme_manuscript

# Supplementary threshold sensitivity for both internal and cross-study filtering.
int_thr <- aggregate(f1_after ~ tool_label + threshold,
                     data = internal_results,
                     FUN = function(x) median(x, na.rm = TRUE))
int_thr$setting <- "Within benchmark"
cross_thr <- aggregate(f1_after ~ tool_label + threshold,
                       data = crossstudy_results,
                       FUN = function(x) median(x, na.rm = TRUE))
cross_thr$setting <- "Cross-study transfer"
thr_plot <- rbind(int_thr, cross_thr)

pS <- ggplot(thr_plot, aes(x = 100 * threshold, y = f1_after,
                           colour = tool_label, group = tool_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.1) +
  facet_wrap(~ setting) +
  scale_colour_manual(values = profiler_cols, name = "Profiler") +
  scale_x_continuous(breaks = 100 * thresholds, labels = paste0(100 * thresholds, "%")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  labs(x = "Artefact threshold", y = "Median filtered F1") +
  theme_manuscript +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())

# Save title-free individual panels.
strip_titles <- function(p) p + labs(title = NULL, subtitle = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())

ggsave(file.path(outdir, "panel_A_internal_artifact_enrichment.pdf"), strip_titles(pA), width = 6.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(outdir, "panel_B_internal_offtarget_reduction.pdf"), strip_titles(pB), width = 8.5, height = 4.2, device = cairo_pdf)
ggsave(file.path(outdir, "panel_C_internal_biomarker_performance.pdf"), strip_titles(pC), width = 9.5, height = 4.8, device = cairo_pdf)
ggsave(file.path(outdir, "panel_D_crossstudy_transfer_summary.pdf"), strip_titles(pD), width = 10.0, height = 4.8, device = cairo_pdf)
ggsave(file.path(outdir, "supplementary_threshold_sensitivity.pdf"), pS, width = 8.5, height = 4.2, device = cairo_pdf)

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- (pA | pB) / pC / pD + patchwork::plot_layout(heights = c(0.95, 1.0, 1.0))
  ggsave(file.path(outdir, "manuscript_internal_and_crossstudy_artefact_filter_validation.pdf"),
         combined, width = 7.1, height = 8.4, device = cairo_pdf)
  ggsave(file.path(outdir, "manuscript_internal_and_crossstudy_artefact_filter_validation.png"),
         combined, width = 7.1, height = 8.4, dpi = 450)
} else {
  message("[WARN] Package 'patchwork' is unavailable; individual panels were saved, but not the composite figure.")
}

notes <- c(
  "Internal and cross-study abundance-artifact filtering validation",
  paste0("Trace: ", trace_file),
  paste0("DA results: ", da_file),
  paste0("Main threshold: ", main_threshold),
  paste0("Sensitivity thresholds: ", paste(thresholds, collapse = ", ")),
  paste0("Minimum samples per artifact estimate: ", min_samples),
  paste0("Minimum median expected abundance: ", min_median_expected),
  paste0("Minimum fraction recurrence: ", min_fraction_recurrence),
  "Panels A-C use pooled within-benchmark artefact lists and therefore represent internal proof-of-concept performance.",
  "Panel D uses frozen artefact lists learned in one cohort and applied to the other.",
  "The script applies exclusion lists to already significant original DA calls; it does not rerun MaAsLin2 after filtering.",
  "Implanted targets and aliases are protected from exclusion by design; therefore target recall is not plotted as an empirical filtering outcome."
)
writeLines(notes, file.path(outdir, "analysis_notes.txt"))

message("[INFO] Done. Outputs written to: ", outdir)
