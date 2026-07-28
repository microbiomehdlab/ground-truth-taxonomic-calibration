#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE)))
  getwd()
}

script_dir <- get_script_dir()
source(file.path(script_dir, "..", "R", "common_utils.R"))

option_list <- list(
  make_option("--metrics_dir", type = "character"),
  make_option("--maaslin_dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--plot_filter_modes", type = "character", default = "original"),
  make_option("--focus_mode", type = "character", default = "independent"),
  make_option("--conditions", type = "character", default = "ALL"),
  make_option("--studies", type = "character", default = "ALL"),
  make_option("--detection_field", type = "character", default = "positive",
              help = "Use 'positive' for member_detected_positive or 'any' for member_detected_any [default %default]"),
  make_option("--min_detection_rate", type = "double", default = 0.5,
              help = "Smallest spike fraction is called detectable when detection rate >= this threshold [default %default]"),
  make_option("--threshold_label_style", type = "character", default = "auto",
              help = "auto, spike_label, or member_taxon [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$metrics_dir) || is.null(opt$maaslin_dir) || is.null(opt$outdir)) {
  stop("Required: --metrics_dir --maaslin_dir --outdir", call. = FALSE)
}

ensure_dir(opt$outdir)

read_optional_table <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  safe_read_table_auto(path)
}

read_metric_file <- function(metrics_dir, filename) {
  root <- file.path(metrics_dir, filename)
  if (file.exists(root)) return(safe_read_table_auto(root))
  subs <- list.dirs(metrics_dir, recursive = FALSE, full.names = TRUE)
  hits <- file.path(subs, filename)
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0) return(tibble::tibble())
  dplyr::bind_rows(lapply(hits, safe_read_table_auto))
}

parse_csv_arg <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character())
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

subset_filter_modes <- function(df, keep_modes, label = "table") {
  if (!nrow(df)) return(df)
  if (!"filter_mode" %in% names(df)) {
    df$filter_mode <- "original"
    return(df)
  }
  keep_all <- length(keep_modes) == 0 || any(tolower(keep_modes) %in% c("all", "*"))
  if (keep_all) return(df)
  out <- df %>% dplyr::filter(.data$filter_mode %in% keep_modes)
  if (!nrow(out)) stop(sprintf("After --plot_filter_modes, no rows remain in %s", label), call. = FALSE)
  out
}

subset_string_values <- function(df, column, requested, label = "values") {
  if (!nrow(df) || !(column %in% names(df))) return(df)
  keep_all <- length(requested) == 0 || any(toupper(requested) == "ALL")
  if (keep_all) return(df)
  out <- df %>% dplyr::filter(.data[[column]] %in% requested)
  if (!nrow(out)) stop(sprintf("After filtering %s, no rows remain", label), call. = FALSE)
  out
}

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", colour = "#D9D9D9"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.3),
      plot.caption = ggplot2::element_text(size = base_size - 2, colour = "#5A5A5A"),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

save_pub <- function(p, path, width = 8, height = 5) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, p, width = width, height = height, units = "in", dpi = 320, bg = "white")
}

clip01 <- function(x, lo = 0, hi = 100) pmin(pmax(x, lo), hi)

score_higher_better <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  lo <- suppressWarnings(as.numeric(stats::quantile(x[ok], probs = 0.05, na.rm = TRUE, names = FALSE)))
  hi <- suppressWarnings(as.numeric(stats::quantile(x[ok], probs = 0.95, na.rm = TRUE, names = FALSE)))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    out[ok] <- 50
    return(out)
  }
  out[ok] <- 100 * (pmin(pmax(x[ok], lo), hi) - lo) / (hi - lo)
  clip01(out)
}

score_lower_better <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  lo <- suppressWarnings(as.numeric(stats::quantile(x[ok], probs = 0.05, na.rm = TRUE, names = FALSE)))
  hi <- suppressWarnings(as.numeric(stats::quantile(x[ok], probs = 0.95, na.rm = TRUE, names = FALSE)))
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    out[ok] <- 50
    return(out)
  }
  out[ok] <- 100 * (hi - pmin(pmax(x[ok], lo), hi)) / (hi - lo)
  clip01(out)
}

score_ideal_one <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  if (!any(ok)) return(out)
  dev <- abs(log2(x[ok]))
  hi <- suppressWarnings(as.numeric(stats::quantile(dev, probs = 0.95, na.rm = TRUE, names = FALSE)))
  if (!is.finite(hi) || hi <= 0) {
    out[ok] <- 100
    return(out)
  }
  dev_clip <- pmin(dev, hi)
  out[ok] <- 100 * (hi - dev_clip) / hi
  clip01(out)
}

format_label <- function(metric, value) {
  out <- rep("", length(value))
  ok <- is.finite(value)
  if (!any(ok)) return(out)
  m <- as.character(metric)
  v <- as.numeric(value)
  
  out[ok & m == "baseline_log_raw"] <- sprintf("%.1f", v[ok & m == "baseline_log_raw"])
  out[ok & m %in% c("tax_detect_raw", "bio_detect_raw")] <- sprintf("%.0f", 100 * v[ok & m %in% c("tax_detect_raw", "bio_detect_raw")])
  out[ok & m == "tax_recovery_raw"] <- sprintf("%.2f", v[ok & m == "tax_recovery_raw"])
  out[ok & m == "recovery_var_raw"] <- sprintf("%.2f", v[ok & m == "recovery_var_raw"])
  out[ok & m == "fp_taxa_raw"] <- sprintf("%.1f", v[ok & m == "fp_taxa_raw"])
  out[ok & m == "q_strength_raw"] <- sprintf("%.1f", v[ok & m == "q_strength_raw"])
  out[ok & m == "rank_raw"] <- sprintf("%.1f", v[ok & m == "rank_raw"])
  out
}

base_heat <- function(dat, x_col, label_col, fill_col, title_txt, show_y = FALSE,
                      palette_type = c("baseline", "favorability")) {
  palette_type <- match.arg(palette_type)
  p <- ggplot(dat, aes(x = .data[[x_col]], y = .data$spike_label, fill = .data[[fill_col]])) +
    geom_tile(width = 0.96, height = 0.96, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = .data[[label_col]]), size = 3.0) +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme_pub() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank()
    )
  
  if (palette_type == "baseline") {
    p <- p +
      scale_fill_gradient(
        low = "#F3EEF7",
        high = "#7A5195",
        name = "Baseline\nabundance",
        guide = guide_colorbar(barheight = unit(3.5, "cm"))
      )
  } else {
    p <- p +
      scale_fill_gradient(
        low = "#F2F2F2",
        high = "#2166AC",
        limits = c(0, 100),
        breaks = c(0, 100),
        labels = c("Less favorable", "More favorable"),
        oob = scales::squish,
        name = "Closer to\nideal",
        guide = guide_colorbar(barheight = unit(3.5, "cm"))
      )
  }
  p
}

short_member_label <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  keep <- !is.na(x) & nzchar(trimws(x))
  if (!any(keep)) return(out)
  vals <- trimws(x[keep])
  parts <- strsplit(vals, "[[:space:]]+")
  short <- vapply(parts, function(tok) {
    if (length(tok) == 0) return(NA_character_)
    genus <- tok[[1]]
    species <- if (length(tok) >= 2) tok[[2]] else tok[[1]]
    paste0(substr(genus, 1, 1), substr(species, 1, 4))
  }, character(1))
  out[keep] <- short
  out
}

choose_row_label <- function(df, label_style = "auto", focus_mode = "independent") {
  if (label_style == "spike_label" && "spike_label" %in% names(df)) return(as.character(df$spike_label))
  if (label_style == "member_taxon" && "member_taxon" %in% names(df)) return(as.character(df$member_taxon))
  if (focus_mode == "independent" && "spike_label" %in% names(df)) return(as.character(df$spike_label))
  if ("member_taxon" %in% names(df)) return(short_member_label(df$member_taxon))
  if ("spike_label" %in% names(df)) return(as.character(df$spike_label))
  rep(NA_character_, nrow(df))
}

frac_to_label <- function(x) {
  ifelse(is.finite(x), scales::label_number(scale = 100, accuracy = 0.01, suffix = "%")(x), "NR")
}

copy_first_existing_col <- function(df, target, candidates) {
  if (target %in% names(df)) return(df)
  hits <- candidates[candidates %in% names(df)]
  if (!length(hits)) {
    stop(sprintf(
      "Could not create '%s'; none of these columns exist: %s",
      target, paste(candidates, collapse = ", ")
    ), call. = FALSE)
  }
  df[[target]] <- df[[hits[[1]]]]
  df
}

copy_or_guess_col <- function(df, target, candidates = character(), pattern = NULL, exclude = character()) {
  if (target %in% names(df)) return(df)
  
  hits <- candidates[candidates %in% names(df)]
  if (!length(hits) && !is.null(pattern) && nzchar(pattern)) {
    hits <- names(df)[grepl(pattern, names(df), ignore.case = TRUE)]
    if (length(exclude)) hits <- setdiff(hits, exclude)
  }
  
  if (!length(hits)) {
    stop(sprintf(
      paste0(
        "Could not create '%s'.\n",
        "Tried explicit candidates: %s\n",
        "Regex used: %s\n",
        "Available columns are: %s"
      ),
      target,
      paste(candidates, collapse = ", "),
      ifelse(is.null(pattern), "<none>", pattern),
      paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }
  
  if (length(hits) > 1) {
    stop(sprintf(
      paste0(
        "Could not create '%s' because multiple candidate columns matched: %s\n",
        "Please edit the script and set the right source column explicitly.\n",
        "Available columns are: %s"
      ),
      target,
      paste(hits, collapse = ", "),
      paste(names(df), collapse = ", ")
    ), call. = FALSE)
  }
  
  df[[target]] <- df[[hits[[1]]]]
  df
}

slugify <- function(x) {
  x <- as.character(x)
  x <- gsub("%", "pct", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

plot_threshold_heatmap <- function(thresh_df, outpath, title_txt, subtitle_txt, row_order_pref = NULL, by_study = FALSE, caption_txt = NULL) {
  if (!nrow(thresh_df)) return(invisible(NULL))
  
  if (!is.null(row_order_pref)) {
    present <- unique(as.character(thresh_df$row_label))
    row_order <- row_order_pref[row_order_pref %in% present]
    missing <- setdiff(sort(present), row_order)
    row_order <- c(row_order, missing)
  } else {
    row_order <- thresh_df %>%
      dplyr::group_by(.data$row_label) %>%
      dplyr::summarise(
        mean_threshold = mean(.data$threshold_fraction, na.rm = TRUE),
        n_detected = sum(is.finite(.data$threshold_fraction)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(mean_threshold = ifelse(is.nan(.data$mean_threshold), Inf, .data$mean_threshold)) %>%
      dplyr::arrange(.data$mean_threshold, dplyr::desc(.data$n_detected), .data$row_label) %>%
      dplyr::pull(.data$row_label)
  }
  
  plot_df <- thresh_df %>%
    dplyr::mutate(
      row_label = factor(.data$row_label, levels = rev(row_order)),
      threshold_label = frac_to_label(.data$threshold_fraction),
      ease_score = ifelse(is.finite(.data$threshold_fraction), -log10(.data$threshold_fraction), NA_real_)
    )
  
  max_ease <- suppressWarnings(max(plot_df$ease_score, na.rm = TRUE))
  if (!is.finite(max_ease)) max_ease <- 1
  
  p <- ggplot(plot_df, aes(x = .data$background_condition, y = .data$row_label, fill = .data$ease_score)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = .data$threshold_label), size = 3.0) +
    facet_grid(if (by_study) background_study ~ tool else . ~ tool, scales = "free_y", space = "free_y") +
    scale_fill_gradient(
      low = "#E8E8E8",
      high = "#2166AC",
      limits = c(0, max_ease),
      breaks = c(0, max_ease),
      labels = c("Harder / higher fraction", "Easier / lower fraction"),
      na.value = "#F7F7F7",
      oob = scales::squish,
      name = "Detection\nthreshold"
    ) +
    labs(
      x = "Background condition",
      y = "Spiked taxon",
      title = title_txt,
      subtitle = subtitle_txt,
      caption = if (is.null(caption_txt)) paste(
        "Numbers show the minimum spike fraction where biomarker detection rate reached the requested threshold.",
        "NR = not reached within tested spike fractions.",
        "Darker blue means detectable at lower spike fractions."
      ) else caption_txt
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_pub(p, outpath, width = if (by_study) 12 else 10.5, height = max(5.8, 0.42 * dplyr::n_distinct(plot_df$row_label) + if (by_study) 3.5 else 2.6))
}

filter_modes_keep <- parse_csv_arg(opt$plot_filter_modes)
condition_keep <- parse_csv_arg(opt$conditions)
study_keep <- parse_csv_arg(opt$studies)

target_member_metrics <- read_metric_file(opt$metrics_dir, "target_member_errors_with_condition.csv")
combined_metrics <- read_metric_file(opt$metrics_dir, "combined_with_condition.csv")
member_detection <- read_optional_table(file.path(opt$maaslin_dir, "maaslin_member_detection_ALLFILTERS.csv"))
sig_features <- read_optional_table(file.path(opt$maaslin_dir, "maaslin_significant_features_ALLFILTERS.csv"))

if (!nrow(target_member_metrics)) stop("No target_member_errors_with_condition.csv found", call. = FALSE)
if (!nrow(combined_metrics)) stop("No combined_with_condition.csv found", call. = FALSE)
if (!nrow(member_detection)) stop("No maaslin_member_detection_ALLFILTERS.csv found", call. = FALSE)

if (!"spike_mode" %in% names(target_member_metrics)) target_member_metrics$spike_mode <- NA_character_
if (!"Study" %in% names(target_member_metrics)) target_member_metrics$Study <- NA_character_
if (!"spike_mode" %in% names(combined_metrics)) combined_metrics$spike_mode <- NA_character_
if (!"Study" %in% names(combined_metrics)) combined_metrics$Study <- NA_character_
if (!"spike_mode" %in% names(member_detection)) member_detection$spike_mode <- NA_character_
if (!"filter_mode" %in% names(member_detection)) member_detection$filter_mode <- "original"
if (!"background_study" %in% names(member_detection)) member_detection$background_study <- NA_character_

member_detection <- subset_filter_modes(member_detection, filter_modes_keep, "maaslin_member_detection_ALLFILTERS.csv")
if (nrow(sig_features)) {
  if (!"spike_mode" %in% names(sig_features)) sig_features$spike_mode <- NA_character_
  if (!"filter_mode" %in% names(sig_features)) sig_features$filter_mode <- "original"
  if (!"background_study" %in% names(sig_features)) sig_features$background_study <- NA_character_
  sig_features <- subset_filter_modes(sig_features, filter_modes_keep, "maaslin_significant_features_ALLFILTERS.csv")
}

target_member_metrics <- subset_string_values(target_member_metrics, "Target_Condition", condition_keep, "conditions")
target_member_metrics <- subset_string_values(target_member_metrics, "Study", study_keep, "studies")
combined_metrics <- subset_string_values(combined_metrics, "Target_Condition", condition_keep, "conditions")
combined_metrics <- subset_string_values(combined_metrics, "Study", study_keep, "studies")
member_detection <- subset_string_values(member_detection, "background_condition", condition_keep, "conditions")
member_detection <- subset_string_values(member_detection, "background_study", study_keep, "studies")
if (nrow(sig_features)) {
  sig_features <- subset_string_values(sig_features, "background_condition", condition_keep, "conditions")
  sig_features <- subset_string_values(sig_features, "background_study", study_keep, "studies")
}

focus_member <- target_member_metrics %>% dplyr::filter(.data$spike_mode == opt$focus_mode)
focus_combined <- combined_metrics %>% dplyr::filter(.data$spike_mode == opt$focus_mode)
focus_maaslin <- member_detection %>% dplyr::filter(.data$spike_mode == opt$focus_mode)
focus_sig <- if (nrow(sig_features)) sig_features %>% dplyr::filter(.data$spike_mode == opt$focus_mode) else tibble::tibble()

if (!nrow(focus_member)) stop(sprintf("No target-member metric rows for focus_mode = %s", opt$focus_mode), call. = FALSE)
if (!nrow(focus_maaslin)) stop(sprintf("No MaAsLin member rows for focus_mode = %s", opt$focus_mode), call. = FALSE)

focus_member <- copy_first_existing_col(focus_member, "background_condition", c("background_condition", "Target_Condition"))
focus_member <- copy_first_existing_col(focus_member, "background_study", c("background_study", "Study"))
focus_member <- copy_or_guess_col(
  focus_member,
  "spike_fraction",
  candidates = c(
    "spike_fraction", "Spike_Fraction", "fraction", "spike_frac", "SpikeFrac",
    "spike_proportion", "Spike_Proportion", "proportion", "prop",
    "SpikeInFraction", "spikein_fraction", "spike_in_fraction",
    "expected_fraction", "Expected_Fraction",
    "spike_relative_abundance", "expected_relative_abundance", "relative_abundance"
  ),
  pattern = "(spike.*(fraction|frac|prop|proportion|abundance)|expected.*(fraction|abundance)|(fraction|frac|prop|proportion|abundance).*spike)",
  exclude = c("spike_mode", "spike_label")
)

focus_combined <- copy_first_existing_col(focus_combined, "background_condition", c("background_condition", "Target_Condition"))
focus_combined <- copy_first_existing_col(focus_combined, "background_study", c("background_study", "Study"))
focus_combined <- copy_or_guess_col(
  focus_combined,
  "spike_fraction",
  candidates = c(
    "spike_fraction", "Spike_Fraction", "fraction", "spike_frac", "SpikeFrac",
    "spike_proportion", "Spike_Proportion", "proportion", "prop",
    "SpikeInFraction", "spikein_fraction", "spike_in_fraction",
    "expected_fraction", "Expected_Fraction",
    "spike_relative_abundance", "expected_relative_abundance", "relative_abundance"
  ),
  pattern = "(spike.*(fraction|frac|prop|proportion|abundance)|expected.*(fraction|abundance)|(fraction|frac|prop|proportion|abundance).*spike)",
  exclude = c("spike_mode", "spike_label")
)

focus_maaslin <- copy_first_existing_col(focus_maaslin, "background_condition", c("background_condition"))
focus_maaslin <- copy_first_existing_col(focus_maaslin, "background_study", c("background_study"))
focus_maaslin <- copy_or_guess_col(
  focus_maaslin,
  "spike_fraction",
  candidates = c(
    "spike_fraction", "Spike_Fraction", "fraction", "spike_frac", "SpikeFrac",
    "spike_proportion", "Spike_Proportion", "proportion", "prop",
    "SpikeInFraction", "spikein_fraction", "spike_in_fraction",
    "expected_fraction", "Expected_Fraction",
    "spike_relative_abundance", "expected_relative_abundance", "relative_abundance"
  ),
  pattern = "(spike.*(fraction|frac|prop|proportion|abundance)|expected.*(fraction|abundance)|(fraction|frac|prop|proportion|abundance).*spike)",
  exclude = c("spike_mode", "spike_label")
)

if (nrow(focus_sig)) {
  focus_sig <- copy_first_existing_col(focus_sig, "background_condition", c("background_condition"))
  focus_sig <- copy_first_existing_col(focus_sig, "background_study", c("background_study"))
  focus_sig <- copy_or_guess_col(
    focus_sig,
    "spike_fraction",
    candidates = c(
      "spike_fraction", "Spike_Fraction", "fraction", "spike_frac", "SpikeFrac",
      "spike_proportion", "Spike_Proportion", "proportion", "prop",
      "SpikeInFraction", "spikein_fraction", "spike_in_fraction",
      "expected_fraction", "Expected_Fraction",
      "spike_relative_abundance", "expected_relative_abundance", "relative_abundance"
    ),
    pattern = "(spike.*(fraction|frac|prop|proportion|abundance)|expected.*(fraction|abundance)|(fraction|frac|prop|proportion|abundance).*spike)",
    exclude = c("spike_mode", "spike_label")
  )
}

detect_col <- if (tolower(opt$detection_field) == "any") "member_detected_any" else "member_detected_positive"
if (!(detect_col %in% names(focus_maaslin))) {
  stop(sprintf("Requested detection field not present: %s", detect_col), call. = FALSE)
}

metric_taxon <- focus_member %>%
  dplyr::group_by(.data$tool, .data$spike_label) %>%
  dplyr::summarise(
    baseline_log_raw = if ("baseline" %in% names(.)) median(log10(.data$baseline + 1e-6), na.rm = TRUE) else NA_real_,
    tax_detect_raw = if ("detected" %in% names(.)) mean(.data$detected, na.rm = TRUE) else NA_real_,
    tax_recovery_raw = if ("observed_over_expected" %in% names(.)) median(.data$observed_over_expected, na.rm = TRUE) else NA_real_,
    recovery_var_raw = if ("observed_over_expected" %in% names(.)) stats::sd(log2(pmax(.data$observed_over_expected, 1e-6)), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

metric_fp <- focus_combined %>%
  dplyr::group_by(.data$tool, .data$spike_label) %>%
  dplyr::summarise(
    fp_taxa_raw = if ("n_false_positive_taxa" %in% names(.)) mean(.data$n_false_positive_taxa, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

bio_summary <- focus_maaslin %>%
  dplyr::group_by(.data$tool, .data$spike_label) %>%
  dplyr::summarise(
    bio_detect_raw = mean(as.numeric(dplyr::coalesce(.data[[detect_col]], FALSE)), na.rm = TRUE),
    q_strength_raw = median(-log10(dplyr::coalesce(.data$min_q_target, 1)), na.rm = TRUE),
    .groups = "drop"
  )

rank_summary <- tibble::tibble(tool = character(), spike_label = character(), rank_raw = numeric())
if (nrow(focus_sig) > 0 && all(c("tool", "spike_label", "feature") %in% names(focus_sig))) {
  q_col <- if ("qval" %in% names(focus_sig)) "qval" else if ("q" %in% names(focus_sig)) "q" else NA_character_
  coef_col <- if ("coef" %in% names(focus_sig)) "coef" else if ("coefficient" %in% names(focus_sig)) "coefficient" else NA_character_
  if (!is.na(q_col) && q_col %in% names(focus_sig)) {
    rank_summary <- focus_sig %>%
      dplyr::mutate(
        q_use = .data[[q_col]],
        coef_use = if (!is.na(coef_col) && coef_col %in% names(.)) .data[[coef_col]] else 1,
        is_positive2 = if ("is_positive" %in% names(.)) .data$is_positive else (.data$coef_use > 0),
        is_target2 = if ("is_target" %in% names(.)) .data$is_target else (.data$feature == .data$spike_label)
      ) %>%
      dplyr::filter(.data$is_positive2) %>%
      dplyr::group_by(.data$background_condition, .data$background_study, .data$tool, .data$spike_mode, .data$spike_label, .data$spike_fraction, .data$filter_mode) %>%
      dplyr::arrange(.data$q_use, dplyr::desc(abs(.data$coef_use)), .by_group = TRUE) %>%
      dplyr::mutate(rank_in_run = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.data$is_target2) %>%
      dplyr::group_by(.data$tool, .data$spike_label) %>%
      dplyr::summarise(rank_raw = median(.data$rank_in_run, na.rm = TRUE), .groups = "drop")
  }
}

summary_wide <- metric_taxon %>%
  dplyr::left_join(metric_fp, by = c("tool", "spike_label")) %>%
  dplyr::left_join(bio_summary, by = c("tool", "spike_label")) %>%
  dplyr::left_join(rank_summary, by = c("tool", "spike_label"))

if (!nrow(summary_wide)) stop("No rows left after summarising", call. = FALSE)

summary_scores <- summary_wide %>%
  dplyr::mutate(
    tax_detect_score = score_higher_better(.data$tax_detect_raw),
    tax_recovery_score = score_ideal_one(.data$tax_recovery_raw),
    recovery_var_score = score_lower_better(.data$recovery_var_raw),
    fp_taxa_score = score_lower_better(.data$fp_taxa_raw),
    bio_detect_score = score_higher_better(.data$bio_detect_raw),
    q_strength_score = score_higher_better(.data$q_strength_raw),
    rank_score = score_lower_better(.data$rank_raw)
  )

row_order <- summary_scores %>%
  dplyr::group_by(.data$spike_label) %>%
  dplyr::summarise(
    baseline_order = mean(.data$baseline_log_raw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$baseline_order), .data$spike_label) %>%
  dplyr::pull(.data$spike_label)

tool_labels <- c(metaphlan4 = "MetaPhlAn 4", kraken2_bracken = "Kraken2 + Bracken")
tool_order <- c("kraken2_bracken", "metaphlan4")
tool_order <- tool_order[tool_order %in% unique(summary_scores$tool)]

summary_scores <- summary_scores %>%
  dplyr::mutate(
    spike_label = factor(.data$spike_label, levels = rev(row_order)),
    tool_pretty = dplyr::recode(as.character(.data$tool), !!!tool_labels, .default = as.character(.data$tool))
  )

baseline_long <- summary_scores %>%
  dplyr::transmute(
    tool = .data$tool,
    tool_pretty = .data$tool_pretty,
    spike_label = .data$spike_label,
    metric = "Baseline\nlog10(x+1e-6)",
    fill_value = .data$baseline_log_raw,
    raw_label = format_label("baseline_log_raw", .data$baseline_log_raw)
  )

taxonomy_long <- dplyr::bind_rows(
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Tax detect\n(100% best)", fill_value = .data$tax_detect_score, raw_label = format_label("tax_detect_raw", .data$tax_detect_raw)),
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Tax recovery\n(1 best)", fill_value = .data$tax_recovery_score, raw_label = format_label("tax_recovery_raw", .data$tax_recovery_raw)),
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Recovery variability\n(0 best)", fill_value = .data$recovery_var_score, raw_label = format_label("recovery_var_raw", .data$recovery_var_raw)),
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Low FP taxa\n(0 best)", fill_value = .data$fp_taxa_score, raw_label = format_label("fp_taxa_raw", .data$fp_taxa_raw))
)

biomarker_long <- dplyr::bind_rows(
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Bio detect\n(100% best)", fill_value = .data$bio_detect_score, raw_label = format_label("bio_detect_raw", .data$bio_detect_raw)),
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Bio -log10(q)\n(high best)", fill_value = .data$q_strength_score, raw_label = format_label("q_strength_raw", .data$q_strength_raw)),
  summary_scores %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, spike_label = .data$spike_label, metric = "Bio rank\n(low best)", fill_value = .data$rank_score, raw_label = format_label("rank_raw", .data$rank_raw))
)

baseline_long$metric <- factor(baseline_long$metric, levels = "Baseline\nlog10(x+1e-6)")
taxonomy_long$metric <- factor(taxonomy_long$metric, levels = c("Tax detect\n(100% best)", "Tax recovery\n(1 best)", "Recovery variability\n(0 best)", "Low FP taxa\n(0 best)"))
biomarker_long$metric <- factor(biomarker_long$metric, levels = c("Bio detect\n(100% best)", "Bio -log10(q)\n(high best)", "Bio rank\n(low best)"))

write_csv_safe(summary_scores, file.path(opt$outdir, "species_driver_summary_wide_v6.csv"))
write_csv_safe(baseline_long, file.path(opt$outdir, "species_driver_baseline_panel_v6.csv"))
write_csv_safe(taxonomy_long, file.path(opt$outdir, "species_driver_taxonomy_panel_v6.csv"))
write_csv_safe(biomarker_long, file.path(opt$outdir, "species_driver_biomarker_panel_v6.csv"))

# ------------------------------------------------------------------------------
# Detailed driver summaries: condition x spike_fraction
# ------------------------------------------------------------------------------

detail_group <- c("tool", "background_condition", "spike_fraction", "spike_label")

metric_taxon_detail <- focus_member %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(detail_group))) %>%
  dplyr::summarise(
    baseline_log_raw = if ("baseline" %in% names(.)) median(log10(.data$baseline + 1e-6), na.rm = TRUE) else NA_real_,
    tax_detect_raw = if ("detected" %in% names(.)) mean(.data$detected, na.rm = TRUE) else NA_real_,
    tax_recovery_raw = if ("observed_over_expected" %in% names(.)) median(.data$observed_over_expected, na.rm = TRUE) else NA_real_,
    recovery_var_raw = if ("observed_over_expected" %in% names(.)) stats::sd(log2(pmax(.data$observed_over_expected, 1e-6)), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

metric_fp_detail <- focus_combined %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(detail_group))) %>%
  dplyr::summarise(
    fp_taxa_raw = if ("n_false_positive_taxa" %in% names(.)) mean(.data$n_false_positive_taxa, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

bio_summary_detail <- focus_maaslin %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(detail_group))) %>%
  dplyr::summarise(
    bio_detect_raw = mean(as.numeric(dplyr::coalesce(.data[[detect_col]], FALSE)), na.rm = TRUE),
    q_strength_raw = median(-log10(dplyr::coalesce(.data$min_q_target, 1)), na.rm = TRUE),
    .groups = "drop"
  )

rank_summary_detail <- tibble::tibble(
  tool = character(),
  background_condition = character(),
  spike_fraction = numeric(),
  spike_label = character(),
  rank_raw = numeric()
)

if (nrow(focus_sig) > 0 && all(c("tool", "spike_label", "feature", "background_condition", "spike_fraction") %in% names(focus_sig))) {
  q_col <- if ("qval" %in% names(focus_sig)) "qval" else if ("q" %in% names(focus_sig)) "q" else NA_character_
  coef_col <- if ("coef" %in% names(focus_sig)) "coef" else if ("coefficient" %in% names(focus_sig)) "coefficient" else NA_character_
  
  if (!is.na(q_col) && q_col %in% names(focus_sig)) {
    rank_summary_detail <- focus_sig %>%
      dplyr::mutate(
        q_use = .data[[q_col]],
        coef_use = if (!is.na(coef_col) && coef_col %in% names(.)) .data[[coef_col]] else 1,
        is_positive2 = if ("is_positive" %in% names(.)) .data$is_positive else (.data$coef_use > 0),
        is_target2 = if ("is_target" %in% names(.)) .data$is_target else (.data$feature == .data$spike_label)
      ) %>%
      dplyr::filter(.data$is_positive2) %>%
      dplyr::group_by(.data$background_condition, .data$background_study, .data$tool, .data$spike_mode, .data$spike_label, .data$spike_fraction, .data$filter_mode) %>%
      dplyr::arrange(.data$q_use, dplyr::desc(abs(.data$coef_use)), .by_group = TRUE) %>%
      dplyr::mutate(rank_in_run = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.data$is_target2) %>%
      dplyr::group_by(.data$tool, .data$background_condition, .data$spike_fraction, .data$spike_label) %>%
      dplyr::summarise(rank_raw = median(.data$rank_in_run, na.rm = TRUE), .groups = "drop")
  }
}

summary_detail_wide <- metric_taxon_detail %>%
  dplyr::left_join(metric_fp_detail, by = c("tool", "background_condition", "spike_fraction", "spike_label")) %>%
  dplyr::left_join(bio_summary_detail, by = c("tool", "background_condition", "spike_fraction", "spike_label")) %>%
  dplyr::left_join(rank_summary_detail, by = c("tool", "background_condition", "spike_fraction", "spike_label"))

summary_scores_detail <- summary_detail_wide %>%
  dplyr::mutate(
    tax_detect_score = score_higher_better(.data$tax_detect_raw),
    tax_recovery_score = score_ideal_one(.data$tax_recovery_raw),
    recovery_var_score = score_lower_better(.data$recovery_var_raw),
    fp_taxa_score = score_lower_better(.data$fp_taxa_raw),
    bio_detect_score = score_higher_better(.data$bio_detect_raw),
    q_strength_score = score_higher_better(.data$q_strength_raw),
    rank_score = score_lower_better(.data$rank_raw),
    spike_label = factor(.data$spike_label, levels = rev(row_order)),
    tool_pretty = dplyr::recode(as.character(.data$tool), !!!tool_labels, .default = as.character(.data$tool))
  )

baseline_long_detail <- summary_scores_detail %>%
  dplyr::transmute(
    tool = .data$tool,
    tool_pretty = .data$tool_pretty,
    background_condition = .data$background_condition,
    spike_fraction = .data$spike_fraction,
    spike_label = .data$spike_label,
    metric = "Baseline\nlog10(x+1e-6)",
    fill_value = .data$baseline_log_raw,
    raw_label = format_label("baseline_log_raw", .data$baseline_log_raw)
  )

taxonomy_long_detail <- dplyr::bind_rows(
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Tax detect\n(100% best)", fill_value = .data$tax_detect_score, raw_label = format_label("tax_detect_raw", .data$tax_detect_raw)),
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Tax recovery\n(1 best)", fill_value = .data$tax_recovery_score, raw_label = format_label("tax_recovery_raw", .data$tax_recovery_raw)),
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Recovery variability\n(0 best)", fill_value = .data$recovery_var_score, raw_label = format_label("recovery_var_raw", .data$recovery_var_raw)),
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Low FP taxa\n(0 best)", fill_value = .data$fp_taxa_score, raw_label = format_label("fp_taxa_raw", .data$fp_taxa_raw))
)

biomarker_long_detail <- dplyr::bind_rows(
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Bio detect\n(100% best)", fill_value = .data$bio_detect_score, raw_label = format_label("bio_detect_raw", .data$bio_detect_raw)),
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Bio -log10(q)\n(high best)", fill_value = .data$q_strength_score, raw_label = format_label("q_strength_raw", .data$q_strength_raw)),
  summary_scores_detail %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Bio rank\n(low best)", fill_value = .data$rank_score, raw_label = format_label("rank_raw", .data$rank_raw))
)

baseline_long_detail$metric <- factor(baseline_long_detail$metric, levels = "Baseline\nlog10(x+1e-6)")
taxonomy_long_detail$metric <- factor(taxonomy_long_detail$metric, levels = c("Tax detect\n(100% best)", "Tax recovery\n(1 best)", "Recovery variability\n(0 best)", "Low FP taxa\n(0 best)"))
biomarker_long_detail$metric <- factor(biomarker_long_detail$metric, levels = c("Bio detect\n(100% best)", "Bio -log10(q)\n(high best)", "Bio rank\n(low best)"))

write_csv_safe(summary_scores_detail, file.path(opt$outdir, "species_driver_summary_by_condition_fraction_v6.csv"))
write_csv_safe(baseline_long_detail, file.path(opt$outdir, "species_driver_baseline_by_condition_fraction_v6.csv"))
write_csv_safe(taxonomy_long_detail, file.path(opt$outdir, "species_driver_taxonomy_by_condition_fraction_v6.csv"))
write_csv_safe(biomarker_long_detail, file.path(opt$outdir, "species_driver_biomarker_by_condition_fraction_v6.csv"))

make_tool_row <- function(tool_name, show_y = FALSE) {
  tpretty <- dplyr::recode(tool_name, !!!tool_labels, .default = tool_name)
  
  p_base <- baseline_long %>%
    dplyr::filter(.data$tool == tool_name) %>%
    base_heat("metric", "raw_label", "fill_value", "Baseline context", show_y = show_y, palette_type = "baseline")
  
  p_tax <- taxonomy_long %>%
    dplyr::filter(.data$tool == tool_name) %>%
    base_heat("metric", "raw_label", "fill_value", "Taxonomy-side drivers", show_y = FALSE, palette_type = "favorability") +
    theme(legend.position = "none")
  
  p_bio <- biomarker_long %>%
    dplyr::filter(.data$tool == tool_name) %>%
    base_heat("metric", "raw_label", "fill_value", "Biomarker outcome", show_y = FALSE, palette_type = "favorability")
  
  (p_base + p_tax + p_bio + patchwork::plot_layout(widths = c(1.3, 4.4, 3.1), guides = "collect")) +
    patchwork::plot_annotation(title = tpretty) &
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

rows <- list()
if (length(tool_order) >= 1) rows[[1]] <- make_tool_row(tool_order[[1]], show_y = TRUE)
if (length(tool_order) >= 2) rows[[2]] <- make_tool_row(tool_order[[2]], show_y = TRUE)

combined_plot <- wrap_plots(rows, ncol = 1, guides = "collect") +
  patchwork::plot_annotation(
    title = sprintf("Which factors explain spike biomarker significance? (%s spikes)", opt$focus_mode),
    subtitle = paste(
      "Baseline is shown separately as context, not as a favorable or unfavorable property.",
      "The taxonomy block summarizes detectability, recovery accuracy, recovery variability, and false-positive burden.",
      "The biomarker block shows detection, q-value strength, and rank."
    ),
    caption = paste(
      "Baseline values are log10(relative abundance + 1e-6).",
      "Tax recovery is observed/expected, where 1 is ideal.",
      "Recovery variability is SD(log2(observed/expected)) across matched observations, where 0 is ideal.",
      "Low FP taxa is the absolute number of false-positive taxa.",
      "For taxonomy and biomarker blocks, color means closer to the metric-specific ideal; the printed numbers are the raw values."
    )
  ) &
  theme_pub() &
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0, margin = margin(b = 10)),
    legend.position = "right"
  )

save_pub(combined_plot, file.path(opt$outdir, "figure_species_driver_heatmap_v6.png"), width = 15.5, height = max(10.5, 0.55 * length(row_order) + 6.2))

make_tool_row_detail <- function(tool_name, condition_name, spike_frac, show_y = FALSE) {
  tpretty <- dplyr::recode(tool_name, !!!tool_labels, .default = tool_name)
  
  base_dat <- baseline_long_detail %>%
    dplyr::filter(.data$tool == tool_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  tax_dat <- taxonomy_long_detail %>%
    dplyr::filter(.data$tool == tool_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  bio_dat <- biomarker_long_detail %>%
    dplyr::filter(.data$tool == tool_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  if (!nrow(base_dat) && !nrow(tax_dat) && !nrow(bio_dat)) return(NULL)
  
  p_base <- base_dat %>%
    base_heat("metric", "raw_label", "fill_value", "Baseline context", show_y = show_y, palette_type = "baseline")
  
  p_tax <- tax_dat %>%
    base_heat("metric", "raw_label", "fill_value", "Taxonomy-side drivers", show_y = FALSE, palette_type = "favorability") +
    theme(legend.position = "none")
  
  p_bio <- bio_dat %>%
    base_heat("metric", "raw_label", "fill_value", "Biomarker outcome", show_y = FALSE, palette_type = "favorability")
  
  (p_base + p_tax + p_bio + patchwork::plot_layout(widths = c(1.3, 4.4, 3.1), guides = "collect")) +
    patchwork::plot_annotation(title = tpretty) &
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

make_detail_plot <- function(condition_name, spike_frac) {
  rows <- lapply(tool_order, function(tt) make_tool_row_detail(tt, condition_name, spike_frac, show_y = TRUE))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  
  wrap_plots(rows, ncol = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Which factors explain spike biomarker significance? (%s spikes)", opt$focus_mode),
      subtitle = paste(
        sprintf("Background condition = %s; spike fraction = %s.", condition_name, frac_to_label(spike_frac)),
        "These summaries are aggregated across studies within that condition and spike level.",
        "This is the condition-level companion to the LoD heatmaps."
      ),
      caption = paste(
        "Baseline values are log10(relative abundance + 1e-6).",
        "Tax recovery is observed/expected, where 1 is ideal.",
        "Recovery variability is SD(log2(observed/expected)) across matched observations, where 0 is ideal.",
        "Low FP taxa is the absolute number of false-positive taxa.",
        "Colors in the taxonomy and biomarker blocks show closeness to the metric-specific ideal; numbers are raw values."
      )
    ) &
    theme_pub() &
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0, margin = margin(b = 10)),
      legend.position = "right"
    )
}

detail_keys <- summary_scores_detail %>%
  dplyr::distinct(.data$background_condition, .data$spike_fraction) %>%
  dplyr::filter(!is.na(.data$background_condition), nzchar(.data$background_condition), is.finite(.data$spike_fraction)) %>%
  dplyr::arrange(.data$background_condition, .data$spike_fraction)

if (nrow(detail_keys)) {
  detail_height <- max(10.5, 0.55 * length(row_order) + 6.2)
  pdf_path <- file.path(opt$outdir, "figure_species_driver_heatmap_by_condition_fraction_v6.pdf")
  
  grDevices::pdf(pdf_path, width = 15.5, height = detail_height, onefile = TRUE)
  for (i in seq_len(nrow(detail_keys))) {
    cond_i <- detail_keys$background_condition[[i]]
    frac_i <- detail_keys$spike_fraction[[i]]
    
    p_i <- make_detail_plot(cond_i, frac_i)
    if (is.null(p_i)) next
    
    print(p_i)
    
    png_path <- file.path(
      opt$outdir,
      sprintf(
        "figure_species_driver_heatmap__cond-%s__spike-%s_v6.png",
        slugify(cond_i),
        slugify(frac_to_label(frac_i))
      )
    )
    save_pub(p_i, png_path, width = 15.5, height = detail_height)
  }
  grDevices::dev.off()
}

# ------------------------------------------------------------------------------
# Detailed driver summaries: study x condition x spike_fraction
# ------------------------------------------------------------------------------

detail_group_study <- c("tool", "background_study", "background_condition", "spike_fraction", "spike_label")

metric_taxon_detail_study <- focus_member %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(detail_group_study))) %>%
  dplyr::summarise(
    baseline_log_raw = if ("baseline" %in% names(.)) median(log10(.data$baseline + 1e-6), na.rm = TRUE) else NA_real_,
    tax_detect_raw = if ("detected" %in% names(.)) mean(.data$detected, na.rm = TRUE) else NA_real_,
    tax_recovery_raw = if ("observed_over_expected" %in% names(.)) median(.data$observed_over_expected, na.rm = TRUE) else NA_real_,
    recovery_var_raw = if ("observed_over_expected" %in% names(.)) stats::sd(log2(pmax(.data$observed_over_expected, 1e-6)), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

metric_fp_detail_study <- focus_combined %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(detail_group_study))) %>%
  dplyr::summarise(
    fp_taxa_raw = if ("n_false_positive_taxa" %in% names(.)) mean(.data$n_false_positive_taxa, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

bio_summary_detail_study <- focus_maaslin %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(detail_group_study))) %>%
  dplyr::summarise(
    bio_detect_raw = mean(as.numeric(dplyr::coalesce(.data[[detect_col]], FALSE)), na.rm = TRUE),
    q_strength_raw = median(-log10(dplyr::coalesce(.data$min_q_target, 1)), na.rm = TRUE),
    .groups = "drop"
  )

rank_summary_detail_study <- tibble::tibble(
  tool = character(),
  background_study = character(),
  background_condition = character(),
  spike_fraction = numeric(),
  spike_label = character(),
  rank_raw = numeric()
)

if (nrow(focus_sig) > 0 && all(c("tool", "spike_label", "feature", "background_study", "background_condition", "spike_fraction") %in% names(focus_sig))) {
  q_col <- if ("qval" %in% names(focus_sig)) "qval" else if ("q" %in% names(focus_sig)) "q" else NA_character_
  coef_col <- if ("coef" %in% names(focus_sig)) "coef" else if ("coefficient" %in% names(focus_sig)) "coefficient" else NA_character_
  
  if (!is.na(q_col) && q_col %in% names(focus_sig)) {
    rank_summary_detail_study <- focus_sig %>%
      dplyr::mutate(
        q_use = .data[[q_col]],
        coef_use = if (!is.na(coef_col) && coef_col %in% names(.)) .data[[coef_col]] else 1,
        is_positive2 = if ("is_positive" %in% names(.)) .data$is_positive else (.data$coef_use > 0),
        is_target2 = if ("is_target" %in% names(.)) .data$is_target else (.data$feature == .data$spike_label)
      ) %>%
      dplyr::filter(.data$is_positive2) %>%
      dplyr::group_by(.data$background_condition, .data$background_study, .data$tool, .data$spike_mode, .data$spike_label, .data$spike_fraction, .data$filter_mode) %>%
      dplyr::arrange(.data$q_use, dplyr::desc(abs(.data$coef_use)), .by_group = TRUE) %>%
      dplyr::mutate(rank_in_run = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.data$is_target2) %>%
      dplyr::group_by(.data$tool, .data$background_study, .data$background_condition, .data$spike_fraction, .data$spike_label) %>%
      dplyr::summarise(rank_raw = median(.data$rank_in_run, na.rm = TRUE), .groups = "drop")
  }
}

summary_detail_wide_study <- metric_taxon_detail_study %>%
  dplyr::left_join(metric_fp_detail_study, by = c("tool", "background_study", "background_condition", "spike_fraction", "spike_label")) %>%
  dplyr::left_join(bio_summary_detail_study, by = c("tool", "background_study", "background_condition", "spike_fraction", "spike_label")) %>%
  dplyr::left_join(rank_summary_detail_study, by = c("tool", "background_study", "background_condition", "spike_fraction", "spike_label"))

summary_scores_detail_study <- summary_detail_wide_study %>%
  dplyr::mutate(
    tax_detect_score = score_higher_better(.data$tax_detect_raw),
    tax_recovery_score = score_ideal_one(.data$tax_recovery_raw),
    recovery_var_score = score_lower_better(.data$recovery_var_raw),
    fp_taxa_score = score_lower_better(.data$fp_taxa_raw),
    bio_detect_score = score_higher_better(.data$bio_detect_raw),
    q_strength_score = score_higher_better(.data$q_strength_raw),
    rank_score = score_lower_better(.data$rank_raw),
    spike_label = factor(.data$spike_label, levels = rev(row_order)),
    tool_pretty = dplyr::recode(as.character(.data$tool), !!!tool_labels, .default = as.character(.data$tool))
  )

baseline_long_detail_study <- summary_scores_detail_study %>%
  dplyr::transmute(
    tool = .data$tool,
    tool_pretty = .data$tool_pretty,
    background_study = .data$background_study,
    background_condition = .data$background_condition,
    spike_fraction = .data$spike_fraction,
    spike_label = .data$spike_label,
    metric = "Baseline
log10(x+1e-6)",
    fill_value = .data$baseline_log_raw,
    raw_label = format_label("baseline_log_raw", .data$baseline_log_raw)
  )

taxonomy_long_detail_study <- dplyr::bind_rows(
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Tax detect
(100% best)", fill_value = .data$tax_detect_score, raw_label = format_label("tax_detect_raw", .data$tax_detect_raw)),
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Tax recovery
(1 best)", fill_value = .data$tax_recovery_score, raw_label = format_label("tax_recovery_raw", .data$tax_recovery_raw)),
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Recovery variability
(0 best)", fill_value = .data$recovery_var_score, raw_label = format_label("recovery_var_raw", .data$recovery_var_raw)),
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Low FP taxa
(0 best)", fill_value = .data$fp_taxa_score, raw_label = format_label("fp_taxa_raw", .data$fp_taxa_raw))
)

biomarker_long_detail_study <- dplyr::bind_rows(
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Bio detect
(100% best)", fill_value = .data$bio_detect_score, raw_label = format_label("bio_detect_raw", .data$bio_detect_raw)),
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Bio -log10(q)
(high best)", fill_value = .data$q_strength_score, raw_label = format_label("q_strength_raw", .data$q_strength_raw)),
  summary_scores_detail_study %>% dplyr::transmute(tool = .data$tool, tool_pretty = .data$tool_pretty, background_study = .data$background_study, background_condition = .data$background_condition, spike_fraction = .data$spike_fraction, spike_label = .data$spike_label, metric = "Bio rank
(low best)", fill_value = .data$rank_score, raw_label = format_label("rank_raw", .data$rank_raw))
)

baseline_long_detail_study$metric <- factor(baseline_long_detail_study$metric, levels = "Baseline
log10(x+1e-6)")
taxonomy_long_detail_study$metric <- factor(taxonomy_long_detail_study$metric, levels = c("Tax detect
(100% best)", "Tax recovery
(1 best)", "Recovery variability
(0 best)", "Low FP taxa
(0 best)"))
biomarker_long_detail_study$metric <- factor(biomarker_long_detail_study$metric, levels = c("Bio detect
(100% best)", "Bio -log10(q)
(high best)", "Bio rank
(low best)"))

write_csv_safe(summary_scores_detail_study, file.path(opt$outdir, "species_driver_summary_by_study_condition_fraction_v6.csv"))
write_csv_safe(baseline_long_detail_study, file.path(opt$outdir, "species_driver_baseline_by_study_condition_fraction_v6.csv"))
write_csv_safe(taxonomy_long_detail_study, file.path(opt$outdir, "species_driver_taxonomy_by_study_condition_fraction_v6.csv"))
write_csv_safe(biomarker_long_detail_study, file.path(opt$outdir, "species_driver_biomarker_by_study_condition_fraction_v6.csv"))

short_study_label <- function(x) {
  x <- as.character(x)
  out <- sub("_.*$", "", x)
  out[!nzchar(out)] <- x[!nzchar(out)]
  out
}

study_order_detail <- summary_scores_detail_study %>%
  dplyr::distinct(.data$background_study) %>%
  dplyr::filter(!is.na(.data$background_study), nzchar(.data$background_study)) %>%
  dplyr::pull(.data$background_study)

study_label_map <- stats::setNames(short_study_label(study_order_detail), study_order_detail)

metric_study_levels <- function(metric_levels, study_order, study_label_map) {
  unlist(lapply(metric_levels, function(mm) {
    paste0(as.character(mm), "
", unname(study_label_map[study_order]))
  }), use.names = FALSE)
}

add_metric_study_column <- function(dat, metric_levels, study_order, study_label_map) {
  dat %>%
    dplyr::mutate(
      study_short = unname(study_label_map[as.character(.data$background_study)]),
      metric_study = factor(
        paste0(as.character(.data$metric), "
", .data$study_short),
        levels = metric_study_levels(metric_levels, study_order, study_label_map)
      )
    )
}

baseline_long_detail_study_split <- add_metric_study_column(
  baseline_long_detail_study,
  metric_levels = levels(baseline_long_detail_study$metric),
  study_order = study_order_detail,
  study_label_map = study_label_map
)

taxonomy_long_detail_study_split <- add_metric_study_column(
  taxonomy_long_detail_study,
  metric_levels = levels(taxonomy_long_detail_study$metric),
  study_order = study_order_detail,
  study_label_map = study_label_map
)

biomarker_long_detail_study_split <- add_metric_study_column(
  biomarker_long_detail_study,
  metric_levels = levels(biomarker_long_detail_study$metric),
  study_order = study_order_detail,
  study_label_map = study_label_map
)

make_tool_row_detail_split_study <- function(tool_name, condition_name, spike_frac, show_y = FALSE) {
  tpretty <- dplyr::recode(tool_name, !!!tool_labels, .default = tool_name)
  
  base_dat <- baseline_long_detail_study_split %>%
    dplyr::filter(.data$tool == tool_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  tax_dat <- taxonomy_long_detail_study_split %>%
    dplyr::filter(.data$tool == tool_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  bio_dat <- biomarker_long_detail_study_split %>%
    dplyr::filter(.data$tool == tool_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  if (!nrow(base_dat) && !nrow(tax_dat) && !nrow(bio_dat)) return(NULL)
  
  p_base <- base_dat %>%
    base_heat("metric_study", "raw_label", "fill_value", "Baseline context", show_y = show_y, palette_type = "baseline") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  
  p_tax <- tax_dat %>%
    base_heat("metric_study", "raw_label", "fill_value", "Taxonomy-side drivers", show_y = FALSE, palette_type = "favorability") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  
  p_bio <- bio_dat %>%
    base_heat("metric_study", "raw_label", "fill_value", "Biomarker outcome", show_y = FALSE, palette_type = "favorability") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  
  (p_base + p_tax + p_bio + patchwork::plot_layout(widths = c(1.8, 6.2, 4.5), guides = "collect")) +
    patchwork::plot_annotation(title = tpretty) &
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

make_detail_plot_split_study <- function(condition_name, spike_frac) {
  rows <- lapply(tool_order, function(tt) make_tool_row_detail_split_study(tt, condition_name, spike_frac, show_y = TRUE))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  
  wrap_plots(rows, ncol = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Which factors explain spike biomarker significance with both studies shown? (%s spikes)", opt$focus_mode),
      subtitle = paste(
        sprintf("Background condition = %s; spike fraction = %s.", condition_name, frac_to_label(spike_frac)),
        "Each metric is repeated once per study, so you can compare both studies within the same panel."
      ),
      caption = paste(
        "Baseline values are log10(relative abundance + 1e-6).",
        "Tax recovery is observed/expected, where 1 is ideal.",
        "Recovery variability is SD(log2(observed/expected)) across matched observations, where 0 is ideal.",
        "Low FP taxa is the absolute number of false-positive taxa.",
        "Colors in the taxonomy and biomarker blocks show closeness to the metric-specific ideal; numbers are raw values."
      )
    ) &
    theme_pub() &
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0, margin = margin(b = 10)),
      legend.position = "right"
    )
}

detail_keys_split_study <- summary_scores_detail_study %>%
  dplyr::distinct(.data$background_condition, .data$spike_fraction) %>%
  dplyr::filter(!is.na(.data$background_condition), nzchar(.data$background_condition), is.finite(.data$spike_fraction)) %>%
  dplyr::arrange(.data$background_condition, .data$spike_fraction)

if (nrow(detail_keys_split_study)) {
  detail_height_split <- max(10.5, 0.55 * length(row_order) + 6.2)
  pdf_path_split <- file.path(opt$outdir, "figure_species_driver_heatmap_by_condition_fraction_both_studies_v6.pdf")
  
  grDevices::pdf(pdf_path_split, width = 19, height = detail_height_split, onefile = TRUE)
  for (i in seq_len(nrow(detail_keys_split_study))) {
    cond_i <- detail_keys_split_study$background_condition[[i]]
    frac_i <- detail_keys_split_study$spike_fraction[[i]]
    
    p_i <- make_detail_plot_split_study(cond_i, frac_i)
    if (is.null(p_i)) next
    
    print(p_i)
    
    png_path <- file.path(
      opt$outdir,
      sprintf(
        "figure_species_driver_heatmap__cond-%s__spike-%s__both-studies_v6.png",
        slugify(cond_i),
        slugify(frac_to_label(frac_i))
      )
    )
    save_pub(p_i, png_path, width = 19, height = detail_height_split)
  }
  grDevices::dev.off()
}

make_tool_row_detail_study <- function(tool_name, study_name, condition_name, spike_frac, show_y = FALSE) {
  tpretty <- dplyr::recode(tool_name, !!!tool_labels, .default = tool_name)
  
  base_dat <- baseline_long_detail_study %>%
    dplyr::filter(.data$tool == tool_name, .data$background_study == study_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  tax_dat <- taxonomy_long_detail_study %>%
    dplyr::filter(.data$tool == tool_name, .data$background_study == study_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  bio_dat <- biomarker_long_detail_study %>%
    dplyr::filter(.data$tool == tool_name, .data$background_study == study_name, .data$background_condition == condition_name, dplyr::near(.data$spike_fraction, spike_frac))
  
  if (!nrow(base_dat) && !nrow(tax_dat) && !nrow(bio_dat)) return(NULL)
  
  p_base <- base_dat %>%
    base_heat("metric", "raw_label", "fill_value", "Baseline context", show_y = show_y, palette_type = "baseline")
  
  p_tax <- tax_dat %>%
    base_heat("metric", "raw_label", "fill_value", "Taxonomy-side drivers", show_y = FALSE, palette_type = "favorability") +
    theme(legend.position = "none")
  
  p_bio <- bio_dat %>%
    base_heat("metric", "raw_label", "fill_value", "Biomarker outcome", show_y = FALSE, palette_type = "favorability")
  
  (p_base + p_tax + p_bio + patchwork::plot_layout(widths = c(1.3, 4.4, 3.1), guides = "collect")) +
    patchwork::plot_annotation(title = tpretty) &
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

make_detail_plot_study <- function(study_name, condition_name, spike_frac) {
  rows <- lapply(tool_order, function(tt) make_tool_row_detail_study(tt, study_name, condition_name, spike_frac, show_y = TRUE))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  
  wrap_plots(rows, ncol = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Which factors explain spike biomarker significance by study? (%s spikes)", opt$focus_mode),
      subtitle = paste(
        sprintf("Study = %s; background condition = %s; spike fraction = %s.", study_name, condition_name, frac_to_label(spike_frac)),
        "This is the study-level companion to the LoD heatmaps."
      ),
      caption = paste(
        "Baseline values are log10(relative abundance + 1e-6).",
        "Tax recovery is observed/expected, where 1 is ideal.",
        "Recovery variability is SD(log2(observed/expected)) across matched observations, where 0 is ideal.",
        "Low FP taxa is the absolute number of false-positive taxa.",
        "Colors in the taxonomy and biomarker blocks show closeness to the metric-specific ideal; numbers are raw values."
      )
    ) &
    theme_pub() &
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0, margin = margin(b = 10)),
      legend.position = "right"
    )
}

detail_keys_study <- summary_scores_detail_study %>%
  dplyr::distinct(.data$background_study, .data$background_condition, .data$spike_fraction) %>%
  dplyr::filter(!is.na(.data$background_study), nzchar(.data$background_study), !is.na(.data$background_condition), nzchar(.data$background_condition), is.finite(.data$spike_fraction)) %>%
  dplyr::arrange(.data$background_study, .data$background_condition, .data$spike_fraction)

if (nrow(detail_keys_study)) {
  detail_height <- max(10.5, 0.55 * length(row_order) + 6.2)
  pdf_path <- file.path(opt$outdir, "figure_species_driver_heatmap_by_study_condition_fraction_v6.pdf")
  
  grDevices::pdf(pdf_path, width = 15.5, height = detail_height, onefile = TRUE)
  for (i in seq_len(nrow(detail_keys_study))) {
    study_i <- detail_keys_study$background_study[[i]]
    cond_i <- detail_keys_study$background_condition[[i]]
    frac_i <- detail_keys_study$spike_fraction[[i]]
    
    p_i <- make_detail_plot_study(study_i, cond_i, frac_i)
    if (is.null(p_i)) next
    
    print(p_i)
    
    png_path <- file.path(
      opt$outdir,
      sprintf(
        "figure_species_driver_heatmap__study-%s__cond-%s__spike-%s_v6.png",
        slugify(study_i),
        slugify(cond_i),
        slugify(frac_to_label(frac_i))
      )
    )
    save_pub(p_i, png_path, width = 15.5, height = detail_height)
  }
  grDevices::dev.off()
}

threshold_source <- focus_maaslin %>%
  dplyr::mutate(
    row_label = choose_row_label(., label_style = opt$threshold_label_style, focus_mode = opt$focus_mode),
    detected_use = as.numeric(dplyr::coalesce(.data[[detect_col]], FALSE))
  ) %>%
  dplyr::filter(!is.na(.data$row_label), nzchar(.data$row_label), is.finite(.data$spike_fraction))

relevance_by_fraction <- threshold_source %>%
  dplyr::group_by(.data$tool, .data$background_condition, .data$background_study, .data$row_label, .data$spike_fraction) %>%
  dplyr::summarise(
    relevant_any = as.numeric(any(.data$detected_use > 0, na.rm = TRUE)),
    detection_rate = mean(.data$detected_use, na.rm = TRUE),
    n_positive = sum(.data$detected_use > 0, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

threshold_overall <- relevance_by_fraction %>%
  dplyr::group_by(.data$tool, .data$background_condition, .data$row_label, .data$spike_fraction) %>%
  dplyr::summarise(
    relevant_any = as.numeric(any(.data$relevant_any > 0, na.rm = TRUE)),
    detection_rate = weighted.mean(.data$detection_rate, w = .data$n, na.rm = TRUE),
    n_positive = sum(.data$n_positive, na.rm = TRUE),
    n_total = sum(.data$n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$tool, .data$background_condition, .data$row_label, .data$spike_fraction) %>%
  dplyr::group_by(.data$tool, .data$background_condition, .data$row_label) %>%
  dplyr::summarise(
    threshold_fraction = {
      ok <- .data$spike_fraction[.data$relevant_any > 0]
      if (length(ok)) min(ok, na.rm = TRUE) else NA_real_
    },
    best_detection_rate = max(.data$detection_rate, na.rm = TRUE),
    n_positive = sum(.data$n_positive, na.rm = TRUE),
    n_total = sum(.data$n_total, na.rm = TRUE),
    .groups = "drop"
  )

threshold_by_study <- relevance_by_fraction %>%
  dplyr::arrange(.data$tool, .data$background_study, .data$background_condition, .data$row_label, .data$spike_fraction) %>%
  dplyr::group_by(.data$tool, .data$background_study, .data$background_condition, .data$row_label) %>%
  dplyr::summarise(
    threshold_fraction = {
      ok <- .data$spike_fraction[.data$relevant_any > 0]
      if (length(ok)) min(ok, na.rm = TRUE) else NA_real_
    },
    best_detection_rate = max(.data$detection_rate, na.rm = TRUE),
    n_positive = sum(.data$n_positive, na.rm = TRUE),
    n_total = sum(.data$n, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(relevance_by_fraction, file.path(opt$outdir, "detection_rate_by_fraction_v6.csv"))
write_csv_safe(threshold_overall, file.path(opt$outdir, "detection_threshold_overall_v6.csv"))
write_csv_safe(threshold_by_study, file.path(opt$outdir, "detection_threshold_by_study_v6.csv"))

plot_threshold_heatmap(
  threshold_overall,
  file.path(opt$outdir, "figure_detection_threshold_heatmap_v6.png"),
  sprintf("Minimum spike fraction where the biomarker was relevant (%s spikes)", opt$focus_mode),
  sprintf("Threshold = smallest spike fraction where at least one %s call was observed across all studies.",
          ifelse(detect_col == "member_detected_positive", "positive biomarker", "target biomarker")),
  row_order_pref = row_order,
  by_study = FALSE,
  caption_txt = paste(
    "Numbers show the first tested spike fraction where the biomarker was called relevant.",
    "NR = never relevant within the tested spike fractions.",
    "Darker blue means relevance appeared at lower spike fractions."
  )
)

if (any(!is.na(threshold_by_study$background_study) & nzchar(threshold_by_study$background_study))) {
  plot_threshold_heatmap(
    threshold_by_study,
    file.path(opt$outdir, "figure_detection_threshold_heatmap_by_study_v6.png"),
    sprintf("Minimum spike fraction where the biomarker was relevant by study (%s spikes)", opt$focus_mode),
    sprintf("Threshold = smallest spike fraction where at least one %s call was observed within each study.",
            ifelse(detect_col == "member_detected_positive", "positive biomarker", "target biomarker")),
    row_order_pref = row_order,
    by_study = TRUE,
    caption_txt = paste(
      "Numbers show the first tested spike fraction where the biomarker was called relevant within each study.",
      "NR = never relevant within the tested spike fractions.",
      "Darker blue means relevance appeared at lower spike fractions."
    )
  )
}

cat(sprintf("[OK] Wrote outputs under %s\n", opt$outdir))
