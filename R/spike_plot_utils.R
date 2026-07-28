sem <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
binom_sem <- function(p, n) sqrt(p * (1 - p) / pmax(1, n))

wilson_ci <- function(k, n, conf = 0.95) {
  if (!is.finite(k) || !is.finite(n) || n <= 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  denom <- 1 + (z^2 / n)
  center <- (p + (z^2 / (2 * n))) / denom
  half <- (z * sqrt((p * (1 - p) + (z^2 / (4 * n))) / n)) / denom
  c(lower = max(0, center - half), upper = min(1, center + half))
}

tool_palette <- function(vals = NULL) {
  pal <- c(
    metaphlan4 = "#1f77b4",
    metaphlan = "#1f77b4",
    kraken2_bracken = "#d62728",
    kraken_bracken = "#d62728",
    diamond_megan = "#2ca02c"
  )
  if (is.null(vals)) return(pal)
  extras <- setdiff(vals, names(pal))
  if (length(extras) > 0) {
    extra_cols <- grDevices::hcl.colors(length(extras), palette = "Dark 3")
    names(extra_cols) <- extras
    pal <- c(pal, extra_cols)
  }
  pal[vals]
}


condition_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}

blend_with_white <- function(col, amount = 0.3) {
  rgb0 <- grDevices::col2rgb(col) / 255
  out <- rgb0 * (1 - amount) + amount
  grDevices::rgb(out[1, ], out[2, ], out[3, ])
}

condition_palette <- function(vals) {
  if (length(vals) == 0) return(character())
  vals_chr <- as.character(vals)
  keys <- condition_key(vals_chr)
  base_map <- c(
    control = "#2C7FB8",
    healthycontrol = "#2C7FB8",
    normal = "#2C7FB8",
    adenoma = "#E69F00",
    advancedadenoma = "#D95F02",
    colorectalcarcinoma = "#C0392B",
    crc = "#C0392B",
    cancer = "#C0392B"
  )
  cols <- unname(base_map[keys])
  miss <- is.na(cols)
  if (any(miss)) {
    extra_cols <- grDevices::hcl.colors(sum(miss), palette = "Dark 3")
    cols[miss] <- extra_cols
  }
  names(cols) <- vals_chr
  cols
}

condition_study_palette <- function(condition_vals, study_vals) {
  if (length(condition_vals) == 0 || length(study_vals) == 0) return(character())
  df <- tibble::tibble(
    condition = as.character(condition_vals),
    study = as.character(study_vals)
  ) %>%
    dplyr::filter(!is.na(.data$condition), nzchar(.data$condition), !is.na(.data$study), nzchar(.data$study)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$condition, .data$study)
  if (!nrow(df)) return(character())
  base_cols <- condition_palette(unique(df$condition))
  pal <- character(nrow(df))
  for (cond in unique(df$condition)) {
    idx <- which(df$condition == cond)
    n <- length(idx)
    amounts <- if (n == 1) 0.18 else seq(0.08, 0.45, length.out = n)
    pal[idx] <- vapply(amounts, function(a) blend_with_white(base_cols[[cond]], a), character(1))
  }
  names(pal) <- paste(df$condition, df$study, sep = "||")
  pal
}

has_study_column <- function(df) {
  is.data.frame(df) && "Study" %in% names(df) && any(!is.na(df$Study) & nzchar(df$Study))
}


theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.25),
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", colour = "#D9D9D9"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5),
      plot.caption = ggplot2::element_text(size = base_size - 2, colour = "#555555"),
      axis.text.x = ggplot2::element_text(colour = "#222222"),
      axis.text.y = ggplot2::element_text(colour = "#222222")
    )
}

save_pub_plot <- function(plot_obj, path, width = 8, height = 5) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot_obj, width = width, height = height, units = "in", dpi = 320, bg = "white")
}

x_breaks_fun <- function(x) sort(unique(x[is.finite(x) & x > 0]))

x_axis_fraction <- function(xvals) {
  ggplot2::scale_x_continuous(
    trans = "log10",
    breaks = x_breaks_fun(xvals),
    labels = scales::label_number(scale = 100, accuracy = 0.01, suffix = "%")
  )
}

make_spike_factor <- function(df, col = "spike_fraction_total") {
  vals <- sort(unique(df[[col]][is.finite(df[[col]]) & df[[col]] > 0]))
  labs <- scales::label_number(scale = 100, accuracy = 0.01, suffix = "%")(vals)
  df[["spike_factor"]] <- factor(df[[col]], levels = vals, labels = labs)
  df
}

save_condition_curve <- function(dat, ycol, semcol, ylab, title, path) {
  if (!nrow(dat) || !"Target_Condition" %in% names(dat)) return(invisible(NULL))
  cols <- condition_palette(sort(unique(dat$Target_Condition)))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_fraction_total, y = .data[[ycol]], colour = .data$Target_Condition, fill = .data$Target_Condition)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[ycol]] - .data[[semcol]], ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.15, colour = NA) +
    x_axis_fraction(dat$spike_fraction_total) +
    ggplot2::facet_grid(tool ~ spike_label, scales = "free_y") +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = "Spike fraction", y = ylab, title = title, colour = "Condition", fill = "Condition") +
    theme_pub()
  save_pub_plot(p, path, width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 6.8)
}

plot_mean_sem_generic <- function(summary_df, y_mean, y_sem, ylab, fname, outdir, hline = NULL, percent_y = FALSE, caption = NULL) {
  if (!nrow(summary_df) || !all(c(y_mean, y_sem) %in% names(summary_df))) return(invisible(NULL))
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$spike_fraction_total, y = .data[[y_mean]], colour = .data$tool, fill = .data$tool)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[y_mean]] - .data[[y_sem]], ymax = .data[[y_mean]] + .data[[y_sem]]), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    x_axis_fraction(summary_df$spike_fraction_total) +
    ggplot2::facet_grid(tool ~ spike_label, scales = if (percent_y) "free_y" else "free_y") +
    ggplot2::scale_colour_manual(values = tool_palette(sort(unique(summary_df$tool))), drop = FALSE) +
    ggplot2::scale_fill_manual(values = tool_palette(sort(unique(summary_df$tool))), drop = FALSE) +
    ggplot2::labs(x = "Spike fraction", y = ylab, colour = "Tool", fill = "Tool", caption = caption) +
    theme_pub()
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, linewidth = 0.35)
  if (percent_y) p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  save_pub_plot(p, file.path(outdir, fname), width = max(10, 2.8 * dplyr::n_distinct(summary_df$spike_label)), height = 6.8)
}

plot_summary_curves <- function(combined_df, outdir, prefix = "") {
  if (!nrow(combined_df)) return(invisible(NULL))
  ensure_dir(outdir)
  dat <- combined_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_fraction_total) %>%
    dplyr::summarise(
      target_recovery_mean = mean(.data$target_recovery, na.rm = TRUE),
      target_recovery_sem = sem(.data$target_recovery),
      target_abs_error_mean = mean(.data$target_abs_error, na.rm = TRUE),
      target_abs_error_sem = sem(.data$target_abs_error),
      false_positive_mass_mean = mean(.data$false_positive_mass, na.rm = TRUE),
      false_positive_mass_sem = sem(.data$false_positive_mass),
      off_target_abs_error_mean = mean(.data$off_target_abs_error, na.rm = TRUE),
      off_target_abs_error_sem = sem(.data$off_target_abs_error),
      detection_fraction_mean = mean(.data$target_detection_fraction, na.rm = TRUE),
      detection_fraction_sem = sem(.data$target_detection_fraction),
      .groups = "drop"
    )
  write_csv_safe(dat, file.path(outdir, paste0(prefix, "per_spike_means_sem.csv")))

  plot_mean_sem_generic(dat, "target_recovery_mean", "target_recovery_sem", "Observed / expected target mass", paste0(prefix, "meansem_target_recovery.png"), outdir, hline = 1, caption = "Dashed line marks perfect recovery.")
  plot_mean_sem_generic(dat, "target_abs_error_mean", "target_abs_error_sem", "Absolute target error", paste0(prefix, "meansem_target_abs_error.png"), outdir, percent_y = TRUE)
  plot_mean_sem_generic(dat, "false_positive_mass_mean", "false_positive_mass_sem", "False-positive abundance mass", paste0(prefix, "meansem_false_positive_mass.png"), outdir, percent_y = TRUE)
  plot_mean_sem_generic(dat, "off_target_abs_error_mean", "off_target_abs_error_sem", "Off-target absolute error", paste0(prefix, "meansem_offtarget_abs_error.png"), outdir, percent_y = TRUE)
  plot_mean_sem_generic(dat, "detection_fraction_mean", "detection_fraction_sem", "Detected target-member fraction", paste0(prefix, "meansem_target_detection_fraction.png"), outdir, hline = 1, percent_y = TRUE)

  if ("Target_Condition" %in% names(combined_df)) {
    cond_dat <- combined_df %>%
      dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
      dplyr::group_by(.data$tool, .data$spike_label, .data$Target_Condition, .data$spike_fraction_total) %>%
      dplyr::summarise(
        target_recovery_mean = mean(.data$target_recovery, na.rm = TRUE),
        target_recovery_sem = sem(.data$target_recovery),
        false_positive_mass_mean = mean(.data$false_positive_mass, na.rm = TRUE),
        false_positive_mass_sem = sem(.data$false_positive_mass),
        .groups = "drop"
      )
    save_condition_curve(cond_dat, "target_recovery_mean", "target_recovery_sem", "Observed / expected target mass", "Target recovery by condition", file.path(outdir, paste0(prefix, "target_recovery_by_condition.png")))
    save_condition_curve(cond_dat, "false_positive_mass_mean", "false_positive_mass_sem", "False-positive abundance mass", "False-positive mass by condition", file.path(outdir, paste0(prefix, "false_positive_mass_by_condition.png")))
  }
  invisible(NULL)
}

plot_target_member_error_boxplots <- function(trace_df, outdir, prefix = "") {
  if (!nrow(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))

  p1 <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_factor, y = .data$relative_error, fill = .data$tool)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#666666") +
    ggplot2::geom_boxplot(outlier.alpha = 0.18) +
    ggplot2::facet_grid(tool ~ spike_label) +
    ggplot2::scale_fill_manual(values = tool_palette(sort(unique(dat$tool))), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Spike fraction", y = "Target-member relative error", title = "Target-member relative error across spike fractions") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
  save_pub_plot(p1, file.path(outdir, paste0(prefix, "target_member_error_boxplot.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 6.8)

  dat_zoom <- dat %>% dplyr::filter(.data$relative_error >= -0.25, .data$relative_error <= 0.25)
  if (nrow(dat_zoom)) {
    p2 <- ggplot2::ggplot(dat_zoom, ggplot2::aes(x = .data$spike_factor, y = .data$relative_error, fill = .data$tool)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#666666") +
      ggplot2::geom_boxplot(outlier.alpha = 0.18) +
      ggplot2::facet_grid(tool ~ spike_label) +
      ggplot2::scale_fill_manual(values = tool_palette(sort(unique(dat_zoom$tool))), drop = FALSE) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.25, 0.25)) +
      ggplot2::labs(x = "Spike fraction", y = "Target-member relative error", title = "Target-member relative error (zoomed to ±25%)") +
      theme_pub() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
    save_pub_plot(p2, file.path(outdir, paste0(prefix, "target_member_error_boxplot_zoom.png")), width = max(10, 2.8 * dplyr::n_distinct(dat_zoom$spike_label)), height = 6.8)
  }
  invisible(NULL)
}

plot_lod_curve <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_fraction_total) %>%
    dplyr::summarise(
      n = sum(!is.na(.data$detected)),
      p_detect = mean(.data$detected, na.rm = TRUE),
      sem = binom_sem(mean(.data$detected, na.rm = TRUE), sum(!is.na(.data$detected))),
      .groups = "drop"
    )
  write_csv_safe(dat, file.path(outdir, paste0(prefix, "lod_detection_rates.csv")))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_fraction_total, y = .data$p_detect, colour = .data$tool, fill = .data$tool)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$p_detect - .data$sem), ymax = pmin(1, .data$p_detect + .data$sem)), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    x_axis_fraction(dat$spike_fraction_total) +
    ggplot2::facet_grid(tool ~ spike_label) +
    ggplot2::scale_colour_manual(values = tool_palette(sort(unique(dat$tool))), drop = FALSE) +
    ggplot2::scale_fill_manual(values = tool_palette(sort(unique(dat$tool))), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Spike fraction", y = "Detection probability", title = "Target-member detection probability", caption = "Detection probability is averaged across samples and target members.") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "lod_detection_probability.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 6.8)
  invisible(NULL)
}

plot_member_recovery_heatmap <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$member_taxon), nzchar(.data$member_taxon), is.finite(.data$observed_over_expected)) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$member_taxon, .data$spike_fraction_total) %>%
    dplyr::summarise(median_recovery = stats::median(.data$observed_over_expected, na.rm = TRUE), .groups = "drop") %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))
  dat <- dat %>% dplyr::mutate(median_recovery = pmax(0, pmin(.data$median_recovery, 2)))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_factor, y = .data$member_taxon, fill = .data$median_recovery)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::facet_grid(tool ~ spike_label, scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 1, limits = c(0, 2), oob = scales::squish) +
    ggplot2::labs(x = "Spike fraction", y = NULL, fill = "Observed / expected", title = "Target-member recovery heatmap") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "target_member_recovery_heatmap.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = max(6, 0.35 * dplyr::n_distinct(dat$member_taxon) + 4))
  invisible(NULL)
}

plot_target_member_outlier_partitions <- function(trace_df, outdir, prefix = "", outlier_thr = 0.10, var_thr = 0.10) {
  if (!nrow(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$baseline), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))

  species_class <- dat %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$member_taxon) %>%
    dplyr::summarise(
      median_baseline = stats::median(.data$baseline, na.rm = TRUE),
      max_abs_rel_error = max(abs(.data$relative_error), na.rm = TRUE),
      sd_rel_error = stats::sd(.data$relative_error, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      abund_bin = cut(100 * .data$median_baseline,
                      breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1, 100),
                      labels = c("0–0.01%", "0.01–0.05%", "0.05–0.1%", "0.1–0.5%", "0.5–1%", "1–100%"),
                      include.lowest = TRUE),
      category = dplyr::case_when(
        .data$max_abs_rel_error > outlier_thr ~ "Strong outlier",
        .data$sd_rel_error > var_thr ~ "High variance",
        TRUE ~ "Stable"
      )
    )
  if (!nrow(species_class)) return(invisible(NULL))

  summ <- species_class %>%
    dplyr::count(.data$tool, .data$spike_label, .data$spike_factor, .data$abund_bin, .data$category, name = "n") %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$abund_bin) %>%
    dplyr::mutate(p = .data$n / sum(.data$n)) %>%
    dplyr::ungroup()
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "outlier_partition_counts_by_spike_abundbin.csv")))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$abund_bin, y = .data$p, fill = .data$category)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Baseline abundance bin", y = "Percent of target members", fill = "Category", title = "Target-member stability across baseline-abundance bins") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "outlier_partition_stacked_relative.png")), width = max(11, 2.8 * dplyr::n_distinct(summ$spike_label)), height = 7.4)
  invisible(NULL)
}

plot_error_variance_scatter <- function(trace_df, outdir, prefix = "") {
  if (!nrow(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0) %>%
    make_spike_factor() %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$member_taxon) %>%
    dplyr::summarise(
      mean_abs_rel_error = mean(abs(.data$relative_error), na.rm = TRUE),
      sd_rel_error = stats::sd(.data$relative_error, na.rm = TRUE),
      .groups = "drop"
    )
  if (!nrow(dat)) return(invisible(NULL))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$mean_abs_rel_error, y = .data$sd_rel_error, colour = .data$tool)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.8) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_colour_manual(values = tool_palette(sort(unique(dat$tool))), drop = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Mean absolute relative error", y = "SD of relative error", colour = "Tool", title = "Target-member error magnitude versus variability") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "error_vs_variance_scatter.png")), width = max(11, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.4)

  dat_zoom <- dat %>% dplyr::filter(.data$mean_abs_rel_error <= 0.25, .data$sd_rel_error <= 0.25)
  if (nrow(dat_zoom)) {
    pz <- ggplot2::ggplot(dat_zoom, ggplot2::aes(x = .data$mean_abs_rel_error, y = .data$sd_rel_error, colour = .data$tool)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.9) +
      ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
      ggplot2::scale_colour_manual(values = tool_palette(sort(unique(dat_zoom$tool))), drop = FALSE) +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.25)) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.25)) +
      ggplot2::labs(x = "Mean absolute relative error", y = "SD of relative error", colour = "Tool", title = "Target-member error magnitude versus variability (zoomed to 25%)") +
      theme_pub()
    save_pub_plot(pz, file.path(outdir, paste0(prefix, "error_vs_variance_scatter_zoom.png")), width = max(11, 2.8 * dplyr::n_distinct(dat_zoom$spike_label)), height = 7.4)
  }
  invisible(NULL)
}

plot_target_dilution_checks <- function(trace_df, outdir, prefix = "", max_taxa_per_label = 12) {
  if (!nrow(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$member_taxon), nzchar(.data$member_taxon))
  if (!nrow(dat)) return(invisible(NULL))

  top_taxa <- dat %>%
    dplyr::group_by(.data$spike_label, .data$member_taxon) %>%
    dplyr::summarise(rank_mass = max(.data$expected, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$spike_label) %>%
    dplyr::slice_max(order_by = .data$rank_mass, n = max_taxa_per_label, with_ties = FALSE) %>%
    dplyr::ungroup()

  keep <- dat %>%
    dplyr::inner_join(top_taxa, by = c("spike_label", "member_taxon")) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$member_taxon, .data$spike_fraction_total) %>%
    dplyr::summarise(
      expected_mean = mean(.data$expected, na.rm = TRUE),
      observed_mean = mean(.data$observed, na.rm = TRUE),
      observed_sem = sem(.data$observed),
      .groups = "drop"
    )
  if (!nrow(keep)) return(invisible(NULL))

  for (lab in unique(keep$spike_label)) {
    dlab <- keep %>% dplyr::filter(.data$spike_label == !!lab)
    p <- ggplot2::ggplot(dlab, ggplot2::aes(x = .data$spike_fraction_total, y = .data$observed_mean, colour = .data$tool, fill = .data$tool)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$expected_mean), linetype = 2, linewidth = 0.8) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$observed_mean - .data$observed_sem, ymax = .data$observed_mean + .data$observed_sem), alpha = 0.15, colour = NA) +
      x_axis_fraction(dlab$spike_fraction_total) +
      ggplot2::facet_wrap(~ member_taxon, scales = "free_y", ncol = 3) +
      ggplot2::scale_colour_manual(values = tool_palette(sort(unique(dlab$tool))), drop = FALSE) +
      ggplot2::scale_fill_manual(values = tool_palette(sort(unique(dlab$tool))), drop = FALSE) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      ggplot2::labs(x = "Spike fraction", y = "Observed member abundance", title = paste0("Target-member dilution checks: ", lab), subtitle = "Dashed lines show expected abundance under perfect dilution.") +
      theme_pub()
    save_pub_plot(p, file.path(outdir, paste0(prefix, "dilution_", sanitize_slug(lab), ".png")), width = 12, height = max(6, 2.3 * ceiling(dplyr::n_distinct(dlab$member_taxon) / 3)))
  }
  invisible(NULL)
}

plot_condition_relerr_violin <- function(trace_df, outdir, prefix = "") {
  if (!nrow(trace_df) || !"Target_Condition" %in% names(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
    make_spike_factor() %>%
    dplyr::mutate(rel_err_plot = log10(abs(.data$relative_error) + 1e-3))
  if (!nrow(dat)) return(invisible(NULL))
  cols <- condition_palette(sort(unique(dat$Target_Condition)))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Target_Condition, y = .data$rel_err_plot, fill = .data$Target_Condition)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.18, outlier.alpha = 0.18) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = "Condition", y = "log10(|relative error| + 1e-3)", title = "Condition-specific target-member error") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
  save_pub_plot(p, file.path(outdir, paste0(prefix, "condition_rel_error_violin_log.png")), width = max(12, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.8)
  invisible(NULL)
}

plot_condition_lod_passrate <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df) || !"Target_Condition" %in% names(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))
  cols <- condition_palette(sort(unique(dat$Target_Condition)))

  summ <- dat %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$Target_Condition) %>%
    dplyr::summarise(
      n = sum(!is.na(.data$detected)),
      k = sum(.data$detected, na.rm = TRUE),
      p_pass = ifelse(.data$n[[1]] > 0, .data$k[[1]] / .data$n[[1]], NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ci = list(wilson_ci(.data$k, .data$n)),
      ci_low = .data$ci[[1]][1],
      ci_high = .data$ci[[1]][2]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$ci)
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "lod_passrate_by_condition_wilson95.csv")))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$Target_Condition, y = .data$p_pass, fill = .data$Target_Condition)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high), width = 0.2) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Condition", y = "Detection pass rate", title = "Target-member detection pass rate by condition") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
  save_pub_plot(p, file.path(outdir, paste0(prefix, "lod_passrate_by_condition.png")), width = max(12, 2.8 * dplyr::n_distinct(summ$spike_label)), height = 7.8)
  invisible(NULL)
}

plot_condition_error_bins <- function(trace_df, outdir, prefix = "") {
  if (!nrow(trace_df) || !"Target_Condition" %in% names(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
    make_spike_factor() %>%
    dplyr::mutate(
      err_bin = dplyr::case_when(
        abs(.data$relative_error) <= 0.10 ~ "Good (≤10%)",
        abs(.data$relative_error) <= 0.50 ~ "Borderline (10–50%)",
        TRUE ~ "Bad (>50%)"
      )
    )
  if (!nrow(dat)) return(invisible(NULL))

  summ <- dat %>%
    dplyr::count(.data$tool, .data$spike_label, .data$spike_factor, .data$Target_Condition, .data$err_bin, name = "n") %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$Target_Condition) %>%
    dplyr::mutate(p = .data$n / sum(.data$n)) %>%
    dplyr::ungroup()
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "relerr_bins_by_condition.csv")))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$Target_Condition, y = .data$p, fill = .data$err_bin)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Condition", y = "Percent of target-member observations", fill = "Error bin", title = "Target-member error bins by condition") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "relerr_bins_stacked_by_condition.png")), width = max(12, 2.8 * dplyr::n_distinct(summ$spike_label)), height = 7.8)
  invisible(NULL)
}

plot_baseline_target_mass_by_condition <- function(baseline_df, outdir, prefix = "") {
  if (!nrow(baseline_df) || !"Target_Condition" %in% names(baseline_df)) return(invisible(NULL))
  dat <- baseline_df %>%
    dplyr::filter(is.finite(.data$target_baseline_mass), !is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
    dplyr::mutate(y = log10(.data$target_baseline_mass + 1e-6))
  if (!nrow(dat)) return(invisible(NULL))
  cols <- condition_palette(sort(unique(dat$Target_Condition)))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Target_Condition, y = .data$y, fill = .data$Target_Condition)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.18, outlier.alpha = 0.2) +
    ggplot2::facet_grid(tool ~ spike_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(x = "Condition", y = "log10(target baseline mass + 1e-6)", title = "Baseline target abundance by condition") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
  save_pub_plot(p, file.path(outdir, paste0(prefix, "baseline_target_mass_by_condition_log.png")), width = max(11, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 6.8)
  invisible(NULL)
}

plot_bias_variance_summary <- function(combined_df, outdir, prefix = "") {
  if (!nrow(combined_df) || !"Target_Condition" %in% names(combined_df)) return(invisible(NULL))
  dat <- combined_df %>%
    dplyr::filter(
      is.finite(.data$target_recovery),
      is.finite(.data$spike_fraction_total),
      .data$spike_fraction_total > 0,
      !is.na(.data$Target_Condition),
      nzchar(.data$Target_Condition)
    ) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))
  cols <- condition_palette(sort(unique(dat$Target_Condition)))
  
  summ <- dat %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$Target_Condition) %>%
    dplyr::summarise(
      n = dplyr::n(),
      bias_median = stats::median(.data$target_recovery, na.rm = TRUE),
      spread_iqr = stats::IQR(.data$target_recovery, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "bias_variance_by_condition.csv")))
  
  p <- ggplot2::ggplot(
    summ,
    ggplot2::aes(
      x = .data$bias_median,
      y = .data$spread_iqr,
      colour = .data$Target_Condition,
      size = .data$n
    )
  ) +
    ggplot2::geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor, scales = "free") +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::labs(
      x = "Bias = median(target recovery)",
      y = "Spread = IQR(target recovery)",
      colour = "Condition",
      size = "n",
      title = "Bias–spread summary by condition",
      subtitle = "Each point summarizes one tool × spike label × spike fraction × condition group.",
      caption = paste(
        "Bias is the median of target_recovery across samples in the group.",
        "Spread is the IQR of target_recovery across samples in the group.",
        "For community spikes, this is based on the community-level target_recovery already stored in combined_df,",
        "not on summing per-member variances.",
        sep = " "
      )
    ) +
    theme_pub()
  
  save_pub_plot(
    p,
    file.path(outdir, paste0(prefix, "bias_variance_scatter_by_condition.png")),
    width = max(12, 2.8 * dplyr::n_distinct(summ$spike_label)),
    height = 7.8
  )
  invisible(NULL)
}

plot_target_recovery_vs_baseline <- function(combined_df, baseline_df, outdir, prefix = "") {
  if (!nrow(combined_df) || !nrow(baseline_df) || !"Target_Condition" %in% names(combined_df)) return(invisible(NULL))
  dat <- combined_df %>%
    dplyr::select(.data$sample_id, .data$tool, .data$spike_label, .data$spike_fraction_total, .data$Target_Condition, .data$target_recovery) %>%
    dplyr::left_join(
      baseline_df %>% dplyr::select(.data$sample_id, .data$tool, .data$spike_label, .data$target_baseline_mass) %>% dplyr::distinct(),
      by = c("sample_id", "tool", "spike_label")
    ) %>%
    dplyr::filter(is.finite(.data$target_recovery), is.finite(.data$target_baseline_mass), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$Target_Condition), nzchar(.data$Target_Condition)) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))
  cols <- condition_palette(sort(unique(dat$Target_Condition)))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$target_baseline_mass, y = .data$target_recovery, colour = .data$Target_Condition)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35) +
    ggplot2::geom_point(alpha = 0.7, size = 1.7) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor, scales = "free") +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), labels = scales::percent_format(accuracy = 0.01)) +
    ggplot2::labs(x = "Baseline target abundance", y = "Observed / expected target mass", colour = "Condition", title = "Target recovery versus baseline target abundance") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "target_recovery_vs_baseline.png")), width = max(12, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.8)
  invisible(NULL)
}

plot_community_member_detection_curves <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$member_taxon), nzchar(.data$member_taxon)) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$member_taxon, .data$spike_fraction_total) %>%
    dplyr::summarise(
      n = sum(!is.na(.data$detected)),
      detection_rate = mean(.data$detected, na.rm = TRUE),
      sem = binom_sem(mean(.data$detected, na.rm = TRUE), sum(!is.na(.data$detected))),
      .groups = "drop"
    )
  if (!nrow(dat)) return(invisible(NULL))
  member_cols <- condition_palette(sort(unique(dat$member_taxon)))
  write_csv_safe(dat, file.path(outdir, paste0(prefix, "community_member_detection_summary.csv")))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_fraction_total, y = .data$detection_rate, colour = .data$member_taxon, fill = .data$member_taxon)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$detection_rate - .data$sem), ymax = pmin(1, .data$detection_rate + .data$sem)), alpha = 0.12, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.9) +
    x_axis_fraction(dat$spike_fraction_total) +
    ggplot2::facet_grid(tool ~ spike_label) +
    ggplot2::scale_colour_manual(values = member_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = member_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Total community spike fraction", y = "Detection probability", colour = "Community member", fill = "Community member", title = "Community-member detection probability across total spike fractions", caption = "Each line tracks a community member separately within the multi-species spike.") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "community_member_detection_probability.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.0)
  invisible(NULL)
}

plot_community_member_recovery_curves <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$member_taxon), nzchar(.data$member_taxon), is.finite(.data$observed_over_expected)) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$member_taxon, .data$spike_fraction_total) %>%
    dplyr::summarise(
      mean_recovery = mean(.data$observed_over_expected, na.rm = TRUE),
      sem = sem(.data$observed_over_expected),
      .groups = "drop"
    )
  if (!nrow(dat)) return(invisible(NULL))
  member_cols <- condition_palette(sort(unique(dat$member_taxon)))
  write_csv_safe(dat, file.path(outdir, paste0(prefix, "community_member_recovery_summary.csv")))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_fraction_total, y = .data$mean_recovery, colour = .data$member_taxon, fill = .data$member_taxon)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35, colour = "#666666") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$mean_recovery - .data$sem), ymax = .data$mean_recovery + .data$sem), alpha = 0.12, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.9) +
    x_axis_fraction(dat$spike_fraction_total) +
    ggplot2::facet_grid(tool ~ spike_label) +
    ggplot2::scale_colour_manual(values = member_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = member_cols, drop = FALSE) +
    ggplot2::labs(x = "Total community spike fraction", y = "Observed / expected abundance", colour = "Community member", fill = "Community member", title = "Community-member recovery across total spike fractions", caption = "A value of 1 indicates perfect recovery for that member inside the spiked community.") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "community_member_recovery_curve.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.0)
  invisible(NULL)
}

plot_community_member_error_boxplots <- function(trace_df, outdir, prefix = "") {
  if (!nrow(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$member_taxon), nzchar(.data$member_taxon)) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))
  member_cols <- condition_palette(sort(unique(dat$member_taxon)))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_factor, y = .data$relative_error, fill = .data$member_taxon)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#666666") +
    ggplot2::geom_boxplot(outlier.alpha = 0.18) +
    ggplot2::facet_grid(tool ~ spike_label) +
    ggplot2::scale_fill_manual(values = member_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Total community spike fraction", y = "Relative error", fill = "Community member", title = "Community-member relative error distributions", caption = "Relative error is computed for each spiked community member separately: (observed − expected) / expected.") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "community_member_error_boxplot.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.0)

  dat_zoom <- dat %>% dplyr::filter(.data$relative_error >= -0.25, .data$relative_error <= 0.25)
  if (nrow(dat_zoom)) {
    pz <- ggplot2::ggplot(dat_zoom, ggplot2::aes(x = .data$spike_factor, y = .data$relative_error, fill = .data$member_taxon)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#666666") +
      ggplot2::geom_boxplot(outlier.alpha = 0.18) +
      ggplot2::facet_grid(tool ~ spike_label) +
      ggplot2::scale_fill_manual(values = member_cols, drop = FALSE) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.25, 0.25)) +
      ggplot2::labs(x = "Total community spike fraction", y = "Relative error", fill = "Community member", title = "Community-member relative error distributions (zoomed to ±25%)") +
      theme_pub() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_pub_plot(pz, file.path(outdir, paste0(prefix, "community_member_error_boxplot_zoom.png")), width = max(10, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.0)
  }
  invisible(NULL)
}

plot_community_member_condition_detection <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df) || !"Target_Condition" %in% names(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0, !is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$member_taxon), nzchar(.data$member_taxon)) %>%
    make_spike_factor() %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$member_taxon, .data$Target_Condition) %>%
    dplyr::summarise(detection_rate = mean(.data$detected, na.rm = TRUE), .groups = "drop")
  if (!nrow(dat)) return(invisible(NULL))
  write_csv_safe(dat, file.path(outdir, paste0(prefix, "community_member_detection_by_condition.csv")))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$spike_factor, y = .data$member_taxon, fill = .data$detection_rate)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::facet_grid(tool + Target_Condition ~ spike_label, scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#084594", limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Total community spike fraction", y = NULL, fill = "Detection", title = "Community-member detection by condition", caption = "Heatmap values show the fraction of samples in which each spiked community member was detected.") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "community_member_detection_by_condition_heatmap.png")), width = max(11, 2.8 * dplyr::n_distinct(dat$spike_label)), height = max(7, 0.32 * dplyr::n_distinct(dat$member_taxon) + 4.5))
  invisible(NULL)
}


plot_summary_curves_by_condition_and_study <- function(combined_df, outdir, prefix = "") {
  if (!nrow(combined_df) || !has_study_column(combined_df) || !"Target_Condition" %in% names(combined_df)) return(invisible(NULL))
  dat <- combined_df %>%
    dplyr::filter(
      is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0,
      !is.na(.data$Target_Condition), nzchar(.data$Target_Condition),
      !is.na(.data$Study), nzchar(.data$Study)
    ) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$Target_Condition, .data$Study, .data$spike_fraction_total) %>%
    dplyr::summarise(
      target_recovery_mean = mean(.data$target_recovery, na.rm = TRUE),
      target_recovery_sem = sem(.data$target_recovery),
      false_positive_mass_mean = mean(.data$false_positive_mass, na.rm = TRUE),
      false_positive_mass_sem = sem(.data$false_positive_mass),
      .groups = "drop"
    ) %>%
    dplyr::mutate(cond_study = paste(.data$Target_Condition, .data$Study, sep = "||"))
  if (!nrow(dat)) return(invisible(NULL))
  pal <- condition_study_palette(dat$Target_Condition, dat$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))

  make_curve <- function(ycol, semcol, ylab, title, fname, hline = NULL, percent_y = FALSE) {
    p <- ggplot2::ggplot(dat, ggplot2::aes(
      x = .data$spike_fraction_total,
      y = .data[[ycol]],
      colour = .data$cond_study,
      fill = .data$cond_study,
      group = interaction(.data$Target_Condition, .data$Study, drop = TRUE),
      linetype = .data$Study
    )) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[ycol]] - .data[[semcol]], ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.12, colour = NA) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 1.8) +
      x_axis_fraction(dat$spike_fraction_total) +
      ggplot2::facet_grid(tool ~ spike_label, scales = "free_y") +
      ggplot2::scale_colour_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
      ggplot2::scale_fill_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
      ggplot2::labs(x = "Spike fraction", y = ylab, colour = "Condition • Study", fill = "Condition • Study", linetype = "Study", title = title) +
      theme_pub()
    if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, linewidth = 0.35)
    if (percent_y) p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1))
    save_pub_plot(p, file.path(outdir, paste0(prefix, fname)), width = max(11, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.0)
  }

  make_curve("target_recovery_mean", "target_recovery_sem", "Observed / expected target mass", "Target recovery by condition and study", "target_recovery_by_condition_and_study.png", hline = 1)
  make_curve("false_positive_mass_mean", "false_positive_mass_sem", "False-positive abundance mass", "False-positive mass by condition and study", "false_positive_mass_by_condition_and_study.png", percent_y = TRUE)
  invisible(NULL)
}

plot_condition_relerr_violin_by_study <- function(trace_df, outdir, prefix = "") {
  if (!nrow(trace_df) || !has_study_column(trace_df) || !"Target_Condition" %in% names(trace_df)) return(invisible(NULL))
  dat <- trace_df %>%
    dplyr::filter(.data$target_flag, is.finite(.data$relative_error), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0,
                  !is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$Study), nzchar(.data$Study)) %>%
    make_spike_factor() %>%
    dplyr::mutate(
      rel_err_plot = log10(abs(.data$relative_error) + 1e-3),
      cond_study = paste(.data$Target_Condition, .data$Study, sep = "||")
    )
  if (!nrow(dat)) return(invisible(NULL))
  pal <- condition_study_palette(dat$Target_Condition, dat$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))
  dodge <- ggplot2::position_dodge(width = 0.82)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Target_Condition, y = .data$rel_err_plot, fill = .data$cond_study, group = .data$cond_study)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.85, position = dodge, linewidth = 0.25) +
    ggplot2::geom_boxplot(width = 0.16, outlier.alpha = 0.15, position = dodge) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_fill_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
    ggplot2::labs(x = "Condition", y = "log10(|relative error| + 1e-3)", fill = "Condition • Study", title = "Condition-specific target-member error by study") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "condition_rel_error_violin_log_by_study.png")), width = max(12, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 8.1)
  invisible(NULL)
}

plot_condition_lod_passrate_by_study <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df) || !has_study_column(target_df) || !"Target_Condition" %in% names(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0,
                  !is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$Study), nzchar(.data$Study)) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))

  summ <- dat %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$Target_Condition, .data$Study) %>%
    dplyr::summarise(
      n = sum(!is.na(.data$detected)),
      k = sum(.data$detected, na.rm = TRUE),
      p_pass = ifelse(.data$n[[1]] > 0, .data$k[[1]] / .data$n[[1]], NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ci = list(wilson_ci(.data$k, .data$n)),
      ci_low = .data$ci[[1]][1],
      ci_high = .data$ci[[1]][2],
      cond_study = paste(.data$Target_Condition, .data$Study, sep = "||")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$ci)
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "lod_passrate_by_condition_and_study_wilson95.csv")))

  pal <- condition_study_palette(summ$Target_Condition, summ$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))
  dodge <- ggplot2::position_dodge(width = 0.82)

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$Target_Condition, y = .data$p_pass, fill = .data$cond_study, group = .data$cond_study)) +
    ggplot2::geom_col(position = dodge) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high), width = 0.18, position = dodge) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor) +
    ggplot2::scale_fill_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = "Condition", y = "Detection pass rate", fill = "Condition • Study", title = "Target-member detection pass rate by condition and study") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "lod_passrate_by_condition_and_study.png")), width = max(12, 2.8 * dplyr::n_distinct(summ$spike_label)), height = 8.0)
  invisible(NULL)
}

plot_baseline_target_mass_by_condition_and_study <- function(baseline_df, outdir, prefix = "") {
  if (!nrow(baseline_df) || !has_study_column(baseline_df) || !"Target_Condition" %in% names(baseline_df)) return(invisible(NULL))
  dat <- baseline_df %>%
    dplyr::filter(is.finite(.data$target_baseline_mass), !is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$Study), nzchar(.data$Study)) %>%
    dplyr::mutate(y = log10(.data$target_baseline_mass + 1e-6), cond_study = paste(.data$Target_Condition, .data$Study, sep = "||"))
  if (!nrow(dat)) return(invisible(NULL))
  pal <- condition_study_palette(dat$Target_Condition, dat$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))
  dodge <- ggplot2::position_dodge(width = 0.82)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Target_Condition, y = .data$y, fill = .data$cond_study, group = .data$cond_study)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.85, position = dodge, linewidth = 0.25) +
    ggplot2::geom_boxplot(width = 0.16, outlier.alpha = 0.15, position = dodge) +
    ggplot2::facet_grid(tool ~ spike_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
    ggplot2::labs(x = "Condition", y = "log10(target baseline mass + 1e-6)", fill = "Condition • Study", title = "Baseline target abundance by condition and study") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "baseline_target_mass_by_condition_and_study_log.png")), width = max(11, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 7.3)
  invisible(NULL)
}

plot_bias_variance_summary_by_condition_and_study <- function(combined_df, outdir, prefix = "") {
  if (!nrow(combined_df) || !has_study_column(combined_df) || !"Target_Condition" %in% names(combined_df)) return(invisible(NULL))
  dat <- combined_df %>%
    dplyr::filter(is.finite(.data$target_recovery), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0,
                  !is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$Study), nzchar(.data$Study)) %>%
    make_spike_factor()
  if (!nrow(dat)) return(invisible(NULL))
  summ <- dat %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$spike_factor, .data$Target_Condition, .data$Study) %>%
    dplyr::summarise(n = dplyr::n(), bias_median = stats::median(.data$target_recovery, na.rm = TRUE), spread_iqr = stats::IQR(.data$target_recovery, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(cond_study = paste(.data$Target_Condition, .data$Study, sep = "||"))
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "bias_variance_by_condition_and_study.csv")))
  pal <- condition_study_palette(summ$Target_Condition, summ$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$bias_median, y = .data$spread_iqr, colour = .data$cond_study, shape = .data$Study, size = .data$n)) +
    ggplot2::geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor, scales = "free") +
    ggplot2::scale_colour_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
    ggplot2::labs(x = "Bias = median(target recovery)", y = "Spread = IQR(target recovery)", colour = "Condition • Study", shape = "Study", size = "n", title = "Bias–spread summary by condition and study") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "bias_variance_scatter_by_condition_and_study.png")), width = max(12, 2.8 * dplyr::n_distinct(summ$spike_label)), height = 8.0)
  invisible(NULL)
}

plot_target_recovery_vs_baseline_by_condition_and_study <- function(combined_df, baseline_df, outdir, prefix = "") {
  if (!nrow(combined_df) || !nrow(baseline_df) || !has_study_column(combined_df) || !"Target_Condition" %in% names(combined_df)) return(invisible(NULL))
  dat <- combined_df %>%
    dplyr::select(.data$sample_id, .data$tool, .data$spike_label, .data$spike_fraction_total, .data$Target_Condition, .data$Study, .data$target_recovery) %>%
    dplyr::left_join(
      baseline_df %>% dplyr::select(.data$sample_id, .data$tool, .data$spike_label, .data$target_baseline_mass) %>% dplyr::distinct(),
      by = c("sample_id", "tool", "spike_label")
    ) %>%
    dplyr::filter(is.finite(.data$target_recovery), is.finite(.data$target_baseline_mass), is.finite(.data$spike_fraction_total), .data$spike_fraction_total > 0,
                  !is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$Study), nzchar(.data$Study)) %>%
    make_spike_factor() %>%
    dplyr::mutate(cond_study = paste(.data$Target_Condition, .data$Study, sep = "||"))
  if (!nrow(dat)) return(invisible(NULL))
  pal <- condition_study_palette(dat$Target_Condition, dat$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$target_baseline_mass, y = .data$target_recovery, colour = .data$cond_study, shape = .data$Study)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35) +
    ggplot2::geom_point(alpha = 0.7, size = 1.8) +
    ggplot2::facet_grid(tool + spike_label ~ spike_factor, scales = "free") +
    ggplot2::scale_colour_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), labels = scales::percent_format(accuracy = 0.01)) +
    ggplot2::labs(x = "Baseline target abundance", y = "Observed / expected target mass", colour = "Condition • Study", shape = "Study", title = "Target recovery versus baseline abundance by condition and study") +
    theme_pub()
  save_pub_plot(p, file.path(outdir, paste0(prefix, "target_recovery_vs_baseline_by_condition_and_study.png")), width = max(12, 2.8 * dplyr::n_distinct(dat$spike_label)), height = 8.0)
  invisible(NULL)
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

plot_community_member_baseline_by_condition <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df) || !"Target_Condition" %in% names(target_df)) return(invisible(NULL))
  
  dat <- target_df %>%
    dplyr::filter(
      !is.na(.data$Target_Condition), nzchar(.data$Target_Condition),
      !is.na(.data$member_taxon), nzchar(.data$member_taxon),
      is.finite(.data$baseline)
    ) %>%
    dplyr::mutate(
      member_short = short_member_label(.data$member_taxon),
      baseline_log = log10(.data$baseline + 1e-6)
    )
  if (!nrow(dat)) return(invisible(NULL))
  
  member_levels <- dat %>%
    dplyr::group_by(.data$member_taxon) %>%
    dplyr::summarise(
      member_short = dplyr::first(.data$member_short),
      baseline_median = stats::median(.data$baseline, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$baseline_median), .data$member_taxon)
  
  dat <- dat %>%
    dplyr::mutate(
      member_taxon = factor(.data$member_taxon, levels = member_levels$member_taxon, labels = member_levels$member_short)
    )
  
  if (exists("condition_palette_semantic", mode = "function")) {
    cols <- condition_palette_semantic(sort(unique(dat$Target_Condition)))
  } else {
    cols <- condition_palette(sort(unique(dat$Target_Condition)))
  }
  
  summ <- dat %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$member_taxon, .data$Target_Condition) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_baseline = stats::median(.data$baseline, na.rm = TRUE),
      iqr_baseline = stats::IQR(.data$baseline, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "community_member_baseline_by_condition_summary.csv")))
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Target_Condition, y = .data$baseline_log, fill = .data$Target_Condition)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.18, outlier.alpha = 0.2) +
    ggplot2::facet_grid(tool + spike_label ~ member_taxon, scales = "free_y") +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(
      x = "Condition",
      y = "log10(member baseline mass + 1e-6)",
      title = "Community-member baseline relative abundance by condition",
      caption = "Baseline abundance is measured in the matched original samples before spiking."
    ) +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )
  
  save_pub_plot(
    p,
    file.path(outdir, paste0(prefix, "community_member_baseline_by_condition_log.png")),
    width = max(12, 1.8 * dplyr::n_distinct(dat$member_taxon) + 3),
    height = max(6.8, 2.3 * dplyr::n_distinct(interaction(dat$tool, dat$spike_label, drop = TRUE)) + 1.5)
  )
  invisible(NULL)
}

plot_community_member_baseline_by_condition_and_study <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df) || !has_study_column(target_df) || !"Target_Condition" %in% names(target_df)) return(invisible(NULL))
  dat <- target_df %>%
    dplyr::filter(!is.na(.data$Target_Condition), nzchar(.data$Target_Condition), !is.na(.data$Study), nzchar(.data$Study), !is.na(.data$member_taxon), nzchar(.data$member_taxon), is.finite(.data$baseline)) %>%
    dplyr::mutate(member_short = short_member_label(.data$member_taxon), baseline_log = log10(.data$baseline + 1e-6), cond_study = paste(.data$Target_Condition, .data$Study, sep = "||"))
  if (!nrow(dat)) return(invisible(NULL))
  member_levels <- dat %>%
    dplyr::group_by(.data$member_taxon) %>%
    dplyr::summarise(member_short = dplyr::first(.data$member_short), baseline_median = stats::median(.data$baseline, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$baseline_median), .data$member_taxon)
  dat <- dat %>% dplyr::mutate(member_taxon = factor(.data$member_taxon, levels = member_levels$member_taxon, labels = member_levels$member_short))
  pal <- condition_study_palette(dat$Target_Condition, dat$Study)
  lbls <- stats::setNames(gsub("||", " • ", names(pal), fixed = TRUE), names(pal))
  dodge <- ggplot2::position_dodge(width = 0.82)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Target_Condition, y = .data$baseline_log, fill = .data$cond_study, group = .data$cond_study)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.85, position = dodge, linewidth = 0.25) +
    ggplot2::geom_boxplot(width = 0.16, outlier.alpha = 0.15, position = dodge) +
    ggplot2::facet_grid(tool + spike_label ~ member_taxon, scales = "free_y") +
    ggplot2::scale_fill_manual(values = pal, breaks = names(pal), labels = lbls, drop = FALSE) +
    ggplot2::labs(x = "Condition", y = "log10(member baseline mass + 1e-6)", fill = "Condition • Study", title = "Community-member baseline abundance by condition and study") +
    theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pub_plot(p, file.path(outdir, paste0(prefix, "community_member_baseline_by_condition_and_study_log.png")), width = max(12, 1.8 * dplyr::n_distinct(dat$member_taxon) + 3), height = max(7.2, 2.3 * dplyr::n_distinct(interaction(dat$tool, dat$spike_label, drop = TRUE)) + 1.5))
  invisible(NULL)
}

generate_community_member_plots <- function(trace_df, target_df, outdir, prefix = "") {
  ensure_dir(outdir)
  if (nrow(target_df) > 0) {
    plot_community_member_baseline_by_condition(target_df, outdir, prefix = prefix)
    plot_community_member_detection_curves(target_df, outdir, prefix = prefix)
    plot_community_member_recovery_curves(target_df, outdir, prefix = prefix)
    plot_community_member_condition_detection(target_df, outdir, prefix = prefix)
    plot_community_member_baseline_by_condition_and_study(target_df, outdir, prefix = prefix)
  }
  if (nrow(trace_df) > 0) {
    plot_community_member_error_boxplots(trace_df, outdir, prefix = prefix)
  }
  invisible(NULL)
}

generate_all_spike_plots <- function(combined_df, baseline_df, trace_df, target_df, outdir, prefix = "") {
  ensure_dir(outdir)
  if (nrow(combined_df) > 0) plot_summary_curves(combined_df, outdir, prefix = prefix)
  if (nrow(target_df) > 0) {
    plot_lod_curve(target_df, outdir, prefix = prefix)
    plot_member_recovery_heatmap(target_df, outdir, prefix = prefix)
  }
  if (nrow(trace_df) > 0) {
    plot_target_member_error_boxplots(trace_df, outdir, prefix = prefix)
    plot_target_member_outlier_partitions(trace_df, outdir, prefix = prefix)
    plot_error_variance_scatter(trace_df, outdir, prefix = prefix)
    plot_target_dilution_checks(trace_df, outdir, prefix = prefix)
    plot_condition_relerr_violin(trace_df, outdir, prefix = prefix)
    plot_condition_error_bins(trace_df, outdir, prefix = prefix)
    plot_condition_relerr_violin_by_study(trace_df, outdir, prefix = prefix)
  }
  if (nrow(target_df) > 0) {
    plot_condition_lod_passrate(target_df, outdir, prefix = prefix)
    plot_condition_lod_passrate_by_study(target_df, outdir, prefix = prefix)
  }
  if (nrow(baseline_df) > 0) {
    plot_baseline_target_mass_by_condition(baseline_df, outdir, prefix = prefix)
    plot_baseline_target_mass_by_condition_and_study(baseline_df, outdir, prefix = prefix)
  }
  if (nrow(combined_df) > 0) {
    plot_bias_variance_summary(combined_df, outdir, prefix = prefix)
    plot_target_recovery_vs_baseline(combined_df, baseline_df, outdir, prefix = prefix)
    plot_summary_curves_by_condition_and_study(combined_df, outdir, prefix = prefix)
    plot_bias_variance_summary_by_condition_and_study(combined_df, outdir, prefix = prefix)
    plot_target_recovery_vs_baseline_by_condition_and_study(combined_df, baseline_df, outdir, prefix = prefix)
  }
  invisible(NULL)
}
