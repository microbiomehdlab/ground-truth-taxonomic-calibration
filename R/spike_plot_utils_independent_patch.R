# Source this AFTER source("R/spike_plot_utils.R")

if (!exists("short_member_label", mode = "function")) {
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
}

get_independent_target_rows <- function(target_df) {
  if (!nrow(target_df)) return(target_df[0, , drop = FALSE])

  if ("spike_mode" %in% names(target_df)) {
    out <- target_df %>% dplyr::filter(.data$spike_mode == "independent")
    if (nrow(out)) return(out)
  }

  one_member_labels <- target_df %>%
    dplyr::filter(
      !is.na(.data$spike_label), nzchar(.data$spike_label),
      !is.na(.data$member_taxon), nzchar(.data$member_taxon)
    ) %>%
    dplyr::distinct(.data$spike_label, .data$member_taxon) %>%
    dplyr::count(.data$spike_label, name = "n_members") %>%
    dplyr::filter(.data$n_members == 1L) %>%
    dplyr::pull(.data$spike_label)

  target_df %>% dplyr::filter(.data$spike_label %in% one_member_labels)
}

plot_independent_member_error_boxplots <- function(target_df, outdir, prefix = "") {
  if (!nrow(target_df)) return(invisible(NULL))

  dat <- get_independent_target_rows(target_df) %>%
    dplyr::filter(
      is.finite(.data$relative_error),
      is.finite(.data$spike_fraction_total),
      .data$spike_fraction_total > 0,
      !is.na(.data$member_taxon), nzchar(.data$member_taxon)
    ) %>%
    dplyr::mutate(
      member_short = short_member_label(.data$member_taxon)
    ) %>%
    make_spike_factor()

  if (!nrow(dat)) return(invisible(NULL))

  member_map <- dat %>%
    dplyr::distinct(.data$member_taxon, .data$member_short) %>%
    dplyr::arrange(.data$member_short)

  dat <- dat %>%
    dplyr::mutate(
      member_short = factor(.data$member_short, levels = member_map$member_short)
    )

  member_cols <- condition_palette(as.character(levels(dat$member_short)))
  dodge <- ggplot2::position_dodge2(width = 0.86, preserve = "single", padding = 0.15)

  summ <- dat %>%
    dplyr::group_by(.data$tool, .data$member_taxon, .data$member_short, .data$spike_fraction_total) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_relative_error = stats::median(.data$relative_error, na.rm = TRUE),
      iqr_relative_error = stats::IQR(.data$relative_error, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(summ, file.path(outdir, paste0(prefix, "independent_member_error_summary.csv")))

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$spike_factor,
      y = .data$relative_error,
      fill = .data$member_short,
      group = interaction(.data$spike_factor, .data$member_short, drop = TRUE)
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#666666") +
    ggplot2::geom_boxplot(outlier.alpha = 0.18, position = dodge, width = 0.78) +
    ggplot2::facet_grid(tool ~ ., scales = "free_y") +
    ggplot2::scale_fill_manual(values = member_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Spike fraction",
      y = "Relative error",
      fill = "Spiked species",
      title = "Independent-spike relative error distributions",
      caption = "Relative error is computed for each independently spiked target species: (observed − expected) / expected."
    ) +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  save_pub_plot(
    p,
    file.path(outdir, paste0(prefix, "independent_member_error_boxplot.png")),
    width = 14,
    height = 7.2
  )

  dat_zoom <- dat %>% dplyr::filter(.data$relative_error >= -0.25, .data$relative_error <= 0.25)
  if (nrow(dat_zoom)) {
    pz <- ggplot2::ggplot(
      dat_zoom,
      ggplot2::aes(
        x = .data$spike_factor,
        y = .data$relative_error,
        fill = .data$member_short,
        group = interaction(.data$spike_factor, .data$member_short, drop = TRUE)
      )
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = "#666666") +
      ggplot2::geom_boxplot(outlier.alpha = 0.18, position = dodge, width = 0.78) +
      ggplot2::facet_grid(tool ~ ., scales = "free_y") +
      ggplot2::scale_fill_manual(values = member_cols, drop = FALSE) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.25, 0.25)) +
      ggplot2::labs(
        x = "Spike fraction",
        y = "Relative error",
        fill = "Spiked species",
        title = "Independent-spike relative error distributions (zoomed to ±25%)"
      ) +
      theme_pub() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    save_pub_plot(
      pz,
      file.path(outdir, paste0(prefix, "independent_member_error_boxplot_zoom.png")),
      width = 14,
      height = 7.2
    )
  }

  invisible(NULL)
}

generate_independent_member_plots <- function(target_df, outdir, prefix = "") {
  ensure_dir(outdir)
  if (nrow(target_df) > 0) {
    plot_independent_member_error_boxplots(target_df, outdir, prefix = prefix)
  }
  invisible(NULL)
}

if (exists("generate_all_spike_plots", mode = "function")) {
  .generate_all_spike_plots_original <- generate_all_spike_plots

  generate_all_spike_plots <- function(combined_df, baseline_df, trace_df, target_df, outdir, prefix = "") {
    .generate_all_spike_plots_original(
      combined_df = combined_df,
      baseline_df = baseline_df,
      trace_df = trace_df,
      target_df = target_df,
      outdir = outdir,
      prefix = prefix
    )
    if (nrow(target_df) > 0) {
      plot_independent_member_error_boxplots(target_df, outdir, prefix = prefix)
    }
    invisible(NULL)
  }
}
