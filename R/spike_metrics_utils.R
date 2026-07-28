get_design_targets_for_sample <- function(design_sample, table_cols) {
  design_sample %>%
    dplyr::mutate(matched_taxon = purrr::map_chr(.data$member_taxon, match_taxon_to_column, table_cols = table_cols))
}

compute_spike_metrics_for_manifest_row <- function(man_row, design, meta = NULL) {
  spiked <- read_abundance_table(man_row$spiked_table[[1]])
  original <- read_abundance_table(man_row$original_table[[1]])

  design_sub <- design %>%
    dplyr::filter(
      .data$spike_label == man_row$spike_label[[1]],
      abs(.data$spike_fraction_total - man_row$spike_fraction[[1]]) < 1e-12
    )

  design_sub <- design_sub %>% dplyr::filter(.data$sample_id %in% spiked$sample_id)
  if (nrow(design_sub) == 0) {
    warning(sprintf("No design rows matched manifest entry: %s | %s | %s", man_row$spike_label[[1]], man_row$tool[[1]], man_row$spike_fraction[[1]]))
    return(NULL)
  }

  baseline_rows <- list()
  trace_rows <- list()
  summary_rows <- list()
  target_rows <- list()
  idx_b <- idx_t <- idx_s <- idx_tm <- 1L

  if (!is.null(meta)) {
    meta_sub <- meta %>% dplyr::filter(.data$sample_id %in% spiked$sample_id)
  } else {
    meta_sub <- design_sub %>%
      dplyr::distinct(.data$sample_id, .data$base_id, .data$Target_Condition, .data$original_id, .data$spike_mode, .data$spike_label, .data$spike_fraction_total)
  }

  for (sid in unique(design_sub$sample_id)) {
    d_s <- design_sub %>% dplyr::filter(.data$sample_id == !!sid)
    m_s <- meta_sub %>% dplyr::filter(.data$sample_id == !!sid)
    if (nrow(m_s) == 0) next
    base_id <- m_s$base_id[[1]]
    frac <- d_s$spike_fraction_total[[1]]

    sp_row <- spiked %>% dplyr::filter(.data$sample_id == !!sid)
    if (nrow(sp_row) == 0) next
    ob_row <- original %>% dplyr::filter(.data$sample_id == !!base_id)
    if (nrow(ob_row) == 0) {
      warning(sprintf("Base sample %s not found in original table %s", base_id, man_row$original_table[[1]]))
      next
    }

    sp_vec <- stats::setNames(as.numeric(sp_row[1, setdiff(names(sp_row), "sample_id"), drop = TRUE]), setdiff(names(sp_row), "sample_id"))
    ob_vec <- stats::setNames(as.numeric(ob_row[1, setdiff(names(ob_row), "sample_id"), drop = TRUE]), setdiff(names(ob_row), "sample_id"))

    universe <- union(names(ob_vec), names(sp_vec))
    target_map <- get_design_targets_for_sample(d_s, universe)
    unmatched_targets <- target_map %>% dplyr::filter(is.na(.data$matched_taxon))
    if (nrow(unmatched_targets) > 0) {
      universe <- union(universe, unmatched_targets$member_taxon)
      target_map <- get_design_targets_for_sample(d_s, universe)
    }

    universe <- unique(universe)
    base_vals <- stats::setNames(rep(0, length(universe)), universe)
    obs_vals <- stats::setNames(rep(0, length(universe)), universe)
    base_vals[names(ob_vec)] <- ob_vec
    obs_vals[names(sp_vec)] <- sp_vec
    exp_vals <- base_vals * (1 - frac)

    for (i in seq_len(nrow(target_map))) {
      tx <- if (!is.na(target_map$matched_taxon[[i]]) && nzchar(target_map$matched_taxon[[i]])) target_map$matched_taxon[[i]] else target_map$member_taxon[[i]]
      exp_vals[[tx]] <- exp_vals[[tx]] + target_map$member_fraction_expected[[i]]
    }

    target_taxa <- ifelse(is.na(target_map$matched_taxon) | !nzchar(target_map$matched_taxon), target_map$member_taxon, target_map$matched_taxon)
    target_taxa <- unique(target_taxa[!is.na(target_taxa)])
    target_lookup <- tibble::tibble(
      taxon = target_taxa,
      member_taxon = target_map$member_taxon[match(target_taxa, ifelse(is.na(target_map$matched_taxon) | !nzchar(target_map$matched_taxon), target_map$member_taxon, target_map$matched_taxon))]
    )

    trace_df <- tibble::tibble(
      sample_id = sid,
      base_id = base_id,
      spike_mode = d_s$spike_mode[[1]],
      spike_label = d_s$spike_label[[1]],
      spike_fraction_total = frac,
      tool = man_row$tool[[1]],
      taxon = universe,
      target_flag = universe %in% target_taxa,
      observed = as.numeric(obs_vals[universe]),
      baseline = as.numeric(base_vals[universe]),
      expected = as.numeric(exp_vals[universe])
    ) %>%
      dplyr::mutate(
        error = .data$observed - .data$expected,
        abs_error = abs(.data$error),
        relative_error = dplyr::if_else(.data$expected > 0, .data$error / .data$expected, NA_real_),
        detected = .data$observed > 0,
        observed_over_expected = dplyr::if_else(.data$expected > 0, .data$observed / .data$expected, NA_real_)
      ) %>%
      dplyr::left_join(target_lookup, by = "taxon")

    if (nrow(m_s) > 0) {
      extra_cols <- setdiff(names(m_s), c("sample_id", "base_id"))
      for (ec in extra_cols) trace_df[[ec]] <- m_s[[ec]][[1]]
    }

    trace_rows[[idx_t]] <- trace_df
    idx_t <- idx_t + 1L

    target_df <- trace_df %>% dplyr::filter(.data$target_flag)
    target_rows[[idx_tm]] <- target_df
    idx_tm <- idx_tm + 1L

    baseline_rows[[idx_b]] <- tibble::tibble(
      sample_id = sid,
      base_id = base_id,
      tool = man_row$tool[[1]],
      spike_label = d_s$spike_label[[1]],
      spike_mode = d_s$spike_mode[[1]],
      spike_fraction_total = frac,
      target_baseline_mass = sum(target_df$baseline, na.rm = TRUE),
      target_expected_mass = sum(target_df$expected, na.rm = TRUE),
      target_observed_mass = sum(target_df$observed, na.rm = TRUE)
    )
    if (nrow(m_s) > 0) {
      extra_cols <- setdiff(names(m_s), c("sample_id", "base_id"))
      for (ec in extra_cols) baseline_rows[[idx_b]][[ec]] <- m_s[[ec]][[1]]
    }
    idx_b <- idx_b + 1L

    off_df <- trace_df %>% dplyr::filter(!.data$target_flag)
    summary_df <- tibble::tibble(
      sample_id = sid,
      base_id = base_id,
      tool = man_row$tool[[1]],
      spike_label = d_s$spike_label[[1]],
      spike_mode = d_s$spike_mode[[1]],
      spike_fraction_total = frac,
      n_target_members = nrow(target_map),
      n_target_detected = sum(target_df$detected, na.rm = TRUE),
      target_detection_fraction = ifelse(nrow(target_map) > 0, sum(target_df$detected, na.rm = TRUE) / nrow(target_map), NA_real_),
      target_expected_mass = sum(target_df$expected, na.rm = TRUE),
      target_observed_mass = sum(target_df$observed, na.rm = TRUE),
      target_signed_error = sum(target_df$error, na.rm = TRUE),
      target_abs_error = sum(abs(target_df$error), na.rm = TRUE),
      target_recovery = if (sum(target_df$expected, na.rm = TRUE) > 0) sum(target_df$observed, na.rm = TRUE) / sum(target_df$expected, na.rm = TRUE) else NA_real_,
      target_mean_abs_rel_error = mean(abs(target_df$relative_error), na.rm = TRUE),
      off_target_abs_error = sum(abs(off_df$error), na.rm = TRUE),
      off_target_observed_mass = sum(off_df$observed, na.rm = TRUE),
      false_positive_mass = sum(off_df$observed[off_df$baseline <= 0], na.rm = TRUE),
      n_false_positive_taxa = sum(off_df$observed > 0 & off_df$baseline <= 0, na.rm = TRUE)
    )
    if (nrow(m_s) > 0) {
      extra_cols <- setdiff(names(m_s), c("sample_id", "base_id"))
      for (ec in extra_cols) summary_df[[ec]] <- m_s[[ec]][[1]]
    }
    summary_rows[[idx_s]] <- summary_df
    idx_s <- idx_s + 1L
  }

  list(
    baseline = dplyr::bind_rows(baseline_rows),
    trace = dplyr::bind_rows(trace_rows),
    combined = dplyr::bind_rows(summary_rows),
    target_member_errors = dplyr::bind_rows(target_rows)
  )
}

compute_lod_summary <- function(target_member_errors, detection_threshold = 0.95) {
  target_member_errors %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$taxon, .data$member_taxon, .data$spike_fraction_total) %>%
    dplyr::summarise(
      detection_rate = mean(.data$detected, na.rm = TRUE),
      mean_observed = mean(.data$observed, na.rm = TRUE),
      mean_expected = mean(.data$expected, na.rm = TRUE),
      mean_recovery = mean(.data$observed_over_expected, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$tool, .data$spike_label, .data$taxon, .data$member_taxon) %>%
    dplyr::summarise(
      lod_fraction = {
        good <- .data$spike_fraction_total[.data$detection_rate >= detection_threshold]
        if (length(good) == 0) NA_real_ else min(good)
      },
      .groups = "drop"
    )
}

summarise_species_errors <- function(trace_df) {
  if (!nrow(trace_df)) return(tibble::tibble())
  trace_df %>%
    dplyr::group_by(.data$tool, .data$spike_mode, .data$spike_label, .data$taxon, .data$member_taxon, .data$target_flag, .data$spike_fraction_total) %>%
    dplyr::summarise(
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      mean_abs_rel_error = mean(abs(.data$relative_error), na.rm = TRUE),
      sd_rel_error = stats::sd(.data$relative_error, na.rm = TRUE),
      median_observed = stats::median(.data$observed, na.rm = TRUE),
      median_expected = stats::median(.data$expected, na.rm = TRUE),
      .groups = "drop"
    )
}
