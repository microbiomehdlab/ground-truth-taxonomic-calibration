require_maaslin2 <- function() {
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("Package 'Maaslin2' is required for this script but is not installed.", call. = FALSE)
  }
}

make_table_matrix <- function(df) {
  mat <- df %>% as.data.frame()
  rownames(mat) <- mat$sample_id
  mat$sample_id <- NULL
  mat
}

compute_prevalence <- function(df) {
  taxa <- setdiff(names(df), "sample_id")
  tibble::tibble(
    taxon = taxa,
    prevalence = vapply(taxa, function(tx) mean(df[[tx]] > 0, na.rm = TRUE), numeric(1)),
    variance = vapply(taxa, function(tx) stats::var(df[[tx]], na.rm = TRUE), numeric(1))
  )
}

maaslin_alias_normalize_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[[:space:]]+", " ", x)
  x
}

maaslin_prepare_alias_table <- function(alias_table = NULL) {
  if (is.null(alias_table) || nrow(alias_table) == 0) {
    return(tibble::tibble(
      canonical = character(),
      alias = character(),
      tool = character(),
      spike_label = character(),
      canonical_key = character(),
      alias_key = character()
    ))
  }
  
  out <- alias_table
  if (!"tool" %in% names(out)) out$tool <- NA_character_
  if (!"spike_label" %in% names(out)) out$spike_label <- NA_character_
  
  out <- out %>%
    dplyr::mutate(
      canonical = trimws(as.character(.data$canonical)),
      alias = trimws(as.character(.data$alias)),
      tool = dplyr::if_else(is.na(.data$tool) | .data$tool == "", NA_character_, tolower(trimws(as.character(.data$tool)))),
      spike_label = dplyr::if_else(is.na(.data$spike_label) | .data$spike_label == "", NA_character_, trimws(as.character(.data$spike_label)))
    )
  
  if (!"canonical_key" %in% names(out)) out$canonical_key <- maaslin_alias_normalize_key(out$canonical)
  if (!"alias_key" %in% names(out)) out$alias_key <- maaslin_alias_normalize_key(out$alias)
  
  self_maps <- out %>%
    dplyr::distinct(.data$canonical, .data$tool, .data$spike_label) %>%
    dplyr::transmute(
      canonical = .data$canonical,
      alias = .data$canonical,
      tool = .data$tool,
      spike_label = .data$spike_label,
      canonical_key = maaslin_alias_normalize_key(.data$canonical),
      alias_key = maaslin_alias_normalize_key(.data$canonical)
    )
  
  dplyr::bind_rows(out, self_maps) %>% dplyr::distinct()
}

maaslin_build_alias_lookup <- function(alias_table) {
  alias_table <- maaslin_prepare_alias_table(alias_table)
  if (nrow(alias_table) == 0) return(setNames(character(), character()))
  
  lut_df <- alias_table %>%
    dplyr::mutate(
      tool_key = dplyr::if_else(is.na(.data$tool), "*", .data$tool),
      spike_key = dplyr::if_else(is.na(.data$spike_label), "*", .data$spike_label),
      specificity =
        dplyr::if_else(.data$tool_key != "*", 1L, 0L) +
        dplyr::if_else(.data$spike_key != "*", 1L, 0L),
      alias_pref = dplyr::if_else(.data$alias_key != .data$canonical_key, 1L, 0L),
      lookup_key = paste(.data$tool_key, .data$spike_key, .data$alias_key, sep = "||")
    ) %>%
    dplyr::arrange(dplyr::desc(.data$specificity), dplyr::desc(.data$alias_pref)) %>%
    dplyr::distinct(.data$lookup_key, .keep_all = TRUE)
  
  stats::setNames(lut_df$canonical, lut_df$lookup_key)
}

maaslin_canonicalize_taxa_vector <- function(x, alias_lookup, tool = NULL, spike_label = NULL) {
  if (is.null(x)) return(x)
  
  x_chr <- as.character(x)
  n <- length(x_chr)
  if (n == 0L || length(alias_lookup) == 0L) return(x_chr)
  
  tool_vec <- if (is.null(tool)) rep(NA_character_, n) else as.character(tool)
  spike_vec <- if (is.null(spike_label)) rep(NA_character_, n) else as.character(spike_label)
  
  if (length(tool_vec) == 1L && n > 1L) tool_vec <- rep(tool_vec, n)
  if (length(spike_vec) == 1L && n > 1L) spike_vec <- rep(spike_vec, n)
  
  if (length(tool_vec) != n) stop("tool must have length 1 or same length as x", call. = FALSE)
  if (length(spike_vec) != n) stop("spike_label must have length 1 or same length as x", call. = FALSE)
  
  missing_x <- is.na(x_chr) | !nzchar(trimws(x_chr))
  out <- x_chr
  
  x_key <- maaslin_alias_normalize_key(ifelse(missing_x, "", x_chr))
  tool_key <- ifelse(is.na(tool_vec) | !nzchar(trimws(tool_vec)), "*", tolower(trimws(tool_vec)))
  spike_key <- ifelse(is.na(spike_vec) | !nzchar(trimws(spike_vec)), "*", trimws(spike_vec))
  
  k1 <- paste(tool_key, spike_key, x_key, sep = "||")
  k2 <- paste(tool_key, "*", x_key, sep = "||")
  k3 <- paste("*", spike_key, x_key, sep = "||")
  k4 <- paste("*", "*", x_key, sep = "||")
  
  hit1 <- unname(alias_lookup[k1])
  hit2 <- unname(alias_lookup[k2])
  hit3 <- unname(alias_lookup[k3])
  hit4 <- unname(alias_lookup[k4])
  
  replace_idx <- !missing_x & !is.na(hit1)
  out[replace_idx] <- hit1[replace_idx]
  
  replace_idx <- !missing_x & is.na(hit1) & !is.na(hit2)
  out[replace_idx] <- hit2[replace_idx]
  
  replace_idx <- !missing_x & is.na(hit1) & is.na(hit2) & !is.na(hit3)
  out[replace_idx] <- hit3[replace_idx]
  
  replace_idx <- !missing_x & is.na(hit1) & is.na(hit2) & is.na(hit3) & !is.na(hit4)
  out[replace_idx] <- hit4[replace_idx]
  
  out
}

build_errorvar_filter <- function(metrics_dir, threshold = 0.01, alias_table = NULL) {
  err_path <- file.path(metrics_dir, "species_errors_with_condition.csv")
  stop_if_missing(err_path)
  err <- safe_read_csv(err_path)
  
  req <- c("tool", "taxon")
  miss <- setdiff(req, names(err))
  if (length(miss) > 0) {
    stop(
      sprintf(
        "species_errors_with_condition.csv is missing required columns: %s",
        paste(miss, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if ("mean_abs_error" %in% names(err)) {
    out <- err %>%
      dplyr::select(.data$tool, .data$taxon, .data$mean_abs_error) %>%
      dplyr::distinct()
  } else if ("abs_error" %in% names(err)) {
    out <- err %>%
      dplyr::group_by(.data$tool, .data$taxon) %>%
      dplyr::summarise(
        mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
        .groups = "drop"
      )
  } else if ("error" %in% names(err)) {
    out <- err %>%
      dplyr::group_by(.data$tool, .data$taxon) %>%
      dplyr::summarise(
        mean_abs_error = mean(abs(.data$error), na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    stop(
      "species_errors_with_condition.csv must contain one of: mean_abs_error, abs_error, or error",
      call. = FALSE
    )
  }
  
  alias_lookup <- maaslin_build_alias_lookup(alias_table)
  if (length(alias_lookup) > 0L) {
    out <- out %>%
      dplyr::mutate(
        tool = as.character(.data$tool),
        taxon = maaslin_canonicalize_taxa_vector(.data$taxon, alias_lookup, tool = .data$tool, spike_label = NULL)
      )
  }
  
  out %>%
    dplyr::filter(is.finite(.data$mean_abs_error)) %>%
    dplyr::filter(.data$mean_abs_error > threshold) %>%
    dplyr::distinct(.data$tool, .data$taxon, .keep_all = TRUE)
}

filter_abundance_table <- function(df, mode = "original", tool = NULL, metrics_dir = NULL,
                                   min_prevalence = 0.05, error_threshold = 0.01,
                                   alias_table = NULL, spike_label = NULL) {
  taxa <- setdiff(names(df), "sample_id")
  keep <- taxa
  audit <- tibble::tibble(taxon = taxa, keep = TRUE, reason = "kept")
  
  if (mode %in% c("std", "errorvar")) {
    prev <- compute_prevalence(df)
    drop_std <- prev %>%
      dplyr::filter(.data$prevalence < min_prevalence | is.na(.data$variance) | .data$variance <= 0) %>%
      dplyr::pull(.data$taxon)
    if (length(drop_std) > 0) {
      audit$keep[audit$taxon %in% drop_std] <- FALSE
      audit$reason[audit$taxon %in% drop_std] <- "std_filter"
    }
  }
  
  if (mode == "errorvar") {
    if (is.null(metrics_dir)) stop("metrics_dir is required for errorvar mode", call. = FALSE)
    
    ev <- build_errorvar_filter(metrics_dir, threshold = error_threshold, alias_table = alias_table)
    if (!is.null(tool)) ev <- ev %>% dplyr::filter(.data$tool == !!tool)
    
    alias_lookup <- maaslin_build_alias_lookup(alias_table)
    taxa_tbl <- tibble::tibble(
      taxon = taxa,
      taxon_canonical = if (length(alias_lookup) > 0L) {
        maaslin_canonicalize_taxa_vector(taxa, alias_lookup, tool = tool, spike_label = spike_label)
      } else {
        as.character(taxa)
      }
    )
    
    drop_canonical <- unique(as.character(ev$taxon))
    drop_ev <- taxa_tbl %>%
      dplyr::filter(.data$taxon_canonical %in% drop_canonical) %>%
      dplyr::pull(.data$taxon)
    
    if (length(drop_ev) > 0) {
      audit$keep[audit$taxon %in% drop_ev] <- FALSE
      audit$reason[audit$taxon %in% drop_ev] <- "errorvar_filter"
    }
  }
  
  keep <- audit %>% dplyr::filter(.data$keep) %>% dplyr::pull(.data$taxon)
  if (length(keep) == 0) stop("All taxa were filtered out", call. = FALSE)
  
  list(
    table = df %>% dplyr::select("sample_id", dplyr::all_of(keep)),
    audit = audit
  )
}

run_one_maaslin <- function(abundance_df, metadata_df, fixed_effect = "Target_Condition",
                            reference_level = "Control", outdir) {
  require_maaslin2()
  ensure_dir(outdir)
  input_data <- make_table_matrix(abundance_df)
  input_meta <- metadata_df %>% as.data.frame()
  
  if (!(fixed_effect %in% colnames(input_meta))) {
    stop(sprintf("Metadata is missing fixed effect column: %s", fixed_effect), call. = FALSE)
  }
  
  vals <- unique(as.character(input_meta[[fixed_effect]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) > 2) {
    if (!(reference_level %in% vals)) {
      stop(sprintf(
        "Reference level '%s' was not found in %s. Available levels: %s",
        reference_level, fixed_effect, paste(sort(vals), collapse = ", ")
      ), call. = FALSE)
    }
    input_meta[[fixed_effect]] <- stats::relevel(factor(input_meta[[fixed_effect]]), ref = reference_level)
    reference_arg <- paste0(fixed_effect, ",", reference_level)
  } else {
    input_meta[[fixed_effect]] <- factor(input_meta[[fixed_effect]])
    reference_arg <- NULL
  }
  
  rownames(input_meta) <- input_meta$sample_id
  input_meta$sample_id <- NULL
  
  maaslin_fun <- getExportedValue("Maaslin2", "Maaslin2")
  args <- list(
    input_data = input_data,
    input_metadata = input_meta,
    output = outdir,
    fixed_effects = fixed_effect,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  if (!is.null(reference_arg)) args$reference <- reference_arg
  
  do.call(maaslin_fun, args)
}

read_maaslin_results <- function(outdir) {
  fp <- file.path(outdir, "all_results.tsv")
  if (!file.exists(fp)) return(tibble::tibble())
  safe_read_tsv(fp)
}
