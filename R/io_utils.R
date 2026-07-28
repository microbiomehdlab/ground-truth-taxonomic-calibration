.abundance_scale_log_env <- new.env(parent = emptyenv())

detect_abundance_scale <- function(dat, max_fraction_value = 1.0, percent_sum_threshold = 1.5) {
  taxa <- setdiff(names(dat), "sample_id")
  if (length(taxa) == 0) {
    return(list(
      inferred_scale = "fraction",
      converted_to_fraction = FALSE,
      n_samples = nrow(dat),
      n_taxa = 0L,
      max_value = NA_real_,
      median_row_sum = NA_real_,
      p95_row_sum = NA_real_,
      reason = "no_taxa_columns"
    ))
  }
  
  mat <- as.matrix(dat[, taxa, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- NA_real_
  
  finite_vals <- as.numeric(mat)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  row_sums <- rowSums(mat, na.rm = TRUE)
  row_sums <- row_sums[is.finite(row_sums)]
  
  max_value <- if (length(finite_vals)) max(finite_vals) else NA_real_
  median_row_sum <- if (length(row_sums)) stats::median(row_sums) else NA_real_
  p95_row_sum <- if (length(row_sums)) as.numeric(stats::quantile(row_sums, probs = 0.95, na.rm = TRUE, names = FALSE)) else NA_real_
  
  inferred_scale <- dplyr::case_when(
    is.finite(max_value) && max_value > max_fraction_value + 1e-8 ~ "percent",
    is.finite(median_row_sum) && median_row_sum > percent_sum_threshold ~ "percent",
    TRUE ~ "fraction"
  )
  
  reason <- dplyr::case_when(
    identical(inferred_scale, "percent") && is.finite(max_value) && max_value > max_fraction_value + 1e-8 ~ "max_value_gt_1",
    identical(inferred_scale, "percent") && is.finite(median_row_sum) && median_row_sum > percent_sum_threshold ~ "median_row_sum_gt_1.5",
    TRUE ~ "already_fraction_like"
  )
  
  list(
    inferred_scale = inferred_scale,
    converted_to_fraction = identical(inferred_scale, "percent"),
    n_samples = nrow(dat),
    n_taxa = length(taxa),
    max_value = max_value,
    median_row_sum = median_row_sum,
    p95_row_sum = p95_row_sum,
    reason = reason
  )
}

maybe_log_abundance_scale <- function(path, info) {
  key <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!exists(key, envir = .abundance_scale_log_env, inherits = FALSE)) {
    msg <- sprintf(
      "[INFO] %s | inferred_scale=%s | converted_to_fraction=%s | max_value=%.6f | median_row_sum=%.6f | reason=%s",
      basename(path),
      info$inferred_scale,
      ifelse(isTRUE(info$converted_to_fraction), "yes", "no"),
      info$max_value %||% NA_real_,
      info$median_row_sum %||% NA_real_,
      info$reason %||% "unknown"
    )
    message(msg)
    assign(key, TRUE, envir = .abundance_scale_log_env)
  }
  invisible(NULL)
}

audit_abundance_table <- function(path, abundance_scale = c("auto", "fraction", "percent")) {
  abundance_scale <- match.arg(abundance_scale)
  stop_if_missing(path, "abundance table")
  ext <- tolower(tools::file_ext(path))
  dat <- if (ext %in% c("tsv", "txt")) safe_read_tsv(path) else safe_read_csv(path)
  if (ncol(dat) < 2) stop(sprintf("Table has <2 columns: %s", path), call. = FALSE)
  
  sample_col <- first_existing(names(dat), c("sample", "sample_id", "Sample", "SampleID", "ID"))
  if (is.na(sample_col)) sample_col <- names(dat)[[1]]
  
  dat <- dat %>% dplyr::rename(sample_id = !!sample_col)
  dat$sample_id <- as.character(dat$sample_id)
  
  for (nm in setdiff(names(dat), "sample_id")) {
    dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
    dat[[nm]][is.na(dat[[nm]])] <- 0
  }
  
  info <- detect_abundance_scale(dat)
  if (abundance_scale != "auto") {
    info$inferred_scale <- abundance_scale
    info$converted_to_fraction <- identical(abundance_scale, "percent")
    info$reason <- paste0("user_override_", abundance_scale)
  }
  
  tibble::tibble(
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    file = basename(path),
    inferred_scale = info$inferred_scale,
    converted_to_fraction = isTRUE(info$converted_to_fraction),
    n_samples = info$n_samples,
    n_taxa = info$n_taxa,
    max_value = info$max_value,
    median_row_sum = info$median_row_sum,
    p95_row_sum = info$p95_row_sum,
    reason = info$reason
  )
}

read_abundance_table <- function(path, abundance_scale = c("auto", "fraction", "percent")) {
  abundance_scale <- match.arg(abundance_scale)
  stop_if_missing(path, "abundance table")
  ext <- tolower(tools::file_ext(path))
  dat <- if (ext %in% c("tsv", "txt")) safe_read_tsv(path) else safe_read_csv(path)
  if (ncol(dat) < 2) stop(sprintf("Table has <2 columns: %s", path), call. = FALSE)
  
  sample_col <- first_existing(names(dat), c("sample", "sample_id", "Sample", "SampleID", "ID"))
  if (is.na(sample_col)) sample_col <- names(dat)[[1]]
  
  dat <- dat %>% dplyr::rename(sample_id = !!sample_col)
  dat$sample_id <- as.character(dat$sample_id)
  
  for (nm in setdiff(names(dat), "sample_id")) {
    dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
    dat[[nm]][is.na(dat[[nm]])] <- 0
  }
  
  info <- detect_abundance_scale(dat)
  if (abundance_scale != "auto") {
    info$inferred_scale <- abundance_scale
    info$converted_to_fraction <- identical(abundance_scale, "percent")
    info$reason <- paste0("user_override_", abundance_scale)
  }
  
  if (isTRUE(info$converted_to_fraction)) {
    taxa <- setdiff(names(dat), "sample_id")
    dat[taxa] <- lapply(dat[taxa], function(x) x / 100)
  }
  
  attr(dat, "abundance_scale_input") <- info$inferred_scale
  attr(dat, "abundance_scale_converted") <- isTRUE(info$converted_to_fraction)
  attr(dat, "abundance_scale_reason") <- info$reason
  maybe_log_abundance_scale(path, info)
  dat
}

read_metadata_table <- function(path) {
  stop_if_missing(path, "metadata")
  ext <- tolower(tools::file_ext(path))
  dat <- if (ext %in% c("tsv", "txt")) safe_read_tsv(path) else safe_read_csv(path)
  
  sid <- first_existing(names(dat), c("sample_id", "sample", "Sample", "SampleID", "ID", "base_id"))
  if (is.na(sid)) stop("Metadata must contain sample_id/sample column", call. = FALSE)
  
  dat <- dat %>%
    dplyr::rename(sample_id = !!sid) %>%
    dplyr::mutate(sample_id = as.character(.data$sample_id))
  
  if (!"base_id" %in% names(dat)) dat$base_id <- dat$sample_id
  if (!"Target_Condition" %in% names(dat)) dat$Target_Condition <- NA_character_
  if (!"original_id" %in% names(dat)) dat$original_id <- NA_character_
  
  dat <- standardize_study_column(
    dat,
    sample_col = "sample_id",
    base_col = "base_id",
    fallback_cols = c("Study", "Dataset", "Cohort", "Project", "original_id")
  )
  
  dat
}

resolve_original_tables <- function(original_tables_arg) {
  tbls <- parse_kv_csv(original_tables_arg)
  if (length(tbls) == 0) return(tibble::tibble())
  tibble::tibble(tool = names(tbls), original_table = unname(tbls))
}

build_run_manifest_from_profile_root <- function(profile_root, original_tables_arg, design = NULL) {
  stop_if_missing(profile_root, "profile_root")
  originals <- resolve_original_tables(original_tables_arg)
  if (nrow(originals) == 0) {
    stop("--original_tables is required when using --profile_root", call. = FALSE)
  }
  
  labels <- list.dirs(profile_root, recursive = FALSE, full.names = FALSE)
  rows <- list()
  idx <- 1L
  
  mode_lookup <- NULL
  if (!is.null(design) && nrow(design) > 0) {
    mode_lookup <- design %>%
      dplyr::distinct(.data$spike_label, .data$spike_mode)
  }
  
  for (lab in labels) {
    lab_dir <- file.path(profile_root, lab)
    tools <- list.dirs(lab_dir, recursive = FALSE, full.names = FALSE)
    for (tool in tools) {
      tool_dir <- file.path(lab_dir, tool)
      files <- list.files(tool_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
      files <- files[grepl("spike_.*closed_world|spike_f|merged_[0-9]", basename(files), ignore.case = TRUE)]
      if (length(files) == 0) next
      
      original_table <- originals %>%
        dplyr::filter(.data$tool == !!tool) %>%
        dplyr::pull(.data$original_table)
      if (length(original_table) == 0) next
      original_table <- original_table[[1]]
      
      spike_mode <- NA_character_
      if (!is.null(mode_lookup)) {
        modes <- mode_lookup %>%
          dplyr::filter(.data$spike_label == !!lab) %>%
          dplyr::pull(.data$spike_mode)
        if (length(unique(stats::na.omit(modes))) == 1) {
          spike_mode <- unique(stats::na.omit(modes))[[1]]
        }
      }
      
      for (f in files) {
        frac <- extract_fraction_from_string(basename(f))
        rows[[idx]] <- tibble::tibble(
          spike_label = lab,
          spike_mode = spike_mode,
          tool = tool,
          spike_fraction = frac,
          spiked_table = normalizePath(f, winslash = "/", mustWork = FALSE),
          original_table = normalizePath(original_table, winslash = "/", mustWork = FALSE)
        )
        idx <- idx + 1L
      }
    }
  }
  
  dplyr::bind_rows(rows) %>%
    dplyr::arrange(.data$spike_mode, .data$spike_label, .data$tool, .data$spike_fraction)
}

resolve_run_manifest <- function(run_manifest = NULL, profile_root = NULL, original_tables_arg = NULL, design = NULL) {
  mode_lookup <- NULL
  if (!is.null(design) && nrow(design) > 0) {
    mode_lookup <- design %>%
      dplyr::distinct(.data$spike_label, .data$spike_mode)
  }
  
  if (!is.null(run_manifest) && nzchar(run_manifest)) {
    man <- if (tolower(tools::file_ext(run_manifest)) %in% c("tsv", "txt")) {
      safe_read_tsv(run_manifest)
    } else {
      safe_read_csv(run_manifest)
    }
    
    req <- c("spike_label", "tool", "spiked_table")
    miss <- setdiff(req, names(man))
    if (length(miss) > 0) {
      stop(sprintf("run_manifest missing columns: %s", paste(miss, collapse = ", ")), call. = FALSE)
    }
    
    if (!"spike_fraction" %in% names(man)) man$spike_fraction <- extract_fraction_from_string(man$spiked_table)
    if (!"original_table" %in% names(man)) {
      originals <- resolve_original_tables(original_tables_arg)
      man <- man %>% dplyr::left_join(originals, by = "tool")
    }
    if (!"spike_mode" %in% names(man)) man$spike_mode <- NA_character_
    
    if (!is.null(mode_lookup)) {
      man <- man %>%
        dplyr::select(-dplyr::any_of("spike_mode")) %>%
        dplyr::left_join(mode_lookup, by = "spike_label")
    }
    return(man)
  }
  
  if (!is.null(profile_root) && nzchar(profile_root)) {
    return(build_run_manifest_from_profile_root(profile_root, original_tables_arg, design = design))
  }
  
  stop("Provide either --run_manifest or --profile_root", call. = FALSE)
}
