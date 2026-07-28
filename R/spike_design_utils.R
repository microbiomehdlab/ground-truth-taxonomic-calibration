parse_spiked_sample_id <- function(sample_id) {
  m <- stringr::str_match(sample_id, "^(.*?)_([^_]+)_(f[0-9p]+)$")
  if (all(is.na(m))) {
    tibble::tibble(
      sample_id = sample_id,
      base_id = sample_id,
      spike_label = NA_character_,
      spike_fraction_total = NA_real_
    )
  } else {
    tibble::tibble(
      sample_id = sample_id,
      base_id = m[, 2],
      spike_label = m[, 3],
      spike_fraction_total = parse_fraction_tag(m[, 4])
    )
  }
}

read_manifest_generic <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  
  dat <- if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
    safe_read_tsv(path)
  } else {
    safe_read_csv(path)
  }
  
  sid <- first_existing(names(dat), c("sample_id", "sample", "Sample", "SampleID"))
  if (is.na(sid)) stop(sprintf("Manifest lacks sample_id column: %s", path), call. = FALSE)
  
  dat <- dat %>%
    dplyr::rename(sample_id = !!sid) %>%
    dplyr::mutate(sample_id = as.character(.data$sample_id))
  
  if (!"original_id" %in% names(dat)) dat$original_id <- NA_character_
  if (!"Target_Condition" %in% names(dat)) dat$Target_Condition <- NA_character_
  if (!"Study" %in% names(dat)) dat$Study <- NA_character_
  
  parsed <- parse_spiked_sample_id(dat$sample_id)
  dat <- dplyr::bind_cols(dat, parsed %>% dplyr::select(-"sample_id"))
  
  # Preserve an existing Study column; only standardize/fill if needed
  dat <- standardize_study_column(
    dat,
    sample_col = "sample_id",
    base_col = "base_id",
    fallback_cols = c("Study", "original_id")
  )
  
  if (!"Study" %in% names(dat)) dat$Study <- NA_character_
  dat$Study <- as.character(dat$Study)
  
  dat
}

read_spike_panel <- function(path) {
  pan <- if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
    safe_read_tsv(path)
  } else {
    safe_read_csv(path)
  }
  
  req <- c("label", "taxon_name")
  miss <- setdiff(req, names(pan))
  if (length(miss) > 0) {
    stop(sprintf("spike_panel missing columns: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }
  
  if (!"weight" %in% names(pan)) pan$weight <- 1
  
  pan %>%
    dplyr::mutate(
      label = as.character(.data$label),
      taxon_name = as.character(.data$taxon_name),
      weight = as.numeric(.data$weight %||% 1)
    )
}

read_community_memberships <- function(path, panel, community_labels = NULL) {
  if (!is.null(path) && nzchar(path)) {
    cm <- if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
      safe_read_tsv(path)
    } else {
      safe_read_csv(path)
    }
    
    req <- c("community_label", "label")
    miss <- setdiff(req, names(cm))
    if (length(miss) > 0) {
      stop(sprintf("community_memberships missing columns: %s", paste(miss, collapse = ", ")), call. = FALSE)
    }
    
    if (!"weight" %in% names(cm)) cm$weight <- 1
    
    cm <- cm %>%
      dplyr::mutate(
        community_label = as.character(.data$community_label),
        label = as.character(.data$label),
        weight = as.numeric(.data$weight)
      ) %>%
      dplyr::left_join(
        panel %>% dplyr::select("label", member_taxon = "taxon_name"),
        by = "label"
      )
    
    if (any(is.na(cm$member_taxon))) {
      bad <- cm %>%
        dplyr::filter(is.na(.data$member_taxon)) %>%
        dplyr::pull("label") %>%
        unique()
      stop(sprintf("Unknown labels in community_memberships: %s", paste(bad, collapse = ", ")), call. = FALSE)
    }
    
    return(cm)
  }
  
  if (is.null(community_labels) || length(community_labels) == 0) {
    return(tibble::tibble(
      community_label = character(),
      label = character(),
      weight = numeric(),
      member_taxon = character()
    ))
  }
  
  dplyr::bind_rows(lapply(community_labels, function(cl) {
    panel %>%
      dplyr::transmute(
        community_label = cl,
        label = .data$label,
        weight = .data$weight %||% 1,
        member_taxon = .data$taxon_name
      )
  }))
}

build_spike_design <- function(independent_manifest = NULL,
                               community_manifest = NULL,
                               spike_panel,
                               community_memberships = NULL) {
  panel <- read_spike_panel(spike_panel)
  ind <- read_manifest_generic(independent_manifest)
  com <- read_manifest_generic(community_manifest)
  
  community_labels <- character()
  if (!is.null(com) && nrow(com) > 0) {
    community_labels <- unique(stats::na.omit(com$spike_label))
  }
  
  cm <- read_community_memberships(community_memberships, panel, community_labels = community_labels)
  
  rows <- list()
  idx <- 1L
  
  if (!is.null(ind) && nrow(ind) > 0) {
    ind2 <- ind %>%
      dplyr::left_join(
        panel %>% dplyr::select("label", member_taxon = "taxon_name", member_weight = "weight"),
        by = c("spike_label" = "label")
      ) %>%
      dplyr::mutate(
        spike_mode = "independent",
        member_label = .data$spike_label,
        member_fraction_expected = .data$spike_fraction_total
      )
    
    rows[[idx]] <- ind2 %>%
      dplyr::transmute(
        sample_id = .data$sample_id,
        base_id = .data$base_id,
        spike_mode = .data$spike_mode,
        spike_label = .data$spike_label,
        spike_fraction_total = .data$spike_fraction_total,
        member_label = .data$member_label,
        member_taxon = .data$member_taxon,
        member_weight = .data$member_weight,
        member_fraction_expected = .data$member_fraction_expected,
        original_id = .data$original_id,
        Target_Condition = .data$Target_Condition,
        Study = .data$Study
      )
    
    idx <- idx + 1L
  }
  
  if (!is.null(com) && nrow(com) > 0) {
    if (nrow(cm) == 0) {
      stop("Community manifest provided but no community memberships available", call. = FALSE)
    }
    
    com2 <- com %>%
      dplyr::left_join(
        cm,
        by = c("spike_label" = "community_label"),
        relationship = "many-to-many"
      ) %>%
      dplyr::group_by(.data$sample_id) %>%
      dplyr::mutate(weight_norm = .data$weight / sum(.data$weight, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        spike_mode = "community",
        member_fraction_expected = .data$spike_fraction_total * .data$weight_norm,
        member_weight = .data$weight
      )
    
    rows[[idx]] <- com2 %>%
      dplyr::transmute(
        sample_id = .data$sample_id,
        base_id = .data$base_id,
        spike_mode = .data$spike_mode,
        spike_label = .data$spike_label,
        spike_fraction_total = .data$spike_fraction_total,
        member_label = .data$label,
        member_taxon = .data$member_taxon,
        member_weight = .data$member_weight,
        member_fraction_expected = .data$member_fraction_expected,
        original_id = .data$original_id,
        Target_Condition = .data$Target_Condition,
        Study = .data$Study
      )
  }
  
  design <- dplyr::bind_rows(rows) %>%
    dplyr::arrange(.data$sample_id, .data$member_taxon)
  
  meta <- design %>%
    dplyr::distinct(
      .data$sample_id,
      .data$base_id,
      .data$original_id,
      .data$Target_Condition,
      .data$Study,
      .data$spike_mode,
      .data$spike_label,
      .data$spike_fraction_total
    ) %>%
    dplyr::arrange(.data$sample_id)
  
  list(
    design = design,
    meta = meta,
    community_memberships = cm
  )
}