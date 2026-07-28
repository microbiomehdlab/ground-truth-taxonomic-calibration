suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

stop_if_missing <- function(path, label = NULL) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("Missing %s: %s", label %||% "file", path %||% "<NULL>"), call. = FALSE)
  }
}

write_csv_safe <- function(x, path) {
  ensure_dir(dirname(path))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    readr::write_tsv(x, path)
  } else {
    readr::write_csv(x, path)
  }
  invisible(path)
}

safe_read_table_auto <- function(path) {
  stop_if_missing(path)
  first_line <- readLines(path, n = 1, warn = FALSE)
  ext <- tolower(tools::file_ext(path))

  guess_delim <- NULL
  if (length(first_line) > 0) {
    has_tab <- grepl("\t", first_line[[1]], fixed = FALSE)
    has_comma <- grepl(",", first_line[[1]], fixed = TRUE)
    if (has_tab && !has_comma) guess_delim <- "tab"
    if (has_comma && !has_tab) guess_delim <- "comma"
  }

  if (is.null(guess_delim)) {
    guess_delim <- if (ext %in% c("tsv", "txt")) "tab" else "comma"
  }

  primary <- if (guess_delim == "tab") {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  }

  if (ncol(primary) > 1) return(primary)

  alternate <- tryCatch(
    if (guess_delim == "tab") {
      readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    } else {
      readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    },
    error = function(e) primary
  )

  if (ncol(alternate) > ncol(primary)) return(alternate)
  primary
}

parse_kv_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character())
  items <- strsplit(x, ",", fixed = TRUE)[[1]]
  items <- trimws(items)
  items <- items[nzchar(items)]
  out <- character(length(items))
  nms <- character(length(items))
  for (i in seq_along(items)) {
    parts <- strsplit(items[[i]], "=", fixed = TRUE)[[1]]
    if (length(parts) < 2) stop(sprintf("Invalid key=value item: %s", items[[i]]), call. = FALSE)
    nms[[i]] <- trimws(parts[[1]])
    out[[i]] <- trimws(paste(parts[-1], collapse = "="))
  }
  names(out) <- nms
  out
}

parse_fraction_tag <- function(x) {
  ifelse(
    is.na(x),
    NA_real_,
    suppressWarnings(as.numeric(gsub("p", ".", sub("^f", "", x))))
  )
}

fraction_to_tag <- function(x) {
  s <- format(round(as.numeric(x), 6), scientific = FALSE, trim = TRUE)
  s <- sub("0+$", "", s)
  s <- sub("\\.$", "", s)
  paste0("f", gsub("\\.", "p", s))
}

extract_fraction_from_string <- function(x) {
  m <- stringr::str_match(x, "(f[0-9p]+)")
  parse_fraction_tag(m[, 2])
}

clean_taxon_text <- function(x) {
  x <- x %>%
    as.character() %>%
    stringr::str_replace_all("[|;]", " ") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace_all("__", " ") %>%
    stringr::str_replace_all("\\b[kpcofgst]__", " ") %>%
    stringr::str_replace_all("[^A-Za-z0-9 ]", " ") %>%
    stringr::str_squish() %>%
    tolower()
  x
}

normalize_taxon_full <- function(x) {
  clean_taxon_text(x)
}

normalize_taxon_core <- function(x) {
  x <- clean_taxon_text(x)
  toks <- strsplit(x, " +")
  vapply(toks, function(tok) {
    tok <- tok[nzchar(tok)]
    if (length(tok) == 0) return("")
    if (length(tok) == 1) return(tok[[1]])
    paste(tok[[1]], tok[[2]])
  }, character(1))
}

match_taxon_to_column <- function(target_taxon, table_cols) {
  if (length(table_cols) == 0) return(NA_character_)
  tfull <- normalize_taxon_full(target_taxon)
  tcore <- normalize_taxon_core(target_taxon)
  cfull <- normalize_taxon_full(table_cols)
  ccore <- normalize_taxon_core(table_cols)

  idx <- which(cfull == tfull)
  if (length(idx) == 1) return(table_cols[[idx]])
  idx <- which(ccore == tcore)
  if (length(idx) >= 1) return(table_cols[[idx[[1]]]])
  NA_character_
}

first_existing <- function(x, choices) {
  hit <- choices[choices %in% x]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_read_tsv <- function(path) {
  safe_read_table_auto(path)
}

safe_read_csv <- function(path) {
  safe_read_table_auto(path)
}

sanitize_slug <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", as.character(x))
}



infer_study_from_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  out <- rep(NA_character_, length(x))
  keep <- !is.na(x) & nzchar(x)
  if (!any(keep)) return(out)

  vals <- x[keep]
  pref <- stringr::str_extract(vals, "^[A-Za-z]+[0-9]+")
  pref[is.na(pref) | !nzchar(pref)] <- stringr::str_extract(vals[is.na(pref) | !nzchar(pref)], "^[A-Za-z]+")
  pref[is.na(pref) | !nzchar(pref)] <- vals[is.na(pref) | !nzchar(pref)]
  out[keep] <- pref
  out
}

standardize_study_column <- function(dat, sample_col = "sample_id", base_col = "base_id", fallback_cols = c("original_id")) {
  if (is.null(dat) || !nrow(dat)) {
    if (!"Study" %in% names(dat)) dat$Study <- character()
    return(dat)
  }

  candidate <- first_existing(names(dat), c("Study", "study", "Dataset", "dataset", "Cohort", "cohort", "Project", "project"))
  if (!is.na(candidate)) {
    dat$Study <- as.character(dat[[candidate]])
  } else {
    dat$Study <- NA_character_
  }

  need <- is.na(dat$Study) | !nzchar(trimws(dat$Study))
  if (any(need) && !is.null(base_col) && base_col %in% names(dat)) {
    dat$Study[need] <- infer_study_from_id(dat[[base_col]][need])
    need <- is.na(dat$Study) | !nzchar(trimws(dat$Study))
  }
  if (any(need) && !is.null(sample_col) && sample_col %in% names(dat)) {
    dat$Study[need] <- infer_study_from_id(dat[[sample_col]][need])
    need <- is.na(dat$Study) | !nzchar(trimws(dat$Study))
  }
  for (fc in fallback_cols) {
    if (!any(need)) break
    if (!is.null(fc) && fc %in% names(dat)) {
      dat$Study[need] <- infer_study_from_id(dat[[fc]][need])
      need <- is.na(dat$Study) | !nzchar(trimws(dat$Study))
    }
  }

  dat$Study <- as.character(dat$Study)
  dat
}
