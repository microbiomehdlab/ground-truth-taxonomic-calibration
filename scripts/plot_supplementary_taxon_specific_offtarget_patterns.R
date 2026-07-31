#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
})

option_list <- list(
  make_option("--significant-file", type = "character",
              default = "RUNS_publication/maaslin_spike/maaslin_significant_features_ALLFILTERS.csv",
              help = "MaAsLin significant-feature summary [default %default]"),
  make_option("--artifact-overlap", type = "character",
              default = "RUNS/revised_manuscript_figures/source_panels/Fig5_community_offtarget_artifacts/abundance_artifact_DA_overlap_all_requested_fractions.csv",
              help = "Community non-target abundance-artifact/DA-overlap table [default %default]"),
  make_option("--run-manifest", type = "character",
              default = "RUNS_publication/maaslin_spike/run_design_manifest_ALLFILTERS.csv",
              help = "Optional run manifest; used to improve denominator counts [default %default]"),
  make_option("--alias-file", type = "character",
              default = "RUNS_publication/maaslin_spike/taxon_aliases.resolved.csv",
              help = "Optional target alias file [default %default]"),
  make_option("--outdir", type = "character",
              default = "RUNS_publication/manuscript_figures/supplementary_offtarget_patterns",
              help = "Output directory [default %default]"),
  make_option("--q-threshold", type = "double", default = 0.10,
              help = "Maximum q value for an enriched off-target call [default %default]"),
  make_option("--filter-mode", type = "character", default = "original",
              help = "Filter mode to retain from the original analysis [default %default]"),
  make_option("--community-label", type = "character", default = "CRCpanel",
              help = "Community spike label [default %default]"),
  make_option("--community-size", type = "integer", default = 10,
              help = "Number of taxa in the community spike [default %default]"),
  make_option("--top-features", type = "integer", default = 20,
              help = "Maximum recurrent off-target features per profiler in the heatmap [default %default]"),
  make_option("--min-set-recurrence", type = "double", default = 0.20,
              help = "Minimum within-set recurrence used for UpSet membership [default %default]"),
  make_option("--top-intersections", type = "integer", default = 18,
              help = "Maximum UpSet intersections shown per profiler [default %default]"),
  make_option("--focus-targets", type = "character", default = "Bfrag,Fnuc,Pmic,Pint",
              help = "Comma-separated targets for focused bubble heatmaps [default %default]"),
  make_option("--width", type = "double", default = 14.0,
              help = "Composite figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 16.0,
              help = "Composite figure height in inches [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}
stop_if_missing(opt$`significant-file`, "Significant-feature file")
stop_if_missing(opt$`artifact-overlap`, "Artifact-overlap file")

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

pick_col <- function(df, candidates, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) return(hit[[1]])
  if (required) stop("Could not identify ", label %||% paste(candidates, collapse = "/"),
                     ". Available columns: ", paste(names(df), collapse = ", "), call. = FALSE)
  NULL
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  y <- str_to_lower(trimws(as.character(x)))
  replace_na(y %in% c("true", "t", "1", "yes", "y"), FALSE)
}

norm_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\|", ";")
  x <- str_replace(x, ".*;s__", "")
  x <- str_replace(x, ".*s__", "")
  x <- str_replace_all(x, "[_]+", " ")
  x <- str_replace_all(x, "[;|]", " ")
  x <- str_squish(x)
  x
}

norm_key <- function(x) {
  x %>% norm_text() %>% str_to_lower() %>% str_replace_all("[^a-z0-9]+", " ") %>% str_squish()
}

display_feature <- function(x, width = 40) {
  z <- norm_text(x)
  # Keep database identifiers recognizable but make taxonomic labels readable.
  z <- if_else(
    str_detect(str_to_lower(z), "^(mgyg|gca|gcf)[0-9]"),
    str_to_upper(z),
    str_to_sentence(z)
  )
  str_trunc(z, width)
}

extract_genus <- function(x) {
  z <- norm_text(x)
  z <- str_replace(z, "^(uncultured|unclassified|unknown|candidate)\\s+", "")
  g <- word(z, 1)
  g[is.na(g) | !str_detect(g, "^[A-Za-z][A-Za-z0-9.-]+$")] <- "Unresolved"
  str_to_sentence(g)
}

pretty_tool <- function(x) {
  recode(as.character(x),
         kraken2_bracken = "Kraken2 + Bracken",
         kraken2 = "Kraken2 + Bracken",
         bracken = "Kraken2 + Bracken",
         metaphlan4 = "MetaPhlAn 4",
         metaphlan = "MetaPhlAn 4",
         .default = as.character(x))
}

fmt_fraction <- function(x) {
  x <- as.numeric(x)
  out <- ifelse(x < 0.0001,
                sprintf("%.4f%%", 100 * x),
                ifelse(x < 0.001, sprintf("%.3f%%", 100 * x), sprintf("%.2f%%", 100 * x)))
  sub("0+$", "", sub("\\.$", "", out))
}

target_map <- tibble::tribble(
  ~spike_label, ~target_taxon, ~target_genus,
  "Bfrag", "Bacteroides fragilis", "Bacteroides",
  "Csym",  "Clostridium symbiosum", "Clostridium",
  "Dpne",  "Dialister pneumosintes", "Dialister",
  "Fnuc",  "Fusobacterium nucleatum", "Fusobacterium",
  "Hhat",  "Hungatella hathewayi", "Hungatella",
  "Pmic",  "Parvimonas micra", "Parvimonas",
  "Pana",  "Peptostreptococcus anaerobius", "Peptostreptococcus",
  "Psto",  "Peptostreptococcus stomatis", "Peptostreptococcus",
  "Porp",  "Porphyromonas asaccharolytica", "Porphyromonas",
  "Pint",  "Prevotella intermedia", "Prevotella"
)
target_keys <- norm_key(target_map$target_taxon)

message("[INFO] Reading significant features: ", opt$`significant-file`)
sig0 <- read_csv(opt$`significant-file`, show_col_types = FALSE, progress = FALSE)

c_tool <- pick_col(sig0, c("tool", "profiler"), label = "tool")
c_mode <- pick_col(sig0, c("spike_mode", "mode", "spike_design"), label = "spike mode")
c_label <- pick_col(sig0, c("spike_label", "target_label", "spike_taxon"), label = "spike label")
c_frac <- pick_col(sig0, c("spike_fraction", "fraction", "spike_fraction_total"), label = "spike fraction")
c_feature <- pick_col(sig0, c("feature_norm", "feature_raw", "feature", "taxon", "metadata", "member_taxon"), label = "feature")
c_q <- pick_col(sig0, c("qval", "q_value", "qvalue", "q"), label = "q value")
c_coef <- pick_col(sig0, c("coef", "coefficient", "effect", "beta"), label = "coefficient")
c_sig <- pick_col(sig0, c("is_significant", "significant"), required = FALSE)
c_pos <- pick_col(sig0, c("is_positive", "positive", "is_enriched"), required = FALSE)
c_target <- pick_col(sig0, c("is_target", "target"), required = FALSE)
c_filter <- pick_col(sig0, c("filter_mode", "filter"), required = FALSE)
c_study <- pick_col(sig0, c("background_study", "Study", "study"), required = FALSE)
c_cond <- pick_col(sig0, c("background_condition", "Target_Condition", "condition"), required = FALSE)
c_run <- pick_col(sig0, c("run_id", "analysis_id", "model_id", "run_name"), required = FALSE)

sig <- sig0 %>%
  transmute(
    tool = pretty_tool(.data[[c_tool]]),
    spike_mode = as.character(.data[[c_mode]]),
    spike_label = as.character(.data[[c_label]]),
    spike_fraction = suppressWarnings(as.numeric(.data[[c_frac]])),
    feature = as.character(.data[[c_feature]]),
    qval = suppressWarnings(as.numeric(.data[[c_q]])),
    coef = suppressWarnings(as.numeric(.data[[c_coef]])),
    is_significant_raw = if (!is.null(c_sig)) as_logical_safe(.data[[c_sig]]) else qval <= opt$`q-threshold`,
    is_positive_raw = if (!is.null(c_pos)) as_logical_safe(.data[[c_pos]]) else coef > 0,
    is_target_raw = if (!is.null(c_target)) as_logical_safe(.data[[c_target]]) else FALSE,
    filter_mode = if (!is.null(c_filter)) as.character(.data[[c_filter]]) else "original",
    study = if (!is.null(c_study)) as.character(.data[[c_study]]) else "ALL",
    condition = if (!is.null(c_cond)) as.character(.data[[c_cond]]) else "ALL",
    supplied_run_id = if (!is.null(c_run)) as.character(.data[[c_run]]) else NA_character_
  ) %>%
  filter(filter_mode == opt$`filter-mode`) %>%
  mutate(
    feature_key = norm_key(feature),
    feature_genus = extract_genus(feature),
    target_by_key = feature_key %in% target_keys,
    is_target = is_target_raw | target_by_key,
    is_fp = !is_target & is_significant_raw & is_positive_raw & !is.na(qval) & qval <= opt$`q-threshold` & coef > 0,
    context_id = if_else(
      !is.na(supplied_run_id) & nzchar(supplied_run_id), supplied_run_id,
      paste(tool, spike_mode, spike_label, spike_fraction, study, condition, filter_mode, sep = "||")
    ),
    effective_fraction = if_else(spike_mode == "community", spike_fraction / opt$`community-size`, spike_fraction),
    effective_fraction_label = fmt_fraction(effective_fraction)
  ) %>%
  left_join(target_map, by = "spike_label") %>%
  mutate(
    relation = case_when(
      spike_mode == "community" ~ "Community spike",
      feature_genus == target_genus ~ "Same genus",
      feature_genus == "Unresolved" ~ "Unresolved",
      TRUE ~ "Different genus"
    )
  )

if (!nrow(sig)) stop("No rows remained after filtering filter_mode=", opt$`filter-mode`, call. = FALSE)

# The significant-feature export contains all tested feature rows, including
# non-significant rows. The run manifest is therefore used to validate complete
# design coverage rather than to infer feature-level model contexts.
if (file.exists(opt$`run-manifest`) && file.info(opt$`run-manifest`)$size > 0) {
  manifest0 <- read_csv(opt$`run-manifest`, show_col_types = FALSE, progress = FALSE)
  manifest_required <- c("tool", "spike_mode", "spike_label", "spike_fraction")
  missing_manifest <- setdiff(manifest_required, names(manifest0))
  if (length(missing_manifest)) {
    stop("Run manifest is missing required columns: ",
         paste(missing_manifest, collapse = ", "), call. = FALSE)
  }
  expected_design <- manifest0 %>%
    transmute(
      tool = pretty_tool(tool),
      spike_mode = as.character(spike_mode),
      spike_label = as.character(spike_label),
      spike_fraction = as.numeric(spike_fraction)
    ) %>%
    distinct()
  observed_design <- sig %>%
    distinct(tool, spike_mode, spike_label, spike_fraction)
  missing_design <- anti_join(
    expected_design, observed_design,
    by = c("tool", "spike_mode", "spike_label", "spike_fraction")
  )
  if (nrow(missing_design)) {
    stop(
      "MaAsLin export is missing ", nrow(missing_design),
      " profiler/design/target/fraction combinations declared in the run manifest.",
      call. = FALSE
    )
  }
  message("[PASS] Run-manifest design coverage validated.")
} else {
  warning("Run manifest was not found; recurrence denominators will be derived from evaluated MaAsLin contexts only.")
}

# Denominators: number of model contexts per tool/design/target/fraction.
context_den <- sig %>%
  distinct(tool, spike_mode, spike_label, spike_fraction, effective_fraction,
           effective_fraction_label, study, condition, context_id) %>%
  count(tool, spike_mode, spike_label, spike_fraction, effective_fraction,
        effective_fraction_label, name = "eligible_contexts")

fp_calls <- sig %>% filter(is_fp)
if (!nrow(fp_calls)) stop("No enriched off-target calls found at q <= ", opt$`q-threshold`, call. = FALSE)

fp_recur <- fp_calls %>%
  group_by(tool, spike_mode, spike_label, spike_fraction, effective_fraction,
           effective_fraction_label, feature_key, feature, feature_genus, relation) %>%
  summarise(
    significant_contexts = n_distinct(context_id),
    median_coef = median(coef, na.rm = TRUE),
    median_q = median(qval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(context_den,
            by = c("tool", "spike_mode", "spike_label", "spike_fraction",
                   "effective_fraction", "effective_fraction_label")) %>%
  mutate(recurrence = significant_contexts / pmax(eligible_contexts, 1L))

write_csv(fp_recur, file.path(opt$outdir, "supplementary_offtarget_recurrence_by_target_fraction.csv"))

# ----------------------------- Panel A ----------------------------------
# Target x recurrent off-target heatmap for independent spikes.
heat_rank <- fp_recur %>%
  filter(spike_mode == "independent") %>%
  group_by(tool, feature_key, feature) %>%
  summarise(max_recurrence = max(recurrence, na.rm = TRUE), .groups = "drop") %>%
  group_by(tool) %>%
  slice_max(max_recurrence, n = opt$`top-features`, with_ties = FALSE) %>%
  ungroup()

heat <- fp_recur %>%
  filter(spike_mode == "independent") %>%
  semi_join(heat_rank, by = c("tool", "feature_key")) %>%
  group_by(tool, spike_label, feature_key, feature, feature_genus, relation) %>%
  summarise(recurrence = max(recurrence, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    spike_label = factor(spike_label, levels = target_map$spike_label),
    feature_short = display_feature(feature, 42)
  )

feature_order <- heat %>%
  filter(is.finite(recurrence)) %>%
  group_by(feature_short) %>%
  summarise(max_recurrence = max(recurrence), .groups = "drop") %>%
  arrange(max_recurrence, feature_short) %>%
  pull(feature_short)
feature_order <- unique(feature_order)
heat <- heat %>%
  mutate(feature_short = factor(feature_short, levels = feature_order))

missing_heat_tools <- setdiff(c("Kraken2 + Bracken", "MetaPhlAn 4"), unique(heat$tool))
heat_subtitle <- "Maximum recurrence across spike fractions; blank cells indicate that the feature was not enriched"
if (length(missing_heat_tools)) {
  heat_subtitle <- paste0(
    heat_subtitle, ". ",
    paste(missing_heat_tools, collapse = " and "),
    " produced no enriched independent-spike off-target calls"
  )
}

pA <- ggplot(heat, aes(x = spike_label, y = feature_short, fill = recurrence)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_wrap(~tool, scales = "free_y", ncol = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "white",
                      labels = percent_format(accuracy = 1),
                      limits = c(0, 1), oob = squish) +
  labs(
    title = "A. Recurrent off-target calls vary with the independently implanted taxon",
    subtitle = heat_subtitle,
    x = "Independently implanted taxon", y = "Recurrent off-target feature", fill = "Recurrence"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# ----------------------------- Panel B ----------------------------------
# Same-genus enrichment among enriched off-target versus non-enriched
# non-target feature-context observations.
# This intentionally stays within profiler because feature vocabularies are not globally harmonised.
all_eval <- sig %>%
  filter(spike_mode == "independent", !is_target) %>%
  distinct(tool, spike_label, context_id, feature_key, feature_genus, target_genus) %>%
  mutate(
    same_genus = feature_genus == target_genus & feature_genus != "Unresolved",
    enriched_offtarget = paste(tool, spike_label, context_id, feature_key, sep = "||") %in%
      paste(fp_calls$tool, fp_calls$spike_label, fp_calls$context_id,
            fp_calls$feature_key, sep = "||")
  )

fisher_one <- function(tool_value, target_value) {
  z <- all_eval %>% filter(tool == tool_value, spike_label == target_value)
  a <- sum(z$enriched_offtarget & z$same_genus, na.rm = TRUE)
  b <- sum(z$enriched_offtarget & !z$same_genus, na.rm = TRUE)
  c <- sum(!z$enriched_offtarget & z$same_genus, na.rm = TRUE)
  d <- sum(!z$enriched_offtarget & !z$same_genus, na.rm = TRUE)
  aa <- a + 0.5
  bb <- b + 0.5
  cc <- c + 0.5
  dd <- d + 0.5
  log_or <- log((aa * dd) / (bb * cc))
  se <- sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
  tibble(
    tool = tool_value, spike_label = target_value,
    odds_ratio = exp(log_or),
    conf_low = exp(log_or - 1.96 * se),
    conf_high = exp(log_or + 1.96 * se),
    n_fp = a + b, same_genus_fp = a,
    n_non_enriched = c + d, same_genus_non_enriched = c
  )
}

enrich <- bind_rows(lapply(c("Kraken2 + Bracken", "MetaPhlAn 4"), function(tt) {
  bind_rows(lapply(target_map$spike_label, function(ss) fisher_one(tt, ss)))
})) %>%
  mutate(
    spike_label = factor(spike_label, levels = rev(target_map$spike_label)),
    plot_or = if_else(n_fp > 0 & is.finite(odds_ratio) & odds_ratio > 0, odds_ratio, NA_real_),
    plot_low = if_else(n_fp > 0 & is.finite(conf_low) & conf_low > 0, conf_low, NA_real_),
    plot_high = if_else(n_fp > 0 & is.finite(conf_high) & conf_high > 0, conf_high, NA_real_)
  )
write_csv(enrich, file.path(opt$outdir, "supplementary_same_genus_enrichment.csv"))

pB_note <- enrich %>%
  group_by(tool) %>%
  summarise(n_fp_total = sum(n_fp), .groups = "drop") %>%
  filter(n_fp_total == 0) %>%
  transmute(tool, spike_label = factor("Hhat", levels = rev(target_map$spike_label)),
            x = 1, label = "No enriched off-target calls")

pB <- ggplot(enrich, aes(x = plot_or, y = spike_label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(xmin = plot_low, xmax = plot_high), width = 0.18,
                orientation = "y", na.rm = TRUE, colour = "#08519C") +
  geom_point(aes(size = n_fp), na.rm = TRUE) +
  geom_text(data = pB_note, aes(x = x, y = spike_label, label = label),
            inherit.aes = FALSE, colour = "grey35", size = 3.3) +
  facet_wrap(~tool, nrow = 1, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_x_log10(limits = c(0.1, NA)) +
  scale_size_continuous(range = c(2, 6), breaks = c(50, 100, 200, 300, 400)) +
  labs(
    title = "B. Same-genus representation differs between enriched and non-enriched observations",
    subtitle = "Enriched off-target calls are compared with non-enriched non-target feature–context observations",
    x = "Descriptive same-genus odds ratio (log scale)", y = "Implanted taxon",
    size = "Enriched off-target calls"
  ) +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold", size = 12), strip.text = element_text(face = "bold"), legend.position = "bottom")

# ----------------------------- Panel C ----------------------------------
# The existing revised Fig. 5 source already contains community-level abundance
# artefact summaries joined to off-target DA calls. Use that table directly.
message("[INFO] Reading community artefact/DA overlap: ", opt$`artifact-overlap`)
art0 <- read_csv(opt$`artifact-overlap`, show_col_types = FALSE, progress = FALSE)
required_art_cols <- c(
  "tool", "effective_fraction", "effective_fraction_label", "taxon", "taxon_key",
  "mean_abs_relative_error", "sd_relative_error", "artifact_score",
  "any_offtarget_enriched_DA", "n_offtarget_enriched_contexts",
  "min_q_offtarget_enriched", "max_coef_offtarget_enriched"
)
missing_art <- setdiff(required_art_cols, names(art0))
if (length(missing_art)) {
  stop("Artifact-overlap table is missing required columns: ",
       paste(missing_art, collapse = ", "), call. = FALSE)
}

art <- art0 %>%
  transmute(
    tool = pretty_tool(tool),
    effective_fraction = as.numeric(effective_fraction),
    effective_fraction_label = as.character(effective_fraction_label),
    feature = as.character(taxon),
    feature_key = as.character(taxon_key),
    feature_genus = extract_genus(taxon),
    mean_abs_error_pct = 100 * as.numeric(mean_abs_relative_error),
    error_sd = as.numeric(sd_relative_error),
    artifact_score = as.numeric(artifact_score),
    any_offtarget_enriched_DA = as_logical_safe(any_offtarget_enriched_DA),
    n_offtarget_enriched_contexts = as.numeric(n_offtarget_enriched_contexts),
    min_q_offtarget_enriched = as.numeric(min_q_offtarget_enriched),
    max_coef_offtarget_enriched = as.numeric(max_coef_offtarget_enriched)
  )

# Community false-positive recurrence from the MaAsLin table. Recurrence is
# calculated within profiler and effective fraction, then joined to the already
# computed abundance-artifact metrics by profiler/fraction/taxon key.
community_recur <- fp_recur %>%
  filter(spike_mode == "community", spike_label == opt$`community-label`) %>%
  group_by(tool, effective_fraction, effective_fraction_label, feature_key) %>%
  summarise(
    recurrence = max(recurrence, na.rm = TRUE),
    significant_contexts = sum(significant_contexts, na.rm = TRUE),
    .groups = "drop"
  )

scatter <- art %>%
  left_join(
    community_recur,
    by = c("tool", "effective_fraction", "effective_fraction_label", "feature_key")
  ) %>%
  mutate(
    recurrence = replace_na(recurrence, 0),
    significant_contexts = replace_na(significant_contexts, 0),
    da_status = if_else(any_offtarget_enriched_DA | recurrence > 0,
                        "Enriched off-target DA", "Not enriched")
  ) %>%
  filter(is.finite(mean_abs_error_pct), is.finite(artifact_score))

write_csv(scatter, file.path(opt$outdir, "supplementary_community_error_recurrence_joined.csv"))

# Label only the strongest recurrent/outlier taxa to keep the panel readable.
label_scatter <- scatter %>%
  filter(da_status == "Enriched off-target DA") %>%
  group_by(tool) %>%
  slice_max(order_by = artifact_score * pmax(recurrence, 0.01), n = 4, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = display_feature(feature, 30))

pC <- ggplot(scatter, aes(x = mean_abs_error_pct, y = recurrence)) +
  geom_point(aes(size = pmax(n_offtarget_enriched_contexts, significant_contexts),
                 shape = da_status, colour = da_status),
             alpha = 0.52,
             position = position_jitter(width = 0, height = 0.012, seed = 1)) +
  geom_text(data = label_scatter, aes(label = label), size = 2.2,
            check_overlap = TRUE, nudge_y = 0.025, show.legend = FALSE) +
  facet_wrap(~tool, nrow = 1, scales = "free_x") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1),
    breaks = c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 20000),
    labels = label_number(big.mark = "\u2009")
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_size_continuous(range = c(1.2, 7)) +
  scale_colour_manual(values = c("Enriched off-target DA" = "#D55E00",
                                 "Not enriched" = "grey70")) +
  scale_shape_manual(values = c("Enriched off-target DA" = 16,
                                "Not enriched" = 1)) +
  labs(
    title = "C. Community off-target recurrence is associated with abundance distortion",
    subtitle = "The pseudo-logarithmic x-axis retains zero while separating low and extreme errors",
    x = "Mean absolute relative error from dilution expectation (%)",
    y = "Off-target recurrence",
    size = "Enriched contexts", shape = "DA status", colour = "DA status"
  ) +
  theme_bw(base_size = 9.5) +
  theme(plot.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold"), legend.position = "bottom")

# ----------------------------- Panel D ----------------------------------
# UpSet-like recurrence matrix, implemented in ggplot to avoid a ComplexUpset dependency.
set_tbl <- fp_recur %>%
  group_by(tool, spike_mode, spike_label, feature_key, feature) %>%
  summarise(set_recurrence = max(recurrence, na.rm = TRUE), .groups = "drop") %>%
  mutate(set_name = if_else(spike_mode == "community", "Community", spike_label)) %>%
  filter(set_recurrence >= opt$`min-set-recurrence`) %>%
  distinct(tool, set_name, feature_key, feature)

make_upset <- function(tool_name) {
  x <- set_tbl %>% filter(tool == tool_name)
  if (!nrow(x)) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = paste("No recurrent intersections for", tool_name)) + theme_void())
  }
  set_order <- c(target_map$spike_label, "Community")
  present_sets <- set_order
  wide <- x %>%
    mutate(value = 1L, set_name = factor(set_name, levels = present_sets)) %>%
    select(feature_key, feature, set_name, value) %>%
    distinct() %>%
    pivot_wider(names_from = set_name, values_from = value, values_fill = 0L)
  signature_cols <- intersect(present_sets, names(wide))
  if (!length(signature_cols)) return(ggplot() + theme_void())
  wide <- wide %>%
    rowwise() %>%
    mutate(signature = paste(signature_cols[which(c_across(all_of(signature_cols)) == 1L)], collapse = " & ")) %>%
    ungroup() %>% filter(nzchar(signature))
  ints <- wide %>% count(signature, sort = TRUE, name = "intersection_size") %>%
    slice_head(n = opt$`top-intersections`) %>% mutate(intersection_id = row_number())
  # Rebuild membership safely without vector recycling.
  long_members <- ints %>% separate_rows(signature, sep = " & ") %>%
    transmute(intersection_id, set_name = signature, member = 1L)
  membership <- expand_grid(intersection_id = ints$intersection_id, set_name = present_sets) %>%
    left_join(long_members, by = c("intersection_id", "set_name")) %>%
    mutate(
      member = factor(replace_na(member, 0L), levels = c(0, 1)),
      set_name = factor(set_name, levels = rev(present_sets))
    )

  profiler_colour <- if (tool_name == "Kraken2 + Bracken") "#009E73" else "#6C5CE7"
  bars <- ggplot(ints, aes(x = factor(intersection_id), y = intersection_size)) +
    geom_col(width = 0.72, fill = profiler_colour) +
    geom_text(aes(label = intersection_size), vjust = -0.25, size = 2.7) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    labs(title = tool_name, x = NULL, y = "Features") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.x = element_blank())
  matrix_plot <- ggplot(membership, aes(x = factor(intersection_id), y = set_name)) +
    geom_point(aes(alpha = member), size = 2.8) +
    geom_line(data = membership %>% filter(member == 1L) %>% group_by(intersection_id) %>%
                filter(n() > 1) %>% arrange(set_name), aes(group = intersection_id), linewidth = 0.35) +
    scale_alpha_manual(values = c(`0` = 0.12, `1` = 1), guide = "none") +
    labs(x = "Top recurrent intersections", y = NULL) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  bars / matrix_plot + plot_layout(heights = c(0.58, 0.42))
}

pD_body <- make_upset("Kraken2 + Bracken") | make_upset("MetaPhlAn 4")
pD <- pD_body +
  plot_annotation(
    title = "D. Recurrent off-target calls shared by individual and community spike designs",
    subtitle = paste0("Membership requires recurrence in at least ", percent(opt$`min-set-recurrence`, accuracy = 1), " of eligible contexts"),
    theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9))
  )

# -------------------------- Focused case studies -------------------------
focus_targets <- str_split(opt$`focus-targets`, ",", simplify = TRUE) %>% as.character() %>% trimws()
focus <- fp_recur %>%
  filter(spike_mode == "independent", spike_label %in% focus_targets) %>%
  group_by(tool, spike_label, effective_fraction, effective_fraction_label,
           feature_key, feature, feature_genus, relation) %>%
  summarise(recurrence = max(recurrence), median_coef = median(median_coef), .groups = "drop") %>%
  group_by(tool, spike_label, feature_key, feature) %>%
  mutate(max_rec = max(recurrence)) %>%
  ungroup() %>%
  group_by(tool, spike_label) %>%
  filter(dense_rank(desc(max_rec)) <= 12) %>%
  ungroup() %>%
  mutate(feature_short = display_feature(feature, 38))

pFocus <- ggplot(focus, aes(x = effective_fraction_label,
                            y = fct_reorder(feature_short, max_rec, .fun = max),
                            size = recurrence, fill = median_coef)) +
  geom_point(shape = 21, colour = "grey25", stroke = 0.25) +
  facet_grid(spike_label ~ tool, scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(1.2, 7), labels = percent_format(accuracy = 1)) +
  scale_fill_gradient2(low = "grey80", mid = "white", high = "black", midpoint = 0) +
  labs(
    title = "Focused target-specific off-target recurrence",
    subtitle = "Point size shows false-positive recurrence; fill shows median positive MaAsLin coefficient",
    x = "Effective per-taxon spike fraction", y = "Off-target feature",
    size = "Recurrence", fill = "Median coefficient"
  ) +
  theme_bw(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold", size = 12), strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 6.5), panel.grid.minor = element_blank(), legend.position = "bottom")

# ---------------------------- Save outputs -------------------------------
save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(opt$outdir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(opt$outdir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 400, bg = "white")
}

save_plot(pA, "panel_A_target_by_offtarget_recurrence_heatmap", 10.5, 10.0)
save_plot(pB, "panel_B_same_genus_enrichment", 10.5, 5.8)
save_plot(pC, "panel_C_error_vs_false_positive_recurrence", 10.5, 5.8)
save_plot(pD, "panel_D_individual_community_upset", 12.0, 7.4)
save_plot(pFocus, "supplementary_focused_target_offtarget_bubble_heatmaps", 12.0, 10.5)

composite <- pA / pB / pC / wrap_elements(full = pD) +
  plot_layout(heights = c(1.35, 0.72, 0.78, 1.10)) +
  plot_annotation(
    title = "Taxon-specific structure and recurrence of off-target differential-abundance calls",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0))
  )
save_plot(composite, "supplementary_taxon_specific_offtarget_patterns", opt$width, opt$height)

manifest <- c(
  "Supplementary taxon-specific off-target analysis",
  paste("Significant feature file:", normalizePath(opt$`significant-file`)),
  paste("Artifact-overlap file:", normalizePath(opt$`artifact-overlap`)),
  paste("q threshold:", opt$`q-threshold`),
  paste("filter mode:", opt$`filter-mode`),
  paste("minimum UpSet recurrence:", opt$`min-set-recurrence`),
  "",
  "Main outputs:",
  "- supplementary_taxon_specific_offtarget_patterns.pdf/png",
  "- supplementary_focused_target_offtarget_bubble_heatmaps.pdf/png",
  "- supplementary_offtarget_recurrence_by_target_fraction.csv",
  "- supplementary_same_genus_enrichment.csv",
  "- supplementary_community_error_recurrence_joined.csv",
  "",
  "Interpretation notes:",
  "- Feature overlap is evaluated within profiler because non-target vocabularies are not globally harmonised.",
  "- Same-genus enrichment is a descriptive/probing analysis and does not prove read-level misassignment.",
  "- Community set membership means recurrent in the mixed-community analysis, not attributable to one member."
)
writeLines(manifest, file.path(opt$outdir, "README_outputs.txt"))
message("[OK] Wrote supplementary off-target analysis under: ", normalizePath(opt$outdir))
