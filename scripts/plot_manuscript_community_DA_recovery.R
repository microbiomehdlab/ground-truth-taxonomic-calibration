#!/usr/bin/env Rscript
# Polished v4: fixed manuscript-wide taxon order and improved label legibility.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && !all(is.na(x)) && !(is.character(x) && !nzchar(x[1]))) x else y
}

manuscript_taxon_order <- c(
  "Bfrag", "Csym", "Dpne", "Fnuc", "Hhat",
  "Pmic", "Pana", "Psto", "Porp", "Pint"
)

# optparse sometimes exposes hyphenated option names awkwardly. Define the two
# custom fraction flags with underscores and normalize hyphenated CLI flags below.
option_list <- list(
  make_option("--indir", type = "character", default = "RUNS/maaslin_spike",
              help = "Input directory containing MaAsLin summary CSVs [default %default]"),
  make_option("--member-file", type = "character", default = NULL,
              help = "Optional direct path to maaslin_member_detection_ALLFILTERS.csv"),
  make_option("--run-summary-file", type = "character", default = NULL,
              help = "Optional direct path to maaslin_run_summary_ALLFILTERS.csv"),
  make_option("--outdir", type = "character", default = "RUNS/plots_manuscript_community_da_recovery",
              help = "Output directory [default %default]"),
  make_option("--community-size", type = "integer", default = 10,
              help = "Number of taxa in the community spike. Effective per-species fraction = total community fraction / this value [default %default]"),
  make_option("--matched_effective_fractions", type = "character", default = "0.0001,0.0005,0.001",
              help = "Effective per-species fractions used for independent-vs-community DA comparison [default %default = 0.01%%,0.05%%,0.10%%]"),
  make_option("--main_effective_fractions", type = "character", default = "0.00001,0.00005,0.0001,0.0005,0.001",
              help = "Effective per-species fractions shown in main community-only Panels B/C [default %default = 0.001%%,0.005%%,0.01%%,0.05%%,0.10%%]"),
  make_option("--filter-mode", type = "character", default = "original",
              help = "Filter mode to use if filter_mode column exists [default %default]"),
  make_option("--spike-label-order", type = "character", default = paste(manuscript_taxon_order, collapse = ","),
              help = "Comma-separated manuscript target order, displayed top-to-bottom [default %default]"),
  make_option("--main-width", type = "double", default = 15.2,
              help = "Main figure width in inches [default %default]"),
  make_option("--main-height", type = "double", default = 18.0,
              help = "Main figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 320,
              help = "PNG resolution [default %default]")
)

raw_args <- commandArgs(trailingOnly = TRUE)
raw_args <- sub("^--matched-effective-fractions$", "--matched_effective_fractions", raw_args)
raw_args <- sub("^--main-effective-fractions$", "--main_effective_fractions", raw_args)
opt <- parse_args(OptionParser(option_list = option_list), args = raw_args)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

member_file <- opt$`member-file` %||% file.path(opt$indir, "maaslin_member_detection_ALLFILTERS.csv")
run_summary_file <- opt$`run-summary-file` %||% file.path(opt$indir, "maaslin_run_summary_ALLFILTERS.csv")

if (!file.exists(member_file)) stop("Could not find member detection file: ", member_file, call. = FALSE)
if (!file.exists(run_summary_file)) stop("Could not find run summary file: ", run_summary_file, call. = FALSE)

parse_num_vector <- function(x) {
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)] %>% as.numeric()
}
parse_chr_vector <- function(x) {
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)]
}
fmt_fraction <- function(x) {
  p <- as.numeric(x) * 100
  ifelse(p < 0.01, sprintf("%.3f%%", p), sprintf("%.2f%%", p))
}
to_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  lx <- tolower(trimws(as.character(x)))
  lx %in% c("true", "t", "1", "yes", "y")
}

matched_effective_fractions <- parse_num_vector(opt$matched_effective_fractions)
main_effective_fractions <- parse_num_vector(opt$main_effective_fractions)
preferred_label_order <- parse_chr_vector(opt$`spike-label-order`)
if (!identical(preferred_label_order, manuscript_taxon_order)) {
  stop(
    "--spike-label-order must match the manuscript-wide order exactly: ",
    paste(manuscript_taxon_order, collapse = ","),
    call. = FALSE
  )
}

tool_labeller <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)
tool_levels <- unname(tool_labeller)

profiler_cols <- c(
  "Kraken2 + Bracken" = "#2A9D8F",
  "MetaPhlAn 4" = "#7B61D1"
)
mode_cols <- c(
  "Independent" = "#4C78A8",
  "Community" = "#D95F02"
)
condition_cols <- c(
  "Adenoma" = "#D8A03A",
  "colorectal carcinoma" = "#D95F02",
  "Control" = "#4C78A8"
)

fraction_cols <- c(
  "0.001%" = "#244D8F",
  "0.005%" = "#3B6EA8",
  "0.01%"  = "#2A9D8F",
  "0.05%"  = "#63B179",
  "0.10%"  = "#D99B2B",
  "0.50%"  = "#B65B84",
  "1.00%"  = "#666666"
)

pub_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.28),
      panel.grid.major.x = element_line(colour = "#F3F3F3", linewidth = 0.22),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "#4D4D4D", linewidth = 0.34),
      axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.30),
      axis.ticks.length = unit(2.0, "pt"),
      strip.background = element_rect(fill = "#F8F8F8", colour = "#D9D9D9", linewidth = 0.42),
      strip.text = element_text(face = "bold", colour = "#222222", margin = margin(3, 3, 3, 3)),
      axis.title = element_text(face = "bold", colour = "#222222"),
      axis.text = element_text(colour = "#444444"),
      legend.title = element_text(face = "bold", colour = "#222222", size = rel(0.92)),
      legend.text = element_text(colour = "#333333", size = rel(0.86)),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "bold", size = base_size + 2.4, hjust = 0, colour = "#111111"),
      plot.subtitle = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )
}

save_plot_set <- function(plot_obj, stem, width, height, dpi = 320) {
  pdf_path <- file.path(opt$outdir, paste0(stem, ".pdf"))
  png_path <- file.path(opt$outdir, paste0(stem, ".png"))
  if (capabilities("cairo")) {
    ggsave(pdf_path, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  } else {
    ggsave(pdf_path, plot_obj, width = width, height = height, bg = "white")
  }
  ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

message("[INFO] Reading member detection: ", member_file)
member_raw <- read_csv(member_file, show_col_types = FALSE)
message("[INFO] Reading run summary: ", run_summary_file)
run_raw <- read_csv(run_summary_file, show_col_types = FALSE)

member_required <- c(
  "background_condition", "background_study", "tool", "spike_mode", "spike_label",
  "spike_fraction", "member_taxon", "member_detected_positive",
  "min_q_target", "max_coef_target"
)
missing_member <- setdiff(member_required, names(member_raw))
if (length(missing_member) > 0) {
  stop("Member detection file is missing required columns: ", paste(missing_member, collapse = ", "), call. = FALSE)
}

run_required <- c(
  "background_condition", "background_study", "tool", "spike_mode", "spike_label",
  "spike_fraction", "n_sig_offtarget_enriched", "n_sig_offtarget",
  "n_sig_target_enriched"
)
missing_run <- setdiff(run_required, names(run_raw))
if (length(missing_run) > 0) {
  stop("Run summary file is missing required columns: ", paste(missing_run, collapse = ", "), call. = FALSE)
}

# Optional filter-mode selection.
if ("filter_mode" %in% names(member_raw)) {
  member_raw <- member_raw %>% filter(filter_mode == opt$`filter-mode`)
}
if ("filter_mode" %in% names(run_raw)) {
  run_raw <- run_raw %>% filter(filter_mode == opt$`filter-mode`)
}

# Map community target members onto the independent-spike short labels.
label_map <- member_raw %>%
  filter(spike_mode == "independent", !is.na(spike_label), spike_label != "CRCpanel") %>%
  distinct(member_taxon, spike_label) %>%
  group_by(member_taxon) %>%
  summarise(spike_label_short = first(spike_label), .groups = "drop")

member_tbl <- member_raw %>%
  filter(spike_mode %in% c("independent", "community")) %>%
  mutate(
    tool = recode(as.character(tool), !!!tool_labeller, .default = as.character(tool)),
    tool = factor(tool, levels = tool_levels),
    spike_fraction = as.numeric(spike_fraction),
    effective_fraction = if_else(
      spike_mode == "community",
      spike_fraction / opt$`community-size`,
      spike_fraction
    ),
    effective_fraction_label = fmt_fraction(effective_fraction),
    member_detected_positive = to_logical_safe(member_detected_positive),
    member_detected_any = if ("member_detected_any" %in% names(.)) to_logical_safe(member_detected_any) else NA,
    min_q_target = as.numeric(min_q_target),
    max_coef_target = as.numeric(max_coef_target),
    neglog10_q = -log10(pmax(min_q_target, 1e-300))
  ) %>%
  left_join(label_map, by = "member_taxon") %>%
  mutate(
    spike_label_clean = case_when(
      spike_mode == "independent" ~ as.character(spike_label),
      spike_mode == "community" ~ spike_label_short,
      TRUE ~ as.character(spike_label)
    )
  ) %>%
  filter(!is.na(spike_label_clean), spike_label_clean != "CRCpanel")

run_tbl <- run_raw %>%
  filter(spike_mode %in% c("independent", "community")) %>%
  mutate(
    tool = recode(as.character(tool), !!!tool_labeller, .default = as.character(tool)),
    tool = factor(tool, levels = tool_levels),
    spike_fraction = as.numeric(spike_fraction),
    effective_fraction = if_else(
      spike_mode == "community",
      spike_fraction / opt$`community-size`,
      spike_fraction
    ),
    effective_fraction_label = fmt_fraction(effective_fraction),
    n_sig_offtarget_enriched = as.numeric(n_sig_offtarget_enriched),
    n_sig_offtarget = as.numeric(n_sig_offtarget),
    n_sig_target_enriched = as.numeric(n_sig_target_enriched)
  )

if (!nrow(member_tbl)) stop("No member-detection rows remained after filtering.", call. = FALSE)
if (!nrow(run_tbl)) stop("No run-summary rows remained after filtering.", call. = FALSE)

message("[INFO] Member rows after filtering: ", nrow(member_tbl))
message("[INFO] Run-summary rows after filtering: ", nrow(run_tbl))
message("[INFO] Matched effective fractions for Panel A: ", paste(fmt_fraction(matched_effective_fractions), collapse = ", "))
message("[INFO] Main community-only effective fractions for Panels B/C: ", paste(fmt_fraction(main_effective_fractions), collapse = ", "))

# ----------------------------
# Summary tables
# ----------------------------
summarise_da <- function(dat, by_context = FALSE) {
  group_cols <- c("spike_mode", "tool", "spike_label_clean", "effective_fraction", "effective_fraction_label")
  if (by_context) group_cols <- c(group_cols, "background_condition", "background_study")
  
  dat %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_tests = n(),
      positive_da_rate = mean(member_detected_positive, na.rm = TRUE),
      any_da_rate = if (all(is.na(member_detected_any))) NA_real_ else mean(member_detected_any, na.rm = TRUE),
      median_neglog10_q = median(neglog10_q, na.rm = TRUE),
      median_positive_coef = median(max_coef_target, na.rm = TRUE),
      .groups = "drop"
    )
}

da_overall <- summarise_da(member_tbl, by_context = FALSE)
da_context <- summarise_da(member_tbl, by_context = TRUE)

write_csv(da_overall, file.path(opt$outdir, "da_member_detection_summary_overall.csv"))
write_csv(da_context, file.path(opt$outdir, "da_member_detection_summary_by_context.csv"))

# Shared manuscript order. ggplot draws the first y-factor level at the bottom,
# so reverse once here to display Bfrag -> Pint from top to bottom.
missing_manuscript_taxa <- setdiff(manuscript_taxon_order, unique(member_tbl$spike_label_clean))
if (length(missing_manuscript_taxa)) {
  stop(
    "Input is missing manuscript taxa required for Figure 5: ",
    paste(missing_manuscript_taxa, collapse = ", "),
    call. = FALSE
  )
}
display_levels <- rev(manuscript_taxon_order)

# ----------------------------
# Panel A: independent versus community DA recovery at matched effective fractions
# ----------------------------
matched_labels <- fmt_fraction(matched_effective_fractions)

mode_rate_tbl <- da_overall %>%
  filter(spike_mode %in% c("independent", "community"),
         effective_fraction %in% matched_effective_fractions) %>%
  mutate(
    spike_mode_pretty = recode(
      spike_mode,
      independent = "Independent",
      community = "Community"
    ),
    effective_fraction_label = factor(effective_fraction_label, levels = matched_labels),
    spike_label_plot = factor(spike_label_clean, levels = display_levels)
  ) %>%
  group_by(tool, spike_label_clean, spike_label_plot, effective_fraction, effective_fraction_label, spike_mode_pretty) %>%
  summarise(positive_da_rate = mean(positive_da_rate, na.rm = TRUE), .groups = "drop")

mode_pair_tbl <- mode_rate_tbl %>%
  tidyr::pivot_wider(
    names_from = spike_mode_pretty,
    values_from = positive_da_rate,
    values_fill = 0
  ) %>%
  mutate(spike_label_plot = factor(as.character(spike_label_plot), levels = display_levels))

write_csv(mode_rate_tbl, file.path(opt$outdir, "panel_A_mode_target_DA_detection_matched_effective_fractions.csv"))

panel_A_delta_summary <- mode_pair_tbl %>%
  mutate(delta_community_minus_independent = Community - Independent) %>%
  group_by(tool, effective_fraction_label) %>%
  summarise(
    n_targets = n(),
    mean_independent = mean(Independent, na.rm = TRUE),
    mean_community = mean(Community, na.rm = TRUE),
    median_delta = median(delta_community_minus_independent, na.rm = TRUE),
    prop_community_higher = mean(delta_community_minus_independent > 0, na.rm = TRUE),
    prop_equal = mean(abs(delta_community_minus_independent) < 1e-12, na.rm = TRUE),
    prop_community_lower = mean(delta_community_minus_independent < 0, na.rm = TRUE),
    spearman_rho = suppressWarnings(cor(Independent, Community, method = "spearman", use = "complete.obs")),
    .groups = "drop"
  )

write_csv(panel_A_delta_summary, file.path(opt$outdir, "panel_A_independent_vs_community_target_DA_detection_delta_summary.csv"))

pA <- ggplot(mode_pair_tbl, aes(y = spike_label_plot)) +
  geom_segment(
    aes(x = Independent, xend = Community, yend = spike_label_plot),
    colour = "#C4C4C4",
    linewidth = 0.52,
    lineend = "round",
    na.rm = TRUE
  ) +
  geom_point(aes(x = Independent, colour = "Independent"), size = 3.0, alpha = 0.96, na.rm = TRUE) +
  geom_point(aes(x = Community, colour = "Community"), size = 3.0, alpha = 0.96, na.rm = TRUE) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#9B9B9B", linewidth = 0.40) +
  facet_grid(tool ~ effective_fraction_label) +
  scale_x_continuous(
    # Expand slightly beyond the biological 0–100% range so endpoint points are fully visible.
    limits = c(-0.035, 1.035),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(values = mode_cols, name = "Spike design") +
  labs(
    title = "A. Community spikes shift target DA detection upward at matched fractions",
    x = "Target significant and enriched (% runs)",
    y = NULL
  ) +
  pub_theme(10.7) +
  theme(
    axis.text.y = element_text(size = 9.6, face = "bold"),
    axis.text.x = element_text(size = 9.1),
    strip.text = element_text(size = 9.8, face = "bold", margin = margin(4, 4, 4, 4)),
    panel.grid.major.y = element_line(colour = "#F1F1F1", linewidth = 0.24),
    panel.grid.major.x = element_line(colour = "#F4F4F4", linewidth = 0.22),
    panel.spacing.x = unit(0.95, "lines"),
    panel.spacing.y = unit(0.68, "lines"),
    legend.position = "right",
    plot.margin = margin(4, 8, 4, 4)
  )

# Background/study-resolved concordance scatter for supplement.
ind_ctx <- da_context %>%
  filter(spike_mode == "independent", effective_fraction %in% matched_effective_fractions) %>%
  select(tool, spike_label_clean, effective_fraction, effective_fraction_label,
         background_condition, background_study, ind_da_rate = positive_da_rate)

comm_ctx <- da_context %>%
  filter(spike_mode == "community", effective_fraction %in% matched_effective_fractions) %>%
  select(tool, spike_label_clean, effective_fraction, effective_fraction_label,
         background_condition, background_study, comm_da_rate = positive_da_rate)

concordance_tbl <- inner_join(
  ind_ctx, comm_ctx,
  by = c("tool", "spike_label_clean", "effective_fraction", "effective_fraction_label",
         "background_condition", "background_study")
) %>%
  mutate(
    effective_fraction_label = factor(effective_fraction_label, levels = matched_labels),
    background_condition = factor(background_condition, levels = c("Adenoma", "colorectal carcinoma", "Control")),
    Study_short = recode(background_study, FengQ_2015 = "FengQ", ZellerG_2014 = "Zeller", .default = background_study)
  )

rho_tbl <- concordance_tbl %>%
  group_by(tool, effective_fraction_label) %>%
  summarise(
    n = n(),
    rho = suppressWarnings(cor(ind_da_rate, comm_da_rate, method = "spearman", use = "complete.obs")),
    median_delta = median(comm_da_rate - ind_da_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(is.finite(rho),
                   paste0("rho=", sprintf("%.2f", rho), "\nΔ=", percent(median_delta, accuracy = 1)),
                   paste0("n=", n))
  )

write_csv(concordance_tbl, file.path(opt$outdir, "supp_panel_A_DA_concordance_by_context.csv"))
write_csv(rho_tbl, file.path(opt$outdir, "supp_panel_A_DA_concordance_spearman_summary.csv"))

pA_concordance_supp <- ggplot(concordance_tbl, aes(x = ind_da_rate, y = comm_da_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8A8A8A", linewidth = 0.45) +
  geom_point(
    aes(fill = background_condition, shape = Study_short),
    size = 2.45,
    stroke = 0.35,
    colour = "#333333",
    alpha = 0.85,
    na.rm = TRUE
  ) +
  geom_text(
    data = rho_tbl,
    aes(x = 0.04, y = 0.96, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    colour = "#333333"
  ) +
  facet_grid(tool ~ effective_fraction_label) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = condition_cols, name = "Background") +
  scale_shape_manual(values = c("FengQ" = 21, "Zeller" = 24), name = "Study") +
  labs(
    title = "Supplementary. Background-resolved comparison of target DA detection",
    x = "Independent spikes: target significant and enriched (% runs)\nbalanced subset",
    y = "Community spikes: target significant and enriched (% runs)\nfull cohort"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(colour = "#EEEEEE", linewidth = 0.25),
    strip.text = element_text(size = 9.2, face = "bold"),
    panel.spacing = unit(0.62, "lines"),
    plot.margin = margin(4, 8, 4, 4)
  )

# ----------------------------
# Panel B: full-cohort community DA recovery across ultra-low and weak fractions
# ----------------------------
community_overall_all <- da_overall %>%
  filter(spike_mode == "community") %>%
  mutate(
    spike_label_plot = factor(spike_label_clean, levels = display_levels),
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(sort(unique(effective_fraction)))
    )
  )

community_overall <- community_overall_all %>%
  filter(effective_fraction %in% main_effective_fractions) %>%
  mutate(
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(main_effective_fractions)
    )
  )

if (!nrow(community_overall)) {
  stop("No community DA rows matched --main-effective-fractions: ",
       paste(main_effective_fractions, collapse = ", "), call. = FALSE)
}

make_tool_dumbbell <- function(dat, rate_col, title_text, x_label) {
  pair_tbl <- dat %>%
    mutate(
      tool_short = case_when(
        tool == "Kraken2 + Bracken" ~ "kraken",
        tool == "MetaPhlAn 4" ~ "metaphlan",
        TRUE ~ as.character(tool)
      )
    ) %>%
    group_by(spike_label_plot, effective_fraction_label, tool_short) %>%
    summarise(rate_value = mean(.data[[rate_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = tool_short,
      values_from = rate_value,
      values_fill = 0
    ) %>%
    mutate(spike_label_plot = factor(as.character(spike_label_plot), levels = display_levels))
  
  ggplot(pair_tbl, aes(y = spike_label_plot)) +
    geom_segment(
      aes(x = kraken, xend = metaphlan, yend = spike_label_plot),
      colour = "#C4C4C4",
      linewidth = 0.50,
      lineend = "round"
    ) +
    geom_point(aes(x = kraken, colour = "Kraken2 + Bracken"), size = 2.85, alpha = 0.96) +
    geom_point(aes(x = metaphlan, colour = "MetaPhlAn 4"), size = 2.85, alpha = 0.96) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "#9B9B9B", linewidth = 0.38) +
    facet_wrap(~ effective_fraction_label, nrow = 1) +
    scale_x_continuous(
      limits = c(-0.035, 1.035),
      breaks = c(0, 0.5, 1),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_colour_manual(values = profiler_cols, name = "Profiler") +
    labs(title = title_text, x = x_label, y = NULL) +
    pub_theme(10.6) +
    theme(
      axis.text.y = element_text(size = 9.4, face = "bold"),
      axis.text.x = element_text(size = 8.9),
      strip.text = element_text(size = 9.3, face = "bold", margin = margin(4, 3, 4, 3)),
      panel.grid.major.y = element_line(colour = "#F1F1F1", linewidth = 0.24),
      panel.grid.major.x = element_line(colour = "#F4F4F4", linewidth = 0.22),
      panel.spacing.x = unit(1.05, "lines"),
      legend.position = "right",
      plot.margin = margin(4, 8, 4, 4)
    )
}

pB <- make_tool_dumbbell(
  community_overall,
  rate_col = "positive_da_rate",
  title_text = "C. Full-cohort community target DA detection across ultra-low and weak effective fractions",
  x_label = "Target significant and enriched (% runs)"
)

pB_all <- make_tool_dumbbell(
  community_overall_all,
  rate_col = "positive_da_rate",
  title_text = "Supplementary. Full-cohort community target DA detection across all effective fractions",
  x_label = "Target significant and enriched (% runs)"
)

# Supplementary q-value intensity heatmap.
q_heat <- community_overall_all %>%
  mutate(
    spike_label_plot = factor(spike_label_clean, levels = display_levels),
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(sort(unique(effective_fraction)))
    )
  )

pB_q_supp <- ggplot(q_heat, aes(x = effective_fraction_label, y = spike_label_plot, fill = median_neglog10_q)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  facet_wrap(~ tool, ncol = 1) +
  scale_fill_gradient(
    low = "#F1F1F1",
    high = "#2166AC",
    name = expression(Median~ -log[10](q))
  ) +
  labs(
    title = "Supplementary. Median differential-abundance significance across community spike fractions",
    x = "Effective per-species fraction",
    y = NULL
  ) +
  pub_theme(10.6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# ----------------------------
# Panel C: study-specific community target DA detection
# ----------------------------
community_by_study <- member_tbl %>%
  filter(spike_mode == "community",
         effective_fraction %in% main_effective_fractions) %>%
  mutate(
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(main_effective_fractions)
    ),
    Study_short = recode(background_study, FengQ_2015 = "FengQ", ZellerG_2014 = "Zeller", .default = background_study)
  ) %>%
  group_by(tool, Study_short, effective_fraction, effective_fraction_label, spike_label_clean) %>%
  summarise(
    target_da_detection = mean(member_detected_positive, na.rm = TRUE),
    .groups = "drop"
  )

study_summary <- community_by_study %>%
  group_by(tool, Study_short, effective_fraction, effective_fraction_label) %>%
  summarise(
    n_targets = n(),
    median_target_detection = median(target_da_detection, na.rm = TRUE),
    q25_target_detection = quantile(target_da_detection, 0.25, na.rm = TRUE, names = FALSE),
    q75_target_detection = quantile(target_da_detection, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  mutate(
    Study_short = factor(Study_short, levels = c("FengQ", "Zeller")),
    tool = factor(tool, levels = tool_levels)
  )

write_csv(community_by_study, file.path(opt$outdir, "panel_C_study_target_DA_detection_by_species_main_effective_fractions.csv"))
write_csv(study_summary, file.path(opt$outdir, "panel_C_study_target_DA_detection_summary_main_effective_fractions.csv"))

study_cols <- c(
  "FengQ" = "#4C78A8",
  "Zeller" = "#D8A03A"
)

pC <- ggplot(study_summary, aes(x = effective_fraction_label, y = median_target_detection,
                                color = Study_short, group = Study_short)) +
  geom_line(linewidth = 0.75, alpha = 0.92) +
  geom_errorbar(
    aes(ymin = q25_target_detection, ymax = q75_target_detection),
    width = 0.10,
    linewidth = 0.45,
    alpha = 0.72
  ) +
  geom_point(size = 2.6, alpha = 0.96) +
  facet_grid(tool ~ .) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_color_manual(values = study_cols, name = "Study") +
  labs(
    title = "C. Study-specific community target DA detection",
    x = "Effective per-species fraction",
    y = "Target significant and enriched\nmedian and IQR across species"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 8.7),
    strip.text.y = element_text(angle = 270, size = 9.0, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#EFEFEF", linewidth = 0.25),
    plot.margin = margin(4, 8, 4, 4)
  )

# Supplementary study-resolved species-level panel.
pC_species_supp <- ggplot(community_by_study, aes(x = target_da_detection, y = factor(spike_label_clean, levels = display_levels))) +
  geom_segment(
    data = community_by_study %>%
      tidyr::pivot_wider(names_from = Study_short, values_from = target_da_detection) %>%
      filter(!is.na(FengQ), !is.na(Zeller)),
    aes(
      x = FengQ,
      xend = Zeller,
      y = factor(spike_label_clean, levels = display_levels),
      yend = factor(spike_label_clean, levels = display_levels)
    ),
    inherit.aes = FALSE,
    colour = "#C4C4C4",
    linewidth = 0.45,
    lineend = "round"
  ) +
  geom_point(aes(color = Study_short), size = 2.3, alpha = 0.95) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#9B9B9B", linewidth = 0.38) +
  facet_grid(tool ~ effective_fraction_label) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_color_manual(values = study_cols, name = "Study") +
  labs(
    title = "Supplementary. Study-resolved community target DA detection by species",
    x = "Target significant and enriched (% runs)",
    y = NULL
  ) +
  pub_theme(10.5) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8.3),
    axis.text.x = element_text(size = 8.2),
    strip.text = element_text(size = 8.4, face = "bold"),
    panel.grid.major.y = element_line(colour = "#F1F1F1", linewidth = 0.24),
    panel.grid.major.x = element_line(colour = "#F4F4F4", linewidth = 0.22),
    panel.spacing.x = unit(0.75, "lines")
  )

# ----------------------------
# Panel D: off-target DA burden
# ----------------------------
run_community_all <- run_tbl %>%
  filter(spike_mode == "community") %>%
  mutate(
    effective_fraction_label = factor(
      effective_fraction_label,
      levels = fmt_fraction(sort(unique(effective_fraction)))
    ),
    background_condition = factor(background_condition, levels = c("Adenoma", "colorectal carcinoma", "Control")),
    Study_short = recode(background_study, FengQ_2015 = "FengQ", ZellerG_2014 = "Zeller", .default = background_study)
  )

run_community_main <- run_community_all %>%
  filter(effective_fraction %in% main_effective_fractions) %>%
  mutate(effective_fraction_label = factor(effective_fraction_label, levels = fmt_fraction(main_effective_fractions)))

offtarget_summary <- run_community_main %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  summarise(
    n_runs = n(),
    median_offtarget_enriched = median(n_sig_offtarget_enriched, na.rm = TRUE),
    q25_offtarget_enriched = quantile(n_sig_offtarget_enriched, 0.25, na.rm = TRUE, names = FALSE),
    q75_offtarget_enriched = quantile(n_sig_offtarget_enriched, 0.75, na.rm = TRUE, names = FALSE),
    median_offtarget_total = median(n_sig_offtarget, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(run_community_main, file.path(opt$outdir, "panel_D_offtarget_run_level_main_effective_fractions.csv"))
write_csv(offtarget_summary, file.path(opt$outdir, "panel_D_offtarget_summary_main_effective_fractions.csv"))

pD <- ggplot(offtarget_summary, aes(x = effective_fraction_label, y = median_offtarget_enriched, color = tool, group = tool)) +
  geom_line(linewidth = 0.75, alpha = 0.92) +
  geom_errorbar(
    aes(ymin = q25_offtarget_enriched, ymax = q75_offtarget_enriched),
    width = 0.10,
    linewidth = 0.45,
    alpha = 0.72
  ) +
  geom_point(size = 2.6, alpha = 0.96) +
  scale_color_manual(values = profiler_cols, name = "Profiler") +
  labs(
    title = "D. Off-target enriched DA burden in community spike-ins",
    x = "Effective per-species fraction",
    y = "Off-target enriched features\nmedian and IQR across runs"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 8.7),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#EFEFEF", linewidth = 0.25),
    plot.margin = margin(4, 8, 4, 4)
  )

offtarget_summary_all <- run_community_all %>%
  group_by(tool, effective_fraction, effective_fraction_label) %>%
  summarise(
    n_runs = n(),
    median_offtarget_enriched = median(n_sig_offtarget_enriched, na.rm = TRUE),
    q25_offtarget_enriched = quantile(n_sig_offtarget_enriched, 0.25, na.rm = TRUE, names = FALSE),
    q75_offtarget_enriched = quantile(n_sig_offtarget_enriched, 0.75, na.rm = TRUE, names = FALSE),
    median_offtarget_total = median(n_sig_offtarget, na.rm = TRUE),
    q25_offtarget_total = quantile(n_sig_offtarget, 0.25, na.rm = TRUE, names = FALSE),
    q75_offtarget_total = quantile(n_sig_offtarget, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

write_csv(offtarget_summary_all, file.path(opt$outdir, "supp_offtarget_summary_all_effective_fractions.csv"))

pD_all <- ggplot(offtarget_summary_all, aes(x = effective_fraction_label, y = median_offtarget_enriched, color = tool, group = tool)) +
  geom_line(linewidth = 0.75, alpha = 0.92) +
  geom_errorbar(aes(ymin = q25_offtarget_enriched, ymax = q75_offtarget_enriched),
                width = 0.10, linewidth = 0.45, alpha = 0.72) +
  geom_point(size = 2.6, alpha = 0.96) +
  scale_color_manual(values = profiler_cols, name = "Profiler") +
  labs(
    title = "Supplementary. Off-target enriched DA burden across all community effective fractions",
    x = "Effective per-species fraction",
    y = "Off-target enriched features\nmedian and IQR across runs"
  ) +
  pub_theme(10.7) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 8.7),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#EFEFEF", linewidth = 0.25)
  )

# ----------------------------
# Main figure and outputs
# ----------------------------
main_title <- "Community spike-ins reproduce differential-abundance trends while increasing target detection"
# Standalone overview: landscape layout rather than a long vertical stack.
# The manuscript assembly script uses the individual exported panels, but this
# makes manuscript_community_DA_recovery_overview easier to inspect as well.
main_fig <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1.05, 1.00), heights = c(1.00, 0.82), guides = "collect") +
  plot_annotation(
    title = main_title,
    theme = theme(plot.title = element_text(face = "bold", size = 18.0, hjust = 0, colour = "#111111"))
  ) &
  theme(legend.position = "bottom")

save_plot_set(pA, "panel_A_independent_vs_community_target_DA_detection_dumbbell", width = 12.0, height = 6.5, dpi = opt$dpi)
save_plot_set(pA_concordance_supp, "supp_panel_A_background_resolved_DA_concordance", width = 11.5, height = 6.8, dpi = opt$dpi)
save_plot_set(pB, "panel_B_community_target_DA_detection_main_effective_fractions", width = 12.5, height = 5.5, dpi = opt$dpi)
save_plot_set(pB_all, "supp_panel_B_community_target_DA_detection_all_effective_fractions", width = 15.0, height = 5.5, dpi = opt$dpi)
save_plot_set(pB_q_supp, "supp_panel_B_community_median_neglog10q_heatmap_all_effective_fractions", width = 10.5, height = 7.2, dpi = opt$dpi)
save_plot_set(pC, "panel_C_study_specific_community_target_DA_detection", width = 8.8, height = 4.7, dpi = opt$dpi)
save_plot_set(pC_species_supp, "supp_panel_C_study_resolved_community_target_DA_detection_by_species", width = 12.5, height = 6.8, dpi = opt$dpi)
save_plot_set(pD, "panel_D_community_offtarget_enriched_DA_burden_main_effective_fractions", width = 8.8, height = 4.4, dpi = opt$dpi)
save_plot_set(pD_all, "supp_panel_D_community_offtarget_enriched_DA_burden_all_effective_fractions", width = 9.5, height = 4.4, dpi = opt$dpi)
save_plot_set(main_fig, "manuscript_community_DA_recovery_overview", width = opt$`main-width`, height = opt$`main-height`, dpi = opt$dpi)

readme <- c(
  "Run example:",
  "Rscript scripts/04_plot_manuscript_community_DA_recovery.R \\",
  "  --indir RUNS/maaslin_spike \\",
  "  --outdir RUNS/plots_manuscript_community_DA_recovery \\",
  "  --community-size 10 \\",
  "  --main-effective-fractions 0.00001,0.00005,0.0001,0.0005,0.001",
  "",
  "Main outputs:",
  "- manuscript_community_DA_recovery_overview.(pdf|png)",
  "- panel_A_independent_vs_community_target_DA_detection_dumbbell.(pdf|png)",
  "- panel_B_community_target_DA_detection_main_effective_fractions.(pdf|png)",
  "- panel_C_study_specific_community_target_DA_detection.(pdf|png)",
  "- panel_D_community_offtarget_enriched_DA_burden_main_effective_fractions.(pdf|png)",
  "",
  "Supplementary outputs:",
  "- supp_panel_A_background_resolved_DA_concordance.(pdf|png)",
  "- supp_panel_B_community_target_DA_detection_all_effective_fractions.(pdf|png)",
  "- supp_panel_B_community_median_neglog10q_heatmap_all_effective_fractions.(pdf|png)",
  "- supp_panel_C_study_resolved_community_target_DA_detection_by_species.(pdf|png)",
  "- supp_panel_D_community_offtarget_enriched_DA_burden_all_effective_fractions.(pdf|png)",
  "",
  "Notes:",
  "- Target DA detection is mean(member_detected_positive), i.e. the target was detected as significant with a positive spike-associated coefficient.",
  "- Community total spike fractions are divided by --community-size to obtain effective per-species fractions.",
  "- Main Panel A uses matched effective fractions shared with independent spikes: 0.01%, 0.05%, 0.10% by default.",
  "- Main Panels B/C/D use --main-effective-fractions, defaulting to 0.001%, 0.005%, 0.01%, 0.05%, 0.10%.",
  "- Panel C summarizes study-specific target DA detection as the median and IQR across the ten spiked species."
)
writeLines(readme, con = file.path(opt$outdir, "README_community_DA_recovery.txt"))

message("[OK] Wrote community DA recovery figure set under: ", opt$outdir)
