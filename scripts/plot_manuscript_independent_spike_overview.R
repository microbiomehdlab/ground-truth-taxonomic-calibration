#!/usr/bin/env Rscript
# Publication-polished Fig. 2 source panels: sorted taxa, boxed spike-fraction modules, clearer separators, endpoint-safe axes.

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

option_list <- list(
  make_option("--indir", type = "character", default = "RUNS/spike_metrics",
              help = "Input directory containing target_member_errors_with_condition.csv [default %default]"),
  make_option("--input-file", type = "character", default = NULL,
              help = "Optional direct path to target_member_errors_with_condition.csv"),
  make_option("--outdir", type = "character", default = "RUNS/plots_manuscript_independent_overview_naturestyle",
              help = "Output directory [default %default]"),
  make_option("--weak-fractions", type = "character", default = "0.0001,0.0005,0.001",
              help = "Comma-separated weak spike fractions to show in the main figure [default %default]"),
  make_option("--spike-labels", type = "character", default = NULL,
              help = "Optional comma-separated spike labels to include; default = all independent labels except CRCpanel"),
  make_option("--sort-labels-by", type = "character", default = "taxon",
              help = "Ordering for species labels: 'taxon' uses the alphabetical taxon-name order used in the manuscript table; 'label' uses label alphabetical order; 'input' preserves --spike-labels order; 'contrast' orders by Kraken2-Bracken minus MetaPhlAn good-recovery contrast [default %default]"),
  make_option("--good-threshold", type = "double", default = 0.10,
              help = "Absolute relative-error threshold for Good class [default %default]"),
  make_option("--average-threshold", type = "double", default = 0.50,
              help = "Absolute relative-error threshold for Intermediate class; larger values are Poor/missed [default %default]"),
  make_option("--main-width", type = "double", default = 13.6,
              help = "Main figure width in inches [default %default]"),
  make_option("--main-height", type = "double", default = 10.2,
              help = "Main figure height in inches [default %default]"),
  make_option("--supp-width", type = "double", default = 13,
              help = "Supplementary figure width in inches [default %default]"),
  make_option("--supp-height", type = "double", default = 9.5,
              help = "Supplementary figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 450,
              help = "PNG resolution [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

input_file <- opt$`input-file` %||% file.path(opt$indir, "target_member_errors_with_condition.csv")
if (!file.exists(input_file)) stop("Could not find target-member table at: ", input_file, call. = FALSE)

parse_num_vector <- function(x) {
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)] %>% as.numeric()
}

parse_chr_vector <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  x %>% strsplit(",") %>% unlist() %>% trimws() %>% .[nzchar(.)]
}

fraction_label <- function(x) sprintf("%.2f%%", x * 100)

weak_fractions <- parse_num_vector(opt$`weak-fractions`)
requested_labels <- parse_chr_vector(opt$`spike-labels`)

tool_labeller <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)
tool_levels <- unname(tool_labeller)

fraction_cols <- c(
  "0.01%" = "#4E79A7",
  "0.05%" = "#2A9D8F",
  "0.10%" = "#E69F00",
  "0.50%" = "#7B61D1",
  "1.00%" = "#B65B84",
  "5.00%" = "#666666"
)

pub_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "#6F6F6F", fill = NA, linewidth = 0.46),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.28),
      panel.grid.major.x = element_line(colour = "#F1F1F1", linewidth = 0.24),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.30),
      axis.ticks.length = unit(2.0, "pt"),
      strip.background = element_rect(fill = "#EDEDED", colour = "#6F6F6F", linewidth = 0.50),
      strip.text = element_text(face = "bold", colour = "#222222", margin = margin(4.8, 6, 4.8, 6)),
      axis.title = element_text(face = "bold", colour = "#222222"),
      axis.text = element_text(colour = "#333333"),
      legend.title = element_text(face = "bold", colour = "#222222", size = rel(0.92)),
      legend.text = element_text(colour = "#333333", size = rel(0.88)),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.height = unit(9.5, "pt"),
      legend.key.width = unit(13, "pt"),
      plot.title = element_text(face = "bold", size = base_size + 1.7, hjust = 0, colour = "#111111", margin = margin(b = 5)),
      plot.subtitle = element_blank(),
      plot.margin = margin(7, 8, 7, 8)
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
  saveRDS(plot_obj, file.path(opt$outdir, paste0(stem, ".rds")))
}

message("[INFO] Reading: ", input_file)
raw_tbl <- read_csv(input_file, show_col_types = FALSE)

required_cols <- c("spike_mode", "spike_label", "tool", "spike_fraction_total", "relative_error", "observed_over_expected")
missing_cols <- setdiff(required_cols, names(raw_tbl))
if (length(missing_cols) > 0) {
  stop("Input table is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

tbl <- raw_tbl %>%
  filter(spike_mode == "independent") %>%
  filter(!is.na(spike_label), spike_label != "CRCpanel") %>%
  mutate(
    tool = recode(tool, !!!tool_labeller),
    spike_fraction_total = as.numeric(spike_fraction_total),
    abs_relative_error = abs(as.numeric(relative_error)),
    observed_over_expected = as.numeric(observed_over_expected)
  )

if ("target_flag" %in% names(tbl)) {
  tbl <- tbl %>% filter(is.na(target_flag) | target_flag)
}

if (!is.null(requested_labels)) {
  tbl <- tbl %>% filter(spike_label %in% requested_labels)
}

available_labels <- sort(unique(tbl$spike_label))
if (length(available_labels) == 0) stop("No independent-spike rows remained after filtering.", call. = FALSE)

# Manuscript order used for the assembly table, sorted by full taxon name:
# Bacteroides, Clostridium, Dialister, Fusobacterium, Hungatella,
# Parvimonas, Peptostreptococcus anaerobius, Peptostreptococcus stomatis,
# Porphyromonas, Prevotella.
taxon_alphabetical_levels <- c("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic", "Pana", "Psto", "Porp", "Pint")
label_alphabetical_levels <- sort(available_labels)
input_levels <- requested_labels %||% available_labels
input_levels <- input_levels[input_levels %in% available_labels]

sort_mode <- tolower(opt$`sort-labels-by` %||% "taxon")
if (sort_mode == "label") {
  label_levels <- label_alphabetical_levels[label_alphabetical_levels %in% available_labels]
} else if (sort_mode == "input") {
  label_levels <- input_levels
} else if (sort_mode == "taxon") {
  label_levels <- c(taxon_alphabetical_levels[taxon_alphabetical_levels %in% available_labels],
                    setdiff(label_alphabetical_levels, taxon_alphabetical_levels))
} else if (sort_mode == "contrast") {
  label_levels <- input_levels
} else {
  stop("--sort-labels-by must be one of: taxon, label, input, contrast", call. = FALSE)
}

message("[INFO] Species label order for plotting: ", paste(label_levels, collapse = ", "))

tbl <- tbl %>%
  mutate(
    spike_label = factor(spike_label, levels = label_levels),
    tool = factor(tool, levels = tool_levels),
    fraction_label = factor(fraction_label(spike_fraction_total), levels = fraction_label(sort(unique(spike_fraction_total))))
  )

classify_recovery <- function(abs_rel_err, ratio, good_thr, avg_thr) {
  case_when(
    is.na(ratio) | is.na(abs_rel_err) ~ "Poor / missed",
    abs_rel_err <= good_thr ~ "Good",
    abs_rel_err <= avg_thr ~ "Intermediate",
    TRUE ~ "Poor / missed"
  )
}

tbl <- tbl %>%
  mutate(
    recovery_class = classify_recovery(abs_relative_error, observed_over_expected,
                                       opt$`good-threshold`, opt$`average-threshold`),
    recovery_class = factor(recovery_class, levels = c("Good", "Intermediate", "Poor / missed"))
  )

main_tbl <- tbl %>% filter(spike_fraction_total %in% weak_fractions)
if (nrow(main_tbl) == 0) {
  stop("No rows matched the requested weak fractions: ", paste(weak_fractions, collapse = ", "), call. = FALSE)
}
main_tbl <- main_tbl %>%
  mutate(fraction_label = factor(fraction_label(spike_fraction_total), levels = fraction_label(weak_fractions)))

good_rate_summary <- function(dat) {
  dat %>%
    group_by(tool, spike_label, spike_fraction_total, fraction_label) %>%
    summarise(
      n = n(),
      good_rate = mean(recovery_class == "Good", na.rm = TRUE),
      avg_rate = mean(recovery_class == "Intermediate", na.rm = TRUE),
      bad_rate = mean(recovery_class == "Poor / missed", na.rm = TRUE),
      .groups = "drop"
    )
}

bias_spread_summary <- function(dat) {
  dat %>%
    group_by(tool, spike_label, spike_fraction_total, fraction_label) %>%
    summarise(
      n = sum(!is.na(observed_over_expected)),
      median_recovery = median(observed_over_expected, na.rm = TRUE),
      iqr_recovery = IQR(observed_over_expected, na.rm = TRUE),
      good_rate = mean(recovery_class == "Good", na.rm = TRUE),
      .groups = "drop"
    )
}

median_recovery_summary <- function(dat) {
  dat %>%
    group_by(tool, spike_label, spike_fraction_total, fraction_label) %>%
    summarise(
      n = sum(!is.na(observed_over_expected)),
      median_recovery = median(observed_over_expected, na.rm = TRUE),
      q1 = quantile(observed_over_expected, 0.25, na.rm = TRUE),
      q3 = quantile(observed_over_expected, 0.75, na.rm = TRUE),
      iqr_recovery = IQR(observed_over_expected, na.rm = TRUE),
      .groups = "drop"
    )
}

class_comp_summary <- function(dat) {
  dat %>%
    group_by(tool, spike_label, spike_fraction_total, fraction_label, recovery_class) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

good_main <- good_rate_summary(main_tbl)
bias_main <- bias_spread_summary(main_tbl)
median_main <- median_recovery_summary(main_tbl)
class_main <- class_comp_summary(main_tbl)

good_all <- good_rate_summary(tbl)
bias_all <- bias_spread_summary(tbl)
median_all <- median_recovery_summary(tbl)
class_all <- class_comp_summary(tbl)

write_csv(good_main, file.path(opt$outdir, "overview_main_good_rate_summary.csv"))
write_csv(bias_main, file.path(opt$outdir, "overview_main_bias_spread_summary.csv"))
write_csv(median_main, file.path(opt$outdir, "overview_main_median_recovery_summary.csv"))
write_csv(class_main, file.path(opt$outdir, "overview_main_class_composition_summary.csv"))
write_csv(good_all, file.path(opt$outdir, "overview_allfractions_good_rate_summary.csv"))
write_csv(bias_all, file.path(opt$outdir, "overview_allfractions_bias_spread_summary.csv"))
write_csv(median_all, file.path(opt$outdir, "overview_allfractions_median_recovery_summary.csv"))
write_csv(class_all, file.path(opt$outdir, "overview_allfractions_class_composition_summary.csv"))

# Main manuscript order: alphabetical by full taxon name by default.
# The old contrast-based order can still be requested with --sort-labels-by contrast.
if (sort_mode == "contrast") {
  species_order <- good_main %>%
    group_by(spike_label, tool) %>%
    summarise(mean_good = mean(good_rate, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = tool, values_from = mean_good, values_fill = 0) %>%
    mutate(
      kraken_good = coalesce(.data[["Kraken2 + Bracken"]], 0),
      metaphlan_good = coalesce(.data[["MetaPhlAn 4"]], 0),
      profiler_gap = kraken_good - metaphlan_good,
      mean_good_all = (kraken_good + metaphlan_good) / 2
    ) %>%
    arrange(desc(profiler_gap), desc(kraken_good), desc(mean_good_all), spike_label) %>%
    pull(spike_label) %>%
    as.character()
} else {
  species_order <- label_levels
}
if (length(species_order) == 0) species_order <- label_levels

# ggplot draws lower factor levels at the bottom; reversing gives alphabetical top-to-bottom order.
display_levels <- rev(species_order)

good_main <- good_main %>%
  mutate(
    spike_label_plot = factor(as.character(spike_label), levels = display_levels),
    fraction_label = factor(fraction_label, levels = fraction_label(weak_fractions))
  )

# ----------------------------
# Panel A: paired dumbbell plot of recovery-class rates
# ----------------------------
profiler_cols <- c(
  "Kraken2 + Bracken" = "#009E73",
  "MetaPhlAn 4" = "#6F5BD3"
)

make_panelA_dumbbell <- function(rate_tbl, rate_col, title_text, x_label) {
  pair_tbl <- rate_tbl %>%
    mutate(
      tool_short = case_when(
        tool == "Kraken2 + Bracken" ~ "kraken",
        tool == "MetaPhlAn 4" ~ "metaphlan",
        TRUE ~ as.character(tool)
      )
    ) %>%
    group_by(spike_label_plot, fraction_label, tool_short) %>%
    summarise(rate_value = mean(.data[[rate_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = tool_short,
      values_from = rate_value,
      values_fill = 0
    ) %>%
    mutate(
      spike_label_plot = factor(as.character(spike_label_plot), levels = display_levels)
    )
  
  ggplot(pair_tbl, aes(y = spike_label_plot)) +
    geom_segment(
      aes(x = kraken, xend = metaphlan, yend = spike_label_plot),
      colour = "#B9B9B9",
      linewidth = 0.62,
      lineend = "round"
    ) +
    geom_point(
      aes(x = kraken, colour = "Kraken2 + Bracken", shape = "Kraken2 + Bracken"),
      size = 2.9,
      alpha = 0.98,
      position = position_nudge(y = 0.10)
    ) +
    geom_point(
      aes(x = metaphlan, colour = "MetaPhlAn 4", shape = "MetaPhlAn 4"),
      size = 2.9,
      alpha = 0.98,
      position = position_nudge(y = -0.10)
    ) +
    geom_vline(
      xintercept = 0.5,
      linetype = 2,
      colour = "#6B6B6B",
      linewidth = 0.50
    ) +
    facet_grid(. ~ fraction_label) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = c(-0.045, 1.045), clip = "off") +
    scale_colour_manual(values = profiler_cols, name = "Profiler") +
    scale_shape_manual(values = c("Kraken2 + Bracken" = 16, "MetaPhlAn 4" = 17), name = "Profiler") +
    labs(
      title = title_text,
      x = x_label,
      y = NULL
    ) +
    pub_theme(11.0) +
    theme(
      axis.text.y = element_text(size = 9.5),
      axis.text.x = element_text(size = 9.0),
      panel.grid.major.y = element_line(colour = "#F1F1F1", linewidth = 0.26),
      panel.grid.major.x = element_line(colour = "#F4F4F4", linewidth = 0.24),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1.45, "lines"),
      strip.background = element_rect(fill = "#EDEDED", colour = "#5F5F5F", linewidth = 0.55),
      strip.text = element_text(margin = margin(5, 7, 5, 7), face = "bold"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(5, 8, 5, 5)
    )
}

pA <- make_panelA_dumbbell(
  good_main,
  rate_col = "good_rate",
  title_text = "A. Reliable weak-signal recovery across species",
  x_label = "% samples with good recovery (|relative error| ≤ 10%)"
)

pA_avg <- make_panelA_dumbbell(
  good_main,
  rate_col = "avg_rate",
  title_text = "Supplementary. Intermediate recovery across species",
  x_label = "% samples with intermediate recovery (10% < |relative error| ≤ 50%)"
)

pA_bad <- make_panelA_dumbbell(
  good_main,
  rate_col = "bad_rate",
  title_text = "Supplementary. Poor or missed recovery across species",
  x_label = "% samples with poor or missed recovery (|relative error| > 50% or undetected)"
)

# ----------------------------
# Panel B: aligned species-level bias and variability dot plots
# ----------------------------
fraction_levels_main <- fraction_label(weak_fractions)
fraction_cols_main <- fraction_cols[fraction_levels_main]
missing_fraction_cols <- setdiff(fraction_levels_main, names(fraction_cols_main))
if (length(missing_fraction_cols) > 0) {
  fraction_cols_main <- c(fraction_cols_main, setNames(hue_pal()(length(missing_fraction_cols)), missing_fraction_cols))
}

offset_map <- setNames(seq(-0.22, 0.22, length.out = length(fraction_levels_main)), fraction_levels_main)

bias_long <- bias_main %>%
  mutate(
    fraction_label = factor(fraction_label, levels = fraction_levels_main),
    spike_label_chr = as.character(spike_label),
    y_base = match(spike_label_chr, display_levels),
    y_plot = y_base + unname(offset_map[as.character(fraction_label)])
  ) %>%
  select(tool, spike_label, spike_label_chr, fraction_label, median_recovery, iqr_recovery, y_base, y_plot) %>%
  pivot_longer(
    cols = c(median_recovery, iqr_recovery),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    median_recovery = "Bias\nmedian observed/expected",
                    iqr_recovery = "Variability\nIQR observed/expected"),
    metric = factor(metric, levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"))
  )

ref_lines <- tibble::tibble(
  metric = factor(c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"),
                  levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected")),
  xintercept = c(1, 0)
)


ideal_bands <- tibble::tibble(
  metric = factor(c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"),
                  levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected")),
  xmin = c(0.9, 0.0),
  xmax = c(1.1, 0.12)
)

pB <- ggplot(bias_long, aes(x = value, y = y_plot, color = fraction_label)) +
  geom_rect(
    data = ideal_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#E6F2EC",
    alpha = 0.85,
    colour = NA
  ) +
  geom_vline(data = ref_lines, aes(xintercept = xintercept), inherit.aes = FALSE,
             linetype = 2, colour = "#555555", linewidth = 0.50) +
  geom_hline(yintercept = seq_along(display_levels), colour = "#E8E8E8", linewidth = 0.28) +
  geom_line(aes(group = interaction(tool, metric, spike_label_chr)),
            colour = "#BFBFBF", linewidth = 0.38, alpha = 0.55, na.rm = TRUE) +
  geom_point(size = 2.85, alpha = 0.96, stroke = 0.25, na.rm = TRUE) +
  facet_grid(tool ~ metric, scales = "free_x") +
  scale_y_continuous(
    breaks = seq_along(display_levels),
    labels = display_levels,
    expand = expansion(mult = c(0.035, 0.035))
  ) +
  scale_color_manual(values = fraction_cols_main, name = "Spike fraction") +
  labs(
    title = "B. Species-level recovery bias and variability",
    x = NULL,
    y = NULL
  ) +
  pub_theme(10.8) +
  theme(
    axis.text.y = element_text(size = 9.3),
    axis.text.x = element_text(size = 8.9),
    strip.text = element_text(size = 9.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.24),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.spacing.x = unit(1.25, "lines"),
    panel.spacing.y = unit(0.85, "lines"),
    strip.background = element_rect(fill = "#EDEDED", colour = "#5F5F5F", linewidth = 0.52),
    panel.border = element_rect(colour = "#6F6F6F", fill = NA, linewidth = 0.46),
    plot.margin = margin(5, 8, 5, 5)
  )

main_title <- "Independent spike-ins reveal strong profiler differences in weak-signal recovery"
main_fig <- (pA / pB) +
  plot_layout(heights = c(0.85, 1.0), guides = "collect") +
  plot_annotation(
    title = main_title,
    theme = theme(plot.title = element_text(face = "bold", size = 18.0, hjust = 0, colour = "#111111"))
  ) &
  theme(legend.position = "right")

# ----------------------------
# Supplementary outputs: retain detailed all-fraction summaries
# ----------------------------
fraction_levels_all <- fraction_label(sort(unique(tbl$spike_fraction_total)))
good_all <- good_all %>%
  mutate(
    fraction_label = factor(fraction_label, levels = fraction_levels_all),
    spike_label_plot = factor(as.character(spike_label), levels = display_levels)
  )
bias_all <- bias_all %>% mutate(fraction_label = factor(fraction_label, levels = fraction_levels_all))
median_all <- median_all %>%
  mutate(
    fraction_label = factor(fraction_label, levels = fraction_levels_all),
    spike_label_plot = factor(as.character(spike_label), levels = display_levels)
  )
class_all <- class_all %>% mutate(fraction_label = factor(fraction_label, levels = fraction_levels_all))

supp_good <- ggplot(good_all, aes(x = fraction_label, y = spike_label_plot)) +
  geom_point(aes(fill = good_rate, size = good_rate), shape = 21, colour = "#3F3F3F", stroke = 0.20, alpha = 0.95) +
  facet_grid(tool ~ ., switch = "y") +
  scale_fill_gradient(low = "#EEF6F2", high = "#009E73", limits = c(0, 1), labels = percent_format(accuracy = 1), name = "Good recovery") +
  scale_size_continuous(range = c(1.8, 7.0), limits = c(0, 1), guide = "none") +
  labs(title = "Good-recovery rate across all independent taxa and spike fractions",
       subtitle = "Good recovery: absolute relative recovery error ≤10%",
       x = "Spike fraction", y = NULL) +
  pub_theme(10.8) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
        axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")

supp_bias_long <- bias_all %>%
  mutate(
    spike_label_chr = as.character(spike_label),
    y_base = match(spike_label_chr, display_levels),
    fraction_index = as.integer(fraction_label),
    fraction_count = nlevels(fraction_label),
    # Separate fraction-specific estimates slightly on the categorical y-axis
    # so coincident values do not conceal one another.
    y_plot = y_base + (fraction_index - (fraction_count + 1) / 2) * 0.055
  ) %>%
  pivot_longer(cols = c(median_recovery, iqr_recovery), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    median_recovery = "Bias\nmedian observed/expected",
                    iqr_recovery = "Variability\nIQR observed/expected"),
    metric = factor(metric, levels = c("Bias\nmedian observed/expected", "Variability\nIQR observed/expected"))
  )

supp_bias <- ggplot(supp_bias_long, aes(x = value, y = y_plot, color = fraction_label)) +
  geom_vline(data = ref_lines, aes(xintercept = xintercept), inherit.aes = FALSE,
             linetype = 2, colour = "#555555", linewidth = 0.50) +
  geom_hline(yintercept = seq_along(display_levels), colour = "#E8E8E8", linewidth = 0.28) +
  geom_point(size = 2.1, alpha = 0.88, na.rm = TRUE) +
  facet_grid(tool ~ metric, scales = "free_x") +
  scale_y_continuous(breaks = seq_along(display_levels), labels = display_levels,
                     expand = expansion(mult = c(0.035, 0.035))) +
  labs(title = "Recovery bias and variability across all independent taxa and spike fractions",
       subtitle = "Dashed reference lines mark no bias (observed/expected = 1) and zero variability",
       x = NULL, y = NULL, color = "Spike fraction") +
  pub_theme(10.8) +
  theme(axis.text.x = element_text(size = 8.8), legend.position = "right")

supp_median <- ggplot(median_all, aes(x = fraction_label, y = spike_label_plot, fill = median_recovery)) +
  geom_tile(color = "white", linewidth = 0.55) +
  facet_wrap(~tool, ncol = 1) +
  scale_fill_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027", midpoint = 1,
                       name = "Median\nobs/exp", breaks = c(0.5, 1, 1.5, 2)) +
  labs(title = "Median observed-to-expected recovery across all independent taxa and spike fractions",
       x = "Spike fraction", y = NULL) +
  pub_theme(10.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

class_cols <- c("Good" = "#009E73", "Intermediate" = "#E69F00", "Poor / missed" = "#D55E00")
class_all_display <- class_all %>%
  mutate(
    recovery_class_display = dplyr::recode(
      as.character(recovery_class),
      "Average" = "Intermediate",
      "Bad / missed" = "Poor / missed"
    ),
    recovery_class_display = factor(
      recovery_class_display,
      levels = c("Good", "Intermediate", "Poor / missed")
    )
  )
class_taxon_blocks <- split(label_levels, ceiling(seq_along(label_levels) / 5))
make_class_block <- function(block) {
  dat <- class_all_display %>%
    filter(as.character(spike_label) %in% block) %>%
    mutate(spike_label = factor(as.character(spike_label), levels = block))

  ggplot(dat, aes(x = fraction_label, y = prop, fill = recovery_class_display)) +
    geom_col(width = 0.85, color = "white", linewidth = 0.25) +
    facet_grid(tool ~ spike_label) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = class_cols, name = "Recovery class", drop = FALSE) +
    labs(x = "Spike fraction", y = "% of samples") +
    pub_theme(11.2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}
supp_class <- wrap_plots(
  lapply(class_taxon_blocks, make_class_block),
  ncol = 1,
  guides = "collect"
) +
  plot_annotation(
    title = "Recovery-class composition across all independent taxa and spike fractions",
    subtitle = "Good ≤10%; Intermediate >10% and ≤50%; Poor / missed >50% or undetected"
  ) &
  theme(legend.position = "bottom")

save_plot_set(pA, "panel_A_good_recovery_dumbbell_naturestyle", width = 12.0, height = 4.85, dpi = opt$dpi)
save_plot_set(pA_avg, "supp_panel_A_average_recovery_dumbbell_naturestyle", width = 12.0, height = 4.85, dpi = opt$dpi)
save_plot_set(pA_bad, "supp_panel_A_bad_recovery_dumbbell_naturestyle", width = 12.0, height = 4.85, dpi = opt$dpi)
save_plot_set(pB, "panel_B_bias_variability_dotplots_naturestyle", width = 8.1, height = 6.25, dpi = opt$dpi)
save_plot_set(main_fig, "manuscript_independent_spike_overview_naturestyle_AB", width = opt$`main-width`, height = opt$`main-height`, dpi = opt$dpi)

save_plot_set(supp_good, "supp_good_recovery_dot_heatmap_allfractions", width = opt$`supp-width`, height = opt$`supp-height`, dpi = opt$dpi)
save_plot_set(supp_bias, "supp_bias_variability_dotplots_allfractions", width = opt$`supp-width`, height = opt$`supp-height`, dpi = opt$dpi)
save_plot_set(supp_median, "supp_median_recovery_heatmap_allfractions", width = opt$`supp-width`, height = opt$`supp-height`, dpi = opt$dpi)
save_plot_set(supp_class, "supp_recovery_class_composition_allfractions", width = 13, height = 11, dpi = opt$dpi)

readme_lines <- c(
  "Run example:",
  "Rscript scripts/02_plot_manuscript_independent_spike_overview_naturestyle.R \\",
  "  --indir RUNS/spike_metrics \\",
  "  --outdir RUNS/plots_manuscript_independent_overview_naturestyle",
  "",
  "Main outputs:",
  "- manuscript_independent_spike_overview_naturestyle_AB.(pdf|png)",
  "- panel_A_good_recovery_dumbbell_naturestyle.(pdf|png)",
  "- panel_B_bias_variability_dotplots_naturestyle.(pdf|png)",
  "",
  "Supplementary outputs:",
  "- supp_good_recovery_dot_heatmap_allfractions.(pdf|png)",
  "- supp_bias_variability_dotplots_allfractions.(pdf|png)",
  "- supp_median_recovery_heatmap_allfractions.(pdf|png)",
  "- supp_recovery_class_composition_allfractions.(pdf|png)",
  "- CSV summaries for main and all-fraction panels",
  "",
  "Design notes:",
  "- Panel A uses a paired dumbbell plot to compare profilers directly for each species and weak spike fraction.",
  "- Panel B separates bias (median observed/expected) from variability (IQR observed/expected).",
  paste0("- Species order used in this run: ", paste(species_order, collapse = ", "), "."),
  "- Default species order is alphabetical by full taxon name; use --sort-labels-by contrast to recover the previous contrast-sorted order.",
  "- Facet strips and panel borders are intentionally stronger to make each spike-fraction block visually explicit."
)
writeLines(readme_lines, con = file.path(opt$outdir, "README_independent_spike_overview_naturestyle.txt"))

message("[OK] Wrote Nature-style independent-spike overview under: ", opt$outdir)
message("[OK] Main figure weak fractions: ", paste(fraction_label(weak_fractions), collapse = ", "))
message("[OK] Included labels (", length(label_levels), "): ", paste(label_levels, collapse = ", "))
