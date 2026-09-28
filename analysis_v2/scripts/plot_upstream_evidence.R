#!/usr/bin/env Rscript
# Supplementary figure: "Upstream cohort and profile completeness".
#
# Every plotted value is read from the exported figure_source_data.tsv that the
# builder writes from the three verified production_seal_v2 seals. Nothing is
# recomputed here, so the figure and the released table cannot disagree.
#
# Completeness is a property of the analytical inputs. It does not validate
# biological accuracy, and the caption in the package README says so.

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
value_of <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}

source_data <- value_of("--source-data")
outdir <- value_of("--outdir")
if (is.null(source_data) || is.null(outdir))
  stop("Required: --source-data FILE --outdir DIR")
if (!file.exists(source_data) || file.info(source_data)$size == 0)
  stop("Missing or empty source data: ", source_data)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rows <- utils::read.delim(source_data, sep = "\t", stringsAsFactors = FALSE,
                          check.names = FALSE)
required <- c("panel", "cohort", "group", "measure", "value")
missing <- setdiff(required, names(rows))
if (length(missing)) stop("Source data lacks: ", paste(missing, collapse = ", "))
if (!nrow(rows)) stop("Source data has no rows")
rows$value <- as.numeric(rows$value)
if (any(is.na(rows$value))) stop("Source data has a nonnumeric value")
for (panel in c("A", "B", "C", "D")) {
  if (!any(rows$panel == panel)) stop("Source data lacks panel ", panel)
}

# Okabe-Ito: colour-blind safe and separable in grayscale.
OKABE <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")
CONDITIONS <- c("Control", "Adenoma", "CRC")
COHORTS <- c("yachida", "feng", "zeller")

ordered_cohort <- function(x) factor(x, levels = COHORTS)
# Counts are integers, so the axis must never show fractional breaks.
integer_breaks <- function(limits) {
  top <- max(1, ceiling(max(limits, na.rm = TRUE)))
  unique(as.integer(round(pretty(c(0, top)))))
}
theme_upstream <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom", legend.key.size = unit(3.2, "mm"),
        plot.title = element_text(face = "bold", size = 9))

counted_bars <- function(data, fill_levels, title, subtitle, fill_name) {
  data$cohort <- ordered_cohort(data$cohort)
  data$group <- factor(data$group, levels = fill_levels)
  ggplot(data, aes(x = cohort, y = value, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             colour = "grey20", linewidth = 0.2) +
    geom_text(aes(label = format(value, big.mark = "", trim = TRUE)),
              position = position_dodge(width = 0.8), vjust = -0.35,
              size = 2.4) +
    scale_fill_manual(values = stats::setNames(
      OKABE[seq_along(fill_levels)], fill_levels), name = fill_name) +
    scale_y_continuous(breaks = integer_breaks,
                       expand = expansion(mult = c(0, 0.16))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "count") +
    theme_upstream
}

panel_a <- counted_bars(
  rows[rows$panel == "A" & rows$measure == "samples", ],
  CONDITIONS, "A  Samples by cohort and condition",
  "exact sample counts", "condition")

panel_b_data <- rows[rows$panel == "B", ]
panel_b_data$cohort <- ordered_cohort(panel_b_data$cohort)
panel_b_data$group <- factor(panel_b_data$group,
                             levels = c("baseline", "independent", "community"))
panel_b_data$measure <- factor(panel_b_data$measure,
                               levels = c("expected_profiles", "observed_profiles"),
                               labels = c("expected", "observed"))
panel_b <- ggplot(panel_b_data, aes(x = group, y = value, fill = measure)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = value), position = position_dodge(width = 0.8),
            vjust = -0.35, size = 2.1) +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = c(expected = OKABE[1], observed = OKABE[3]),
                    name = NULL) +
  scale_y_continuous(breaks = integer_breaks,
                     expand = expansion(mult = c(0, 0.18))) +
  labs(title = "B  Expected and observed profiles by cohort and design",
       subtitle = "exact profile counts", x = NULL, y = "profiles") +
  theme_upstream +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

panel_c_data <- rows[rows$panel == "C" & rows$measure == "observed_samples", ]
panel_c <- counted_bars(panel_c_data, CONDITIONS,
                        "C  Independent-subset samples by cohort and condition",
                        "dashed line marks the frozen expectation of 10",
                        "condition") +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.35,
             colour = "grey25")

panel_d_data <- rows[rows$panel == "D" & rows$measure == "missing", ]
panel_d_data$cohort <- ordered_cohort(panel_d_data$cohort)
fields <- sort(unique(panel_d_data$group))
panel_d <- ggplot(panel_d_data, aes(x = cohort, y = value, fill = group)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78,
           colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = value), position = position_dodge(width = 0.85),
            vjust = -0.35, size = 2.1) +
  scale_fill_manual(values = stats::setNames(
    rep(OKABE, length.out = length(fields)), fields), name = "field") +
  scale_y_continuous(breaks = integer_breaks,
                     expand = expansion(mult = c(0, 0.18))) +
  labs(title = "D  Aggregate covariate missingness by cohort and field",
       subtitle = "aggregate counts only; no individual covariate values",
       x = NULL, y = "missing values") +
  theme_upstream

panels <- list(panel_a, panel_b, panel_c, panel_d)
HEIGHTS <- c(1, 1.15, 1, 1)
WIDTH <- 7.2
HEIGHT <- 10.4

# Compose with patchwork or gridExtra when present, otherwise with base grid,
# which ships with R. The figure must not depend on an optional package.
draw_figure <- function() {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(panels, ncol = 1, heights = HEIGHTS))
    return(invisible(NULL))
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(grobs = panels, ncol = 1, heights = HEIGHTS)
    return(invisible(NULL))
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(length(panels), 1,
                               heights = grid::unit(HEIGHTS, "null"))))
  for (index in seq_along(panels)) {
    print(panels[[index]], vp = grid::viewport(layout.pos.row = index,
                                               layout.pos.col = 1))
  }
  grid::popViewport()
  invisible(NULL)
}

open_device <- function(name) {
  path <- file.path(outdir, name)
  if (grepl("[.]pdf$", name)) {
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(path, width = WIDTH, height = HEIGHT)
    } else {
      grDevices::pdf(path, width = WIDTH, height = HEIGHT)
    }
  } else if (grepl("[.]svg$", name)) {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(path, width = WIDTH, height = HEIGHT)
    } else if (capabilities("cairo")) {
      grDevices::svg(path, width = WIDTH, height = HEIGHT)
    } else {
      stop("No SVG device is available; install svglite or a cairo-enabled R.")
    }
  } else {
    grDevices::png(path, width = WIDTH, height = HEIGHT, units = "in",
                   res = 300, bg = "white",
                   type = if (capabilities("cairo")) "cairo" else "Xlib")
  }
  path
}

for (name in c("upstream_completeness.pdf", "upstream_completeness.svg",
               "upstream_completeness.png")) {
  path <- open_device(name)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  draw_figure()
  grDevices::dev.off()
  if (!file.exists(path) || file.info(path)$size == 0)
    stop("Device produced no output: ", name)
}
cat(sprintf("[PASS] Upstream completeness figure: %s\n", outdir))
