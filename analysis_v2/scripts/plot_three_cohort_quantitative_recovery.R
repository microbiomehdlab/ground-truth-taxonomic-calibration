#!/usr/bin/env Rscript
# Development figure: community-spike quantitative recovery on the selected
# profiler scale. Zero observations have undefined log2 recovery and are
# reported in the source summary, never silently assigned a pseudocount.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}
input <- value("--endpoints")
outdir <- value("--outdir")
reference_scale <- if ("--reference-scale" %in% args) value("--reference-scale") else "profiler_scale"
if (!reference_scale %in% c("profiler_scale", "read_proportional"))
  stop("--reference-scale must be profiler_scale or read_proportional")
if (dir.exists(outdir)) stop("Output directory already exists: ", outdir)
if (!file.exists(input)) stop("Missing endpoints: ", input)
x <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("cohort", "condition", "analysis_population", "assembly_arm",
              "profiler", "target_label", "spike_fraction_target", "sample_id",
              "reference_type", "observed_abundance_fraction",
              "expected_abundance_profiler_scale")
if (reference_scale == "read_proportional")
  required <- c(required, "read_proportional_reference",
                "baseline_abundance_fraction", "spike_fraction_total")
missing <- setdiff(required, names(x))
if (length(missing)) stop("Missing endpoint columns: ", paste(missing, collapse = ", "))
x <- x[x$analysis_population == "community" & x$assembly_arm == "original", ]
if (!nrow(x)) stop("No original-assembly community endpoints")
x$spike_fraction_target <- as.numeric(x$spike_fraction_target)
x$observed_abundance_fraction <- as.numeric(x$observed_abundance_fraction)
x$expected_abundance_profiler_scale <- as.numeric(x$expected_abundance_profiler_scale)
if (reference_scale == "read_proportional") {
  for (field in c("read_proportional_reference", "baseline_abundance_fraction",
                  "spike_fraction_total")) x[[field]] <- as.numeric(x[[field]])
  if (any(!is.finite(x$read_proportional_reference)) ||
      any(!is.finite(x$baseline_abundance_fraction)) ||
      any(!is.finite(x$spike_fraction_total)) ||
      any(abs(x$read_proportional_reference -
              ((1 - x$spike_fraction_total) * x$baseline_abundance_fraction +
                 x$spike_fraction_target)) > 1e-8))
    stop("Read-proportional reference contradicts the original equation")
}
if (any(!is.finite(x$spike_fraction_target)) ||
    any(!is.finite(x$observed_abundance_fraction)) ||
    any(!is.finite(x$expected_abundance_profiler_scale)) ||
    any(x$observed_abundance_fraction < 0) ||
    any(x$expected_abundance_profiler_scale <= 0)) stop("Invalid endpoint numbers")
cohorts <- c("feng", "yachida", "zeller")
profilers <- c("kraken2_bracken", "metaphlan4")
conditions <- c("Control", "Adenoma", "CRC")
doses <- c(.0001, .0005, .001)
labels <- c("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic",
            "Pana", "Psto", "Porp", "Pint")
if (any(!x$cohort %in% cohorts) || any(!x$profiler %in% profilers) ||
    any(!x$condition %in% conditions) || any(!x$target_label %in% labels))
  stop("Unexpected cohort, profiler, condition, or target label")
x <- x[vapply(x$spike_fraction_target, function(v)
  any(abs(v - doses) <= doses * .05), logical(1)), ]
if (!nrow(x)) stop("No rows at selected doses")
x$dose <- doses[vapply(x$spike_fraction_target, function(v)
  which.min(abs(v - doses)), integer(1))]
expected_type <- ifelse(x$profiler == "metaphlan4", "genome_equivalent", "read_proportional")
if (any(x$reference_type != expected_type)) stop("Wrong endpoint reference type")
key <- paste(x$cohort, x$condition, x$profiler, x$target_label,
             x$dose, x$sample_id, sep = "\037")
if (anyDuplicated(key)) stop("Duplicate cohort/condition/profiler/target/dose/sample")
contexts <- expand.grid(cohort = cohorts, condition = conditions,
                        profiler = profilers, target_label = labels, dose = doses,
                        stringsAsFactors = FALSE)
group <- function(d) paste(d$cohort, d$condition, d$profiler,
                           d$target_label, d$dose, sep = "\037")
if (!setequal(group(x), group(contexts))) stop("Incomplete figure context grid")
x$recovery <- ifelse(x$observed_abundance_fraction > 0,
                     log2(x$observed_abundance_fraction /
                          if (reference_scale == "read_proportional")
                            x$read_proportional_reference else
                            x$expected_abundance_profiler_scale), NA_real_)
summaries <- lapply(split(x, group(x)), function(part) {
  positive <- part$recovery[is.finite(part$recovery)]
  data.frame(cohort = part$cohort[1], condition = part$condition[1],
             profiler = part$profiler[1], target_label = part$target_label[1],
             dose = part$dose[1], n_total = nrow(part),
             n_positive = length(positive), n_zero = sum(is.na(part$recovery)),
             q1 = if (length(positive)) unname(quantile(positive, .25)) else NA_real_,
             median = if (length(positive)) median(positive) else NA_real_,
             q3 = if (length(positive)) unname(quantile(positive, .75)) else NA_real_)
})
s <- do.call(rbind, summaries)
if (!nrow(s) || all(s$n_positive == 0)) stop("No positive observations for recovery plot")
dir.create(outdir, recursive = TRUE)
write.table(s, file.path(outdir, "quantitative_recovery_box_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c(paste0("reference_scale\t", reference_scale),
             "read_proportional_equation\tE = (1 - F) b + f",
             "F\ttotal implanted fraction", "b\tbaseline abundance fraction",
             "f\timplanted target fraction"),
           file.path(outdir, "reference_definition.tsv"))
hashes <- system2("sha256sum",
  c(shQuote(normalizePath(input)),
    shQuote(normalizePath(file.path(outdir, "quantitative_recovery_box_summary.tsv"))),
    shQuote(normalizePath(file.path(outdir, "reference_definition.tsv")))),
  stdout = TRUE, stderr = TRUE)
if (!is.null(attr(hashes, "status")) && attr(hashes, "status") != 0)
  stop("Could not hash recovery figure inputs")
writeLines(hashes, file.path(outdir, "source.sha256"))
s$cohort <- factor(s$cohort, levels = cohorts, labels = c("Feng", "Yachida", "Zeller"))
s$condition <- factor(s$condition, levels = conditions)
s$profiler <- factor(s$profiler, levels = profilers,
                     labels = c("Kraken2 + Bracken", "MetaPhlAn 4"))
s$target_label <- factor(s$target_label, levels = labels)
s$dose_label <- factor(sprintf("%.3f%%", 100 * s$dose),
                       levels = sprintf("%.3f%%", 100 * doses))
s$x <- as.numeric(s$target_label) + (as.numeric(s$cohort) - 2) * .24 +
  (as.numeric(s$condition) - 2) * .055
colors <- c(Feng = "#E66B50", Yachida = "#209985", Zeller = "#417A9C")
shapes <- c(Control = 16, Adenoma = 17, CRC = 15)
p <- ggplot(s, aes(x = x, color = cohort)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_rect(aes(xmin = x - .0275, xmax = x + .0275, ymin = q1, ymax = q3,
                fill = cohort), alpha = .18, color = NA, na.rm = TRUE) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = .055, linewidth = .35,
                na.rm = TRUE) +
  geom_point(aes(y = median, shape = condition), size = 1.6, na.rm = TRUE) +
  facet_grid(profiler ~ dose_label) +
  scale_x_continuous(breaks = seq_along(labels), labels = labels,
                     limits = c(.5, length(labels) + .5)) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  labs(title = "Quantitative recovery is profiler- and taxon-dependent",
       subtitle = if (reference_scale == "read_proportional")
         "Original read-proportional expectation for both profilers; medians and sample Q1-Q3 among positive reports" else
         "Community spike-ins: points are medians; shaded boxes span Q1-Q3 across positive samples",
       x = "Implanted taxon", y = expression(log[2]("observed / expected")),
       color = "Cohort", fill = "Cohort", shape = "Condition",
       caption = if (reference_scale == "read_proportional")
         "Original expectation: E = (1 - F)b + f. Zeros excluded from log ratio and counted in source TSV. Development-only sensitivity." else NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave(file.path(outdir, "three_cohort_quantitative_recovery_boxes.pdf"), p,
       width = 16, height = 9)
ggsave(file.path(outdir, "three_cohort_quantitative_recovery_boxes.png"), p,
       width = 16, height = 9, dpi = 250)
writeLines("DEVELOPMENT_ONLY", file.path(outdir, "DEVELOPMENT_ONLY.txt"))
writeLines("status\tPASS", file.path(outdir, "SUCCESS"))
message("[PASS] quantitative recovery boxes: ", outdir)
