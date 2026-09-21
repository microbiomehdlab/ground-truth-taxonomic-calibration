#!/usr/bin/env Rscript
# Descriptive observed-detection heatmap from corrected paired endpoints.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}
input <- value("--endpoints")
outdir <- value("--outdir")
if (dir.exists(outdir)) stop("Output directory already exists")
x <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("cohort", "condition", "analysis_population", "assembly_arm",
              "profiler", "target_label", "spike_fraction_target", "sample_id",
              "reference_type", "observed_abundance_fraction",
              "observed_detected_native_nonzero")
missing <- setdiff(required, names(x))
if (length(missing)) stop("Missing endpoint columns: ", paste(missing, collapse = ", "))
x <- x[x$analysis_population == "community" & x$assembly_arm == "original", ]
if (!nrow(x)) stop("No community original-assembly endpoints")
cohorts <- c("feng", "yachida", "zeller")
conditions <- c("Control", "Adenoma", "CRC")
profilers <- c("kraken2_bracken", "metaphlan4")
targets <- c("Bfrag", "Csym", "Dpne", "Fnuc", "Hhat", "Pmic",
             "Pana", "Psto", "Porp", "Pint")
doses <- c(.00001, .00005, .0001, .0005, .001, .005, .01, .05)
x$spike_fraction_target <- as.numeric(x$spike_fraction_target)
x$observed_abundance_fraction <- as.numeric(x$observed_abundance_fraction)
if (any(!is.finite(x$spike_fraction_target)) ||
    any(!is.finite(x$observed_abundance_fraction)) ||
    any(x$observed_abundance_fraction < 0)) stop("Invalid endpoint numbers")
if (any(!x$cohort %in% cohorts) || any(!x$condition %in% conditions) ||
    any(!x$profiler %in% profilers) || any(!x$target_label %in% targets))
  stop("Unexpected figure context")
# The complete design's dose grid is discovered from the corrected endpoint
# table rather than assumed from the older figure. Exact observed fractions
# are assigned to the nearest nominal level only within 5%.
nominal <- function(v) {
  hit <- which.min(abs(v - doses) / doses)
  if (abs(v - doses[hit]) / doses[hit] > .05) stop("Dose outside frozen grid: ", v)
  doses[hit]
}
x$dose <- vapply(x$spike_fraction_target, nominal, numeric(1))
expected_type <- ifelse(x$profiler == "metaphlan4", "genome_equivalent", "read_proportional")
if (any(x$reference_type != expected_type)) stop("Wrong endpoint reference type")
if (any(!x$observed_detected_native_nonzero %in% c(0, 1)) ||
    any(as.integer(x$observed_detected_native_nonzero) !=
        as.integer(x$observed_abundance_fraction > 0)))
  stop("Detection flag does not match native nonzero abundance")
key <- paste(x$cohort, x$condition, x$profiler, x$target_label,
             x$dose, x$sample_id, sep = "\037")
if (anyDuplicated(key)) stop("Duplicate sample/context")
context <- function(d) paste(d$cohort, d$condition, d$profiler,
                             d$target_label, d$dose, sep = "\037")
grid <- expand.grid(cohort = cohorts, condition = conditions, profiler = profilers,
                    target_label = targets, dose = sort(unique(x$dose)),
                    stringsAsFactors = FALSE)
if (length(unique(x$dose)) < 3L) stop("Fewer than three shared dose levels")
if (!setequal(context(x), context(grid))) stop("Incomplete detection context grid")
wilson <- function(k, n, z = 1.95996398454005) {
  p <- k / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, center - half), min(1, center + half))
}
summary <- do.call(rbind, lapply(split(x, context(x)), function(part) {
  n <- nrow(part)
  k <- sum(as.integer(part$observed_detected_native_nonzero))
  ci <- wilson(k, n)
  data.frame(cohort = part$cohort[1], condition = part$condition[1],
             profiler = part$profiler[1], target_label = part$target_label[1],
             dose = part$dose[1], detected = k, samples = n,
             prevalence = k / n, wilson_low = ci[1], wilson_high = ci[2])
}))
dir.create(outdir, recursive = TRUE)
write.table(summary, file.path(outdir, "detection_prevalence_wilson.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
hashes <- system2("sha256sum", c(shQuote(normalizePath(input)),
  shQuote(normalizePath(file.path(outdir, "detection_prevalence_wilson.tsv")))),
  stdout = TRUE, stderr = TRUE)
if (!is.null(attr(hashes, "status")) && attr(hashes, "status") != 0)
  stop("Could not hash detection sources")
writeLines(hashes, file.path(outdir, "source.sha256"))
summary$cohort <- factor(summary$cohort, levels = cohorts,
                         labels = c("Feng", "Yachida", "Zeller"))
summary$profiler <- factor(summary$profiler, levels = profilers,
                           labels = c("Kraken2 + Bracken", "MetaPhlAn 4"))
summary$condition <- factor(summary$condition, levels = conditions)
summary$target_label <- factor(summary$target_label, levels = rev(targets))
summary$dose_label <- factor(sprintf("%.3f%%", 100 * summary$dose),
                             levels = sprintf("%.3f%%", 100 * sort(unique(x$dose))))
p <- ggplot(summary, aes(dose_label, target_label, fill = prevalence)) +
  geom_tile(color = "white", linewidth = .25) +
  facet_grid(profiler + cohort ~ condition) +
  scale_fill_gradient(low = "#f2f5f9", high = "#156b91", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Observed detection of implanted targets",
       subtitle = "Community spikes; cell color is fraction of samples with positive native abundance",
       x = "Nominal implanted target fraction", y = "Implanted taxon",
       fill = "Detected") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(outdir, "three_cohort_detection_heatmap.pdf"), p,
       width = 18, height = 18)
ggsave(file.path(outdir, "three_cohort_detection_heatmap.png"), p,
       width = 18, height = 18, dpi = 220)
writeLines("DEVELOPMENT_ONLY", file.path(outdir, "DEVELOPMENT_ONLY.txt"))
writeLines("status\tPASS", file.path(outdir, "SUCCESS"))
message("[PASS] three-cohort detection heatmap: ", outdir)
