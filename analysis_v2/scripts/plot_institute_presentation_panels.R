#!/usr/bin/env Rscript
# Readable, development-only slides from frozen baseline calls and corrected endpoints.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Required argument: ", flag)
  args[[i + 1L]]
}
baseline_path <- value("--baseline")
endpoints_path <- value("--endpoints")
outdir <- value("--outdir")
if (dir.exists(outdir)) stop("Output directory already exists: ", outdir)
baseline <- read.delim(baseline_path, check.names = FALSE, stringsAsFactors = FALSE)
endpoints <- read.delim(endpoints_path, check.names = FALSE, stringsAsFactors = FALSE)
needed_baseline <- c("cohort", "sample_id", "condition", "profiler", "target_label",
                     "abundance_fraction")
needed_endpoints <- c("cohort", "sample_id", "condition", "analysis_population",
                      "assembly_arm", "profiler", "target_label", "spike_fraction_target",
                      "reference_type", "observed_abundance_fraction",
                      "expected_abundance_profiler_scale")
if (length(setdiff(needed_baseline, names(baseline))) ||
    length(setdiff(needed_endpoints, names(endpoints)))) stop("Required figure columns missing")
cohorts <- c("feng", "yachida", "zeller")
profilers <- c("kraken2_bracken", "metaphlan4")
conditions <- c("Control", "Adenoma", "CRC")
if (!setequal(unique(baseline$cohort), cohorts) ||
    !setequal(unique(baseline$profiler), profilers)) stop("Incomplete baseline strata")
if (any(!baseline$condition %in% conditions)) stop("Unexpected baseline condition")
baseline$abundance_fraction <- as.numeric(baseline$abundance_fraction)
if (any(!is.finite(baseline$abundance_fraction) |
        baseline$abundance_fraction < 0 | baseline$abundance_fraction > 1))
  stop("Invalid baseline abundance")
endpoints <- endpoints[endpoints$analysis_population == "community" &
                       endpoints$assembly_arm == "original", , drop = FALSE]
if (!nrow(endpoints) || !setequal(unique(endpoints$cohort), cohorts) ||
    !setequal(unique(endpoints$profiler), profilers)) stop("Incomplete community endpoints")
expected_ref <- ifelse(endpoints$profiler == "metaphlan4",
                       "genome_equivalent", "read_proportional")
if (any(endpoints$reference_type != expected_ref)) stop("Wrong endpoint reference type")
for (field in c("spike_fraction_target", "observed_abundance_fraction",
                "expected_abundance_profiler_scale")) {
  endpoints[[field]] <- as.numeric(endpoints[[field]])
  if (any(!is.finite(endpoints[[field]]))) stop("Nonfinite ", field)
}
if (any(endpoints$observed_abundance_fraction < 0) ||
    any(endpoints$expected_abundance_profiler_scale <= 0))
  stop("Invalid observed or expected abundance")
if (any(!endpoints$condition %in% conditions)) stop("Unexpected endpoint condition")
dose_grid <- c(.0001, .0005, .001)
near_dose <- function(v) {
  closest <- which.min(abs(dose_grid - v) / dose_grid)
  if (abs(dose_grid[[closest]] - v) / dose_grid[[closest]] > .05) return(NA_real_)
  dose_grid[[closest]]
}
endpoints$dose <- vapply(endpoints$spike_fraction_target, near_dose, numeric(1))
endpoints <- endpoints[is.finite(endpoints$dose), , drop = FALSE]
if (!setequal(unique(endpoints$dose), dose_grid)) stop("Weak-dose grid incomplete")
endpoints$dose_label <- factor(sprintf("%.2f%%", 100 * endpoints$dose),
                               levels = sprintf("%.2f%%", 100 * dose_grid))
dir.create(outdir, recursive = TRUE)
dir.create(file.path(outdir, "sources"))
write_source <- function(data, filename) {
  write.table(data, file.path(outdir, "sources", filename), sep = "\t",
              row.names = FALSE, quote = FALSE, na = "NA")
}
save_slide <- function(plot, stem, width = 14, height = 7.5) {
  ggsave(file.path(outdir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 250, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".pdf")), plot,
         width = width, height = height)
}
wilson <- function(k, n, z = qnorm(.975)) {
  p <- k / n
  d <- 1 + z^2 / n
  center <- (p + z^2 / (2*n)) / d
  radius <- z * sqrt(p*(1-p)/n + z^2/(4*n*n)) / d
  c(max(0, center-radius), min(1, center+radius))
}
colors <- c("Kraken2 + Bracken" = "#009E73", "MetaPhlAn 4" = "#6655CC")
cohort_colors <- c(Feng = "#E46F50", Yachida = "#208E80", Zeller = "#3A759D")
profiler_label <- function(x) factor(x, levels = profilers,
                                     labels = names(colors))
cohort_label <- function(x) factor(x, levels = cohorts,
                                   labels = names(cohort_colors))
theme_slide <- theme_bw(base_size = 17) + theme(
  panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey94"),
  strip.text = element_text(face = "bold"), legend.position = "bottom",
  plot.title = element_text(face = "bold"))

# Slide 3: baseline prevalence, restricted to two interpretable taxa.
b <- baseline[baseline$target_label %in% c("Bfrag", "Fnuc"), , drop = FALSE]
key <- paste(b$cohort, b$sample_id, b$condition, b$profiler, b$target_label, sep = "\037")
if (anyDuplicated(key)) stop("Duplicate baseline biological sample")
base_parts <- split(b, interaction(b$cohort, b$condition, b$profiler,
                                   b$target_label, drop = TRUE))
bs <- do.call(rbind, lapply(base_parts, function(part) {
  n <- nrow(part); k <- sum(part$abundance_fraction > 0)
  ci <- wilson(k, n)
  data.frame(cohort = part$cohort[1], condition = part$condition[1],
             profiler = part$profiler[1], target_label = part$target_label[1],
             positive = k, samples = n, prevalence = k/n,
             lower_95 = ci[1], upper_95 = ci[2])
}))
if (nrow(bs) != 3*3*2*2) stop("Baseline slide contexts incomplete")
write_source(bs, "slide3_baseline_prevalence.tsv")
bs$cohort <- cohort_label(bs$cohort)
bs$profiler <- profiler_label(bs$profiler)
bs$condition <- factor(bs$condition, levels = conditions)
bs$target_label <- factor(bs$target_label, levels = c("Bfrag", "Fnuc"),
                          labels = c("B. fragilis", "F. nucleatum"))
p3 <- ggplot(bs, aes(condition, prevalence, fill = profiler)) +
  geom_col(position = position_dodge(width = .7), width = .62) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                position = position_dodge(width = .7), width = .13) +
  facet_grid(cohort ~ target_label) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(round(100*v), "%")) +
  labs(title = "Profilers disagree before any signal is added",
       subtitle = "Unspiked samples; bars show detected fraction, lines show 95% Wilson intervals",
       x = NULL, y = "Samples with reported species", fill = NULL) + theme_slide
save_slide(p3, "slide3_baseline_visibility")

# Slide 4: corrected quantitative recovery, three contrasting taxa.
r <- endpoints[endpoints$target_label %in% c("Bfrag", "Dpne", "Fnuc"), , drop = FALSE]
key <- paste(r$cohort, r$sample_id, r$condition, r$profiler,
             r$target_label, r$dose, sep = "\037")
if (anyDuplicated(key)) stop("Duplicate recovery sample/context")
r$recovery <- ifelse(r$observed_abundance_fraction > 0,
                     log2(r$observed_abundance_fraction /
                          r$expected_abundance_profiler_scale), NA_real_)
parts <- split(r, interaction(r$cohort, r$profiler, r$target_label, r$dose, drop = TRUE))
rs <- do.call(rbind, lapply(parts, function(part) {
  positive <- part$recovery[is.finite(part$recovery)]
  data.frame(cohort = part$cohort[1], profiler = part$profiler[1],
             target_label = part$target_label[1], dose = part$dose[1],
             n_total = nrow(part), n_positive = length(positive),
             n_zero = nrow(part) - length(positive),
             q1 = if (length(positive)) unname(quantile(positive, .25)) else NA_real_,
             median = if (length(positive)) median(positive) else NA_real_,
             q3 = if (length(positive)) unname(quantile(positive, .75)) else NA_real_)
}))
if (nrow(rs) != 3*2*3*3) stop("Recovery slide contexts incomplete")
write_source(rs, "slide4_recovery_iqr.tsv")
rs$cohort <- cohort_label(rs$cohort)
rs$profiler <- profiler_label(rs$profiler)
rs$target_label <- factor(rs$target_label, levels = c("Bfrag", "Dpne", "Fnuc"),
                          labels = c("B. fragilis", "D. pneumoniae", "F. nucleatum"))
rs$dose_label <- factor(sprintf("%.2f%%", 100 * rs$dose),
                        levels = sprintf("%.2f%%", 100 * dose_grid))
p4 <- ggplot(rs, aes(dose_label, median, color = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                position = position_dodge(width = .45), width = .13, linewidth = .7,
                na.rm = TRUE) +
  geom_point(position = position_dodge(width = .45), size = 3, na.rm = TRUE) +
  facet_grid(profiler ~ target_label) +
  scale_color_manual(values = cohort_colors) +
  labs(title = "Weak-signal recovery depends on taxon and profiler",
       subtitle = "Community spikes; median and sample IQR among positive reports; zeros counted separately",
       x = "Implanted target fraction", y = expression(log[2](observed/expected)),
       color = "Cohort") + theme_slide
save_slide(p4, "slide4_focused_recovery")

# Slide 5: one complete cohort; detection is not conflated with accuracy.
f <- endpoints[endpoints$cohort == "yachida" &
               endpoints$target_label == "Fnuc", , drop = FALSE]
key <- paste(f$sample_id, f$condition, f$profiler, f$dose, sep = "\037")
if (anyDuplicated(key)) stop("Duplicate F. nucleatum sample/context")
fparts <- split(f, interaction(f$condition, f$profiler, f$dose, drop = TRUE))
fs <- do.call(rbind, lapply(fparts, function(part) {
  n <- nrow(part); k <- sum(part$observed_abundance_fraction > 0)
  ci <- wilson(k, n)
  data.frame(condition = part$condition[1], profiler = part$profiler[1],
             dose = part$dose[1], positive = k, samples = n,
             prevalence = k/n, lower_95 = ci[1], upper_95 = ci[2])
}))
if (nrow(fs) != 3*2*3) stop("F. nucleatum slide contexts incomplete")
write_source(fs, "slide5_fnuc_detection.tsv")
fs$condition <- factor(fs$condition, levels = conditions)
fs$profiler <- profiler_label(fs$profiler)
fs$dose_label <- factor(sprintf("%.2f%%", 100 * fs$dose),
                        levels = sprintf("%.2f%%", 100 * dose_grid))
p5 <- ggplot(fs, aes(dose_label, prevalence, color = profiler, group = profiler)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = .10, linewidth = .65) +
  facet_wrap(~condition, nrow = 1) +
  scale_color_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(round(100*v), "%")) +
  labs(title = "Weak F. nucleatum additions are reported differently by profiler",
       subtitle = "Yachida community samples; positive reported abundance and 95% Wilson intervals",
       x = "Implanted target fraction", y = "Samples reporting F. nucleatum",
       color = NULL) + theme_slide
save_slide(p5, "slide5_fnuc_detection")

inputs <- c(baseline_path, endpoints_path)
outputs <- list.files(outdir, pattern = "\\.(tsv|png|pdf)$", recursive = TRUE,
                      full.names = TRUE)
hashes <- system2("sha256sum", shQuote(normalizePath(c(inputs, outputs))),
                  stdout = TRUE, stderr = TRUE)
if (!is.null(attr(hashes, "status")) && attr(hashes, "status") != 0)
  stop("Could not checksum figure inputs and outputs")
writeLines(hashes, file.path(outdir, "source.sha256"))
writeLines("status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO",
           file.path(outdir, "DEVELOPMENT_ONLY.txt"))
writeLines("status\tPASS", file.path(outdir, "SUCCESS"))
message("[PASS] Three focused presentation plots: ", outdir)
