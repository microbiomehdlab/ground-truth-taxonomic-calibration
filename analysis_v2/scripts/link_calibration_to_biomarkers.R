#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args); if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
endpoints_path <- value("--endpoints"); metrics_path <- value("--biomarker-metrics")
secondary_metrics_path <- value("--secondary-biomarker-metrics")
outdir <- value("--outdir"); analysis_status <- value("--analysis-status")
if (is.null(endpoints_path) || is.null(metrics_path) || is.null(outdir) ||
    !analysis_status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE"))
  stop("Required: --endpoints FILE --biomarker-metrics FILE --outdir DIR --analysis-status DEVELOPMENT_ONLY|DEFINITIVE")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("tables", "figure_source", "figures", "diagnostics", "provenance"))
  dir.create(file.path(outdir, subdir), showWarnings = FALSE)

endpoints <- read.delim(endpoints_path, check.names=FALSE, stringsAsFactors=FALSE)
metric_paths <- c(metrics_path, if(!is.null(secondary_metrics_path)) secondary_metrics_path)
metrics <- do.call(rbind,lapply(metric_paths,function(path)
  read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings="NA")))
endpoint_required <- c("cohort", "study", "sample_id", "condition", "analysis_population",
  "target_label", "assembly_arm", "profiler", "spike_fraction_target", "response_ratio",
  "recovered_spike_signal", "signed_reference_error", "absolute_reference_error")
metric_required <- c("cohort", "study", "analysis_population", "target_label", "assembly_arm",
  "profiler", "contrast", "spike_fraction_target", "q_threshold", "target_called",
  "target_effect", "target_q_value", "precision", "off_target_enriched_calls")
if (length(setdiff(endpoint_required, names(endpoints)))) stop("Endpoint columns missing.")
if (length(setdiff(metric_required, names(metrics)))) stop("Biomarker metric columns missing.")
for (field in c("spike_fraction_target", "response_ratio", "recovered_spike_signal",
                "signed_reference_error", "absolute_reference_error")) endpoints[[field]] <- as.numeric(endpoints[[field]])
for (field in c("spike_fraction_target", "q_threshold", "target_called", "target_effect",
                "target_q_value", "precision", "off_target_enriched_calls")) metrics[[field]] <- as.numeric(metrics[[field]])
if (!nrow(endpoints) || !nrow(metrics) || anyNA(endpoints[c("spike_fraction_target", "response_ratio")]) ||
    anyNA(metrics[c("spike_fraction_target", "q_threshold", "target_called")])) stop("Invalid linkage inputs.")

nominal <- c(.0001, .0005, .001, .005, .01, .05)
map_dose <- function(x) {
  mapped <- vapply(x, function(d) nominal[which.min(abs(nominal-d))], numeric(1))
  if (any(abs(x-mapped) > pmax(1e-10, mapped*.001))) stop("Dose outside frozen nominal grid.")
  mapped
}
endpoints$dose_fraction_nominal <- map_dose(endpoints$spike_fraction_target)
metrics$dose_fraction_nominal <- map_dose(metrics$spike_fraction_target)
endpoints$dose_rank <- match(endpoints$dose_fraction_nominal, nominal)
metrics$dose_rank <- match(metrics$dose_fraction_nominal, nominal)

endpoint_key <- c("cohort", "study", "analysis_population", "condition", "target_label",
                  "assembly_arm", "profiler", "dose_rank", "dose_fraction_nominal")
endpoint_groups <- split(endpoints, interaction(endpoints[endpoint_key], drop=TRUE, lex.order=TRUE))
quantitative <- do.call(rbind, lapply(endpoint_groups, function(x) data.frame(
  x[1, endpoint_key, drop=FALSE], biological_samples=nrow(x),
  median_achieved_fraction=median(x$spike_fraction_target),
  median_response_ratio=median(x$response_ratio), mean_response_ratio=mean(x$response_ratio),
  response_ratio_iqr=IQR(x$response_ratio),
  median_recovered_spike_signal=median(x$recovered_spike_signal),
  median_signed_reference_error=median(x$signed_reference_error),
  median_absolute_reference_error=median(x$absolute_reference_error), stringsAsFactors=FALSE)))
quantitative$calibration_class <- ifelse(quantitative$median_response_ratio < .8, "under_response",
  ifelse(quantitative$median_response_ratio > 1.2, "over_response", "read_proportional_band"))

pooled_key <- setdiff(endpoint_key,"condition")
pooled_groups <- split(endpoints,interaction(endpoints[pooled_key],drop=TRUE,lex.order=TRUE))
pooled_quantitative <- do.call(rbind,lapply(pooled_groups,function(x) data.frame(
  x[1,pooled_key,drop=FALSE],condition="ALL",biological_samples=nrow(x),
  median_achieved_fraction=median(x$spike_fraction_target),median_response_ratio=median(x$response_ratio),
  mean_response_ratio=mean(x$response_ratio),response_ratio_iqr=IQR(x$response_ratio),
  median_recovered_spike_signal=median(x$recovered_spike_signal),
  median_signed_reference_error=median(x$signed_reference_error),
  median_absolute_reference_error=median(x$absolute_reference_error),stringsAsFactors=FALSE)))
pooled_quantitative$calibration_class <- ifelse(pooled_quantitative$median_response_ratio<.8,"under_response",
  ifelse(pooled_quantitative$median_response_ratio>1.2,"over_response","read_proportional_band"))
quantitative$analysis_scope <- "phenotype_stratified_secondary"
pooled_quantitative$analysis_scope <- "pooled_primary"
pooled_quantitative <- pooled_quantitative[names(quantitative)]
quantitative <- rbind(quantitative,pooled_quantitative)
metrics$analysis_scope <- ifelse(metrics$contrast=="spiked_vs_matched_baseline__pooled",
  "pooled_primary","phenotype_stratified_secondary")
metrics$condition <- ifelse(metrics$analysis_scope=="pooled_primary","ALL",
  sub("^spiked_vs_matched_baseline__background_", "", metrics$contrast))
if (any(metrics$analysis_scope=="phenotype_stratified_secondary" & metrics$condition==metrics$contrast))
  stop("Could not recover phenotype background from contrast.")
quantitative <- quantitative[quantitative$analysis_scope %in% unique(metrics$analysis_scope),,drop=FALSE]
join_key <- c("cohort", "study", "analysis_population", "condition", "target_label",
              "assembly_arm", "profiler", "dose_rank", "dose_fraction_nominal", "analysis_scope")
if (anyDuplicated(quantitative[join_key])) stop("Duplicate quantitative linkage context.")
if (anyDuplicated(metrics[c(join_key, "q_threshold")])) stop("Duplicate biomarker linkage context.")
linked <- merge(metrics, quantitative, by=join_key, all=TRUE, sort=FALSE)
if (nrow(linked) != nrow(metrics) || anyNA(linked$biological_samples))
  stop("Quantitative and biomarker contexts do not form a complete many-to-one join.")
linked$dose_percent_nominal <- 100*linked$dose_fraction_nominal
order_fields <- c(join_key, "q_threshold")
linked <- linked[do.call(order, linked[order_fields]), ]
write.table(linked, file.path(outdir, "figure_source", "calibration_biomarker_linkage.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="NA")

primary <- linked[abs(linked$q_threshold-.05)<1e-12, ]
if(any(primary$analysis_scope=="pooled_primary")) primary<-primary[primary$analysis_scope=="pooled_primary",]
association_key <- c("analysis_scope", "cohort", "analysis_population", "assembly_arm", "profiler")
association_groups <- split(primary, interaction(primary[association_key], drop=TRUE, lex.order=TRUE))
safe_cor <- function(x,y) if (length(x) >= 3 && sd(x)>0 && sd(y)>0) suppressWarnings(cor(x,y,method="spearman")) else NA_real_
associations <- do.call(rbind, lapply(association_groups, function(x) data.frame(
  x[1, association_key, drop=FALSE], contexts=nrow(x),
  spearman_response_ratio_vs_target_effect=safe_cor(x$median_response_ratio,x$target_effect),
  spearman_abs_calibration_error_vs_off_target_calls=safe_cor(abs(x$median_response_ratio-1),x$off_target_enriched_calls),
  target_recall=mean(x$target_called), median_response_ratio=median(x$median_response_ratio),
  median_off_target_calls=median(x$off_target_enriched_calls), stringsAsFactors=FALSE)))
write.table(associations, file.path(outdir,"tables","calibration_biomarker_associations.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE,na="NA")

class_key <- c("analysis_scope","cohort","analysis_population","assembly_arm","profiler","calibration_class")
class_groups <- split(primary, interaction(primary[class_key],drop=TRUE,lex.order=TRUE))
class_summary <- do.call(rbind,lapply(class_groups,function(x) data.frame(
  x[1,class_key,drop=FALSE], contexts=nrow(x), target_recall=mean(x$target_called),
  mean_precision=mean(x$precision), median_target_effect=median(x$target_effect),
  context_sum_off_target_calls=sum(x$off_target_enriched_calls), stringsAsFactors=FALSE)))
write.table(class_summary,file.path(outdir,"tables","biomarker_outcomes_by_calibration_class.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE,na="NA")

labels <- c(kraken2_bracken="Kraken2 + Bracken",metaphlan4="MetaPhlAn 4")
primary$profiler_display <- unname(labels[primary$profiler])
theme_report <- theme_bw(base_size=10)+theme(legend.position="bottom",panel.grid.minor=element_blank())
p1 <- ggplot(primary,aes(median_response_ratio,target_effect,color=profiler_display,shape=factor(target_called)))+
  geom_vline(xintercept=1,linetype=2,color="grey50")+geom_point(size=2)+
  facet_grid(target_label~assembly_arm)+labs(x="Median read-perturbation response ratio",
  y="Paired artificial-target log2 effect",color="Profiler",shape="Target called")+theme_report
p2 <- ggplot(primary,aes(abs(median_response_ratio-1),off_target_enriched_calls,color=profiler_display))+
  geom_point(size=2)+facet_grid(target_label~assembly_arm)+labs(x="Absolute deviation of median response ratio from 1",
  y="Off-target enriched calls",color="Profiler")+theme_report
p3 <- ggplot(primary,aes(dose_percent_nominal,median_response_ratio,color=profiler_display,
  group=interaction(profiler_display,condition)))+geom_hline(yintercept=1,linetype=2,color="grey50")+
  geom_line(alpha=.5)+geom_point(size=1)+facet_grid(target_label~assembly_arm)+
  labs(x="Implanted target fraction (%)",y="Median read-perturbation response ratio",color="Profiler")+theme_report
for (item in list(list("calibration_vs_target_effect",p1),list("calibration_error_vs_off_target_burden",p2),
                  list("response_ratio_dose_response",p3))) {
  name<-item[[1]]; plot<-item[[2]]
  ggsave(file.path(outdir,"figures",paste0(name,".pdf")),plot,width=8,height=6,units="in")
  ggsave(file.path(outdir,"figures",paste0(name,".png")),plot,width=8,height=6,units="in",dpi=300)
}
diagnostics <- data.frame(metric=c("endpoint_rows","biomarker_metric_rows","linked_rows","primary_rows",
  "quantitative_contexts","cohorts","targets","profilers","analysis_status"),value=c(nrow(endpoints),
  nrow(metrics),nrow(linked),nrow(primary),nrow(quantitative),length(unique(linked$cohort)),
  length(unique(linked$target_label)),length(unique(linked$profiler)),analysis_status))
write.table(diagnostics,file.path(outdir,"diagnostics","linkage_diagnostics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("# Draft figure captions","",
  "## Calibration versus artificial-target effect","Context-level paired biomarker effect against median read-perturbation response ratio. The ratio evaluates response to implanted sequencing evidence on each profiler's native scale; it is not cellular abundance accuracy.","",
  "## Calibration error versus off-target burden","Off-target enriched calls against absolute deviation from read-proportional response. Points are descriptive experimental contexts, not independent biological replicates.","",
  "## Response-ratio dose response","Median response ratio across biological samples at each frozen nominal implanted fraction. Exact achieved fractions are retained in the source table."),file.path(outdir,"captions.md"))
manifest <- data.frame(field=c("status","analysis","created_at","endpoints","biomarker_metrics","secondary_biomarker_metrics"),
 value=c(analysis_status,"calibration_to_artificial_biomarker_linkage",format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
 normalizePath(endpoints_path),normalizePath(metrics_path),if(is.null(secondary_metrics_path)) "NONE" else normalizePath(secondary_metrics_path)))
write.table(manifest,file.path(outdir,"provenance","run_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
if (analysis_status=="DEVELOPMENT_ONLY") writeLines(c("status\tDEVELOPMENT_ONLY","use_for_manuscript\tNO"),file.path(outdir,"DEVELOPMENT_ONLY.txt"))
files <- c(normalizePath(endpoints_path),normalizePath(metric_paths),list.files(outdir,recursive=TRUE,full.names=TRUE))
status <- system2("sha256sum",files,stdout=file.path(outdir,"provenance","linkage.sha256"))
if (!identical(status,0L)) stop("Could not seal linkage package.")
writeLines(c("analysis\tcalibration_to_artificial_biomarker_linkage",paste0("analysis_status\t",analysis_status),"status\tPASS"),file.path(outdir,"SUCCESS"))
message("[PASS] Calibration-to-biomarker linkage completed: ",outdir)
