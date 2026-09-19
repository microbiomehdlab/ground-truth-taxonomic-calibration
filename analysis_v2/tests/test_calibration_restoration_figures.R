#!/usr/bin/env Rscript
arg<-sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))]);root<-normalizePath(file.path(dirname(arg),"../.."))
t<-tempfile("calibration_figures.");dir.create(t);transfer<-file.path(t,"transfers","feng__to__zeller");dir.create(file.path(transfer,"restoration"),recursive=TRUE)
ledger<-data.frame(transfer="feng__to__zeller",training_cohorts="feng",validation_cohort="zeller",directory=transfer,active_coefficients=2,validation_profiles=4)
write.table(ledger,file.path(t,"transfer_ledger.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
c<-expand.grid(training_cohorts="feng",validation_cohort="zeller",analysis_population=c("community","independent"),profiler=c("kraken2_bracken","metaphlan4"),perturbation=c("Target","COMMUNITY_MIXTURE"),feature=c("Bystander A","Bystander B"))
c$slope_fraction_per_implanted_fraction=seq(.05,.2,length.out=nrow(c));c$observations=20;c$training_samples=12;c$positive_residual_fraction=.8;c$correction_active=1
write.table(c,file.path(transfer,"calibration_coefficients.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
a<-expand.grid(validation_cohort="zeller",analysis_population=c("community","independent"),profiler=c("kraken2_bracken","metaphlan4"),sample_id=paste0("s",1:3),spike_fraction_total=c(.001,.01))
a$training_cohorts="feng";a$source_profile=paste0("p",seq_len(nrow(a)));a$perturbation="fixture";a$features_evaluated=10;a$features_altered=2;a$features_improved=2;a$mae_before=.1;a$mae_after=.08;a$relative_mae_restoration=.2;a$direct_targets_protected=1
write.table(a,file.path(transfer,"abundance_restoration_metrics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
b<-expand.grid(validation_cohort="zeller",analysis_population=c("community","independent"),profiler=c("kraken2_bracken","metaphlan4"),spike_fraction_target=c(.001,.01))
b$cohort="zeller";b$study="study";b$target_label="target";b$assembly_arm="original";b$contrast="CRC_vs_Control";b$dose_level="dose_01";b$features_tested=10;b$baseline_biomarkers=3;b$uncorrected_biomarkers=2;b$corrected_biomarkers=3;b$baseline_biomarkers_rescued=1;b$baseline_biomarkers_harmed=0;b$induced_calls_removed=1;b$induced_calls_remaining=0;b$induced_calls_created=0;b$fraction_effect_errors_improved=.7;b$mean_absolute_effect_error_before=.2;b$mean_absolute_effect_error_after=.1;b$relative_mean_effect_restoration=.5
write.table(b,file.path(transfer,"restoration","calibration_restoration_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
out<-file.path(t,"report");rc<-system2("Rscript",c(file.path(root,"analysis_v2/scripts/plot_calibration_restoration.R"),"--calibration-root",t,"--outdir",out))
stopifnot(rc==0,file.exists(file.path(out,"SUCCESS")),file.exists(file.path(out,"figures/distortion_coefficient_atlas.pdf")),file.exists(file.path(out,"figures/abundance_error_restoration.pdf")),file.exists(file.path(out,"figures/biomarker_call_restoration.pdf")),file.exists(file.path(out,"figures/disease_effect_restoration.pdf")))
message("[PASS] calibration-restoration figure fixture")
