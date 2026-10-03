source('analysis_v2/scripts/time_maaslin_context.R')
set.seed(746)
scratch <- tempfile('maaslin-timing-test-'); dir.create(scratch)
a <- matrix(runif(200,1e-6,.01),20,dimnames=list(paste0('s',1:20),paste0('taxon ',1:10)))
a[,9] <- 0; a[,10] <- rep(c(.0001,.001),each=10)
meta <- data.frame(observation_id=rownames(a),group=rep(0:1,each=10))
fixture <- function(name) {
  folder <- file.path(scratch,name); dir.create(folder)
  write.table(data.frame(observation_id=rownames(a),a,check.names=FALSE),file.path(folder,'abundance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(meta,file.path(folder,'metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  folder
}
# This double tests wrapper/files/comparison, not availability of actual MaAsLin2.
backend <- function(input_data,input_metadata,output,...) {
  args <- list(...)
  stopifnot(args$normalization=='NONE',args$transform=='NONE',args$analysis_method=='LM',
    identical(args$fixed_effects,'group'),is.null(args$random_effects),args$cores==1,args$save_models,
    args$min_abundance==0,args$min_prevalence==0,args$min_variance==0,!args$standardize)
  dir.create(file.path(output,'fits'),recursive=TRUE)
  models <- lapply(input_data,function(y) lm(y~group,data=input_metadata))
  rows <- do.call(rbind,lapply(seq_along(models),function(i) {
    coef <- summary(models[[i]])$coefficients
    data.frame(feature=names(input_data)[i],metadata='group',coef=coef[2,1],stderr=coef[2,2],pval=coef[2,4])
  }))
  rows$qval <- p.adjust(rows$pval,'BH')
  write.table(rows,file.path(output,'all_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  saveRDS(models,file.path(output,'fits','models.rds'))
}
folder <- fixture('double')
comparison <- time_maaslin_context(folder,getwd(),backend=backend)
stopifnot(nrow(comparison)==10,all(comparison$agreement),
  comparison$actual_status[9]=='NON_ESTIMABLE_CONSTANT',
  comparison$actual_status[10]=='NON_ESTIMABLE_PERFECT_FIXED_FIT')
actual <- read.delim(file.path(folder,'direct_maaslin','context_results.tsv'),check.names=FALSE)
fast <- unpaired_null_results(log2(1+a/1e-8),meta$group)
actual$beta[1] <- actual$beta[1]+1
stopifnot(!compare_maaslin_timing(actual,fast)$agreement[1])
actual$beta[1] <- fast$beta[1]
actual$positive_discovery[1] <- !actual$positive_discovery[1]
stopifnot(!compare_maaslin_timing(actual,fast)$discovery_equal[1])
stopifnot(inherits(try(compare_maaslin_timing(actual[-1,],fast),silent=TRUE),'try-error'))
cat('[PASS] Timing/comparison fixture; changes in coefficients, decisions and identities detected\n')
if (requireNamespace('Maaslin2',quietly=TRUE)) {
  comparison <- time_maaslin_context(fixture('actual'),getwd())
  stopifnot(all(comparison$agreement))
  cat('[PASS] Actual pinned MaAsLin2 timing/equivalence fixture\n')
} else {
  if (Sys.getenv('REQUIRE_MAASLIN_TIMING_BACKEND')=='1') stop('Pinned MaAsLin2 backend required')
  cat('[SKIP] MaAsLin2 unavailable locally; mandatory real-backend test in cluster prepare job\n')
}
unlink(scratch,recursive=TRUE)
