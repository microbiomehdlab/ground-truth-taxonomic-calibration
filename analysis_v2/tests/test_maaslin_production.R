source('analysis_v2/scripts/fit_maaslin_batch.R')
source('analysis_v2/lib/maaslin_contract.R')
source('analysis_v2/lib/maaslin_context.R')
set.seed(263)
scratch <- tempfile('maaslin-batch-fixture-'); dir.create(scratch)
input <- file.path(scratch,'input'); dir.create(input)
contexts <- data.frame(context_id=c('da3_0000000','da3_0000001'),cohort='feng',
    background='Adenoma',profiler='kraken2_bracken',n=5,arm=c('U','N'),
    anchor=c('0.0001',''),allocation_id='allocation_0000',stringsAsFactors=FALSE)
write.table(contexts,file.path(input,'contexts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
a <- matrix(runif(120,.00001,.01),10,dimnames=list(paste0('sample',1:10),paste0('taxon ',1:12)))
a[,11] <- 0; a[,12] <- rep(c(.0001,.001),each=5)
m <- data.frame(group=rep(0:1,each=5),row.names=rownames(a))
for (id in contexts$context_id) {
  folder <- file.path(input,id); dir.create(folder)
  write.table(data.frame(observation_id=rownames(a),a,check.names=FALSE),file.path(folder,'abundance.tsv'),
      sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(observation_id=rownames(m),m),file.path(folder,'metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(feature=colnames(a)[1:10],target_label=paste0('T',1:10)),
      file.path(folder,'targets.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
mock <- function(input_data,input_metadata,output,...) {
  opts <- list(...)
  stopifnot(opts$save_models,opts$cores==1,opts$analysis_method=='LM',
            opts$normalization=='NONE',opts$transform=='NONE',identical(opts$fixed_effects,'group'))
  dir.create(file.path(output,'fits'),recursive=TRUE)
  models <- lapply(input_data,function(y) lm(y~group,data=input_metadata))
  result <- do.call(rbind,lapply(seq_along(models),function(i) {
    s <- summary(models[[i]])$coefficients['group',]
    data.frame(feature=names(input_data)[i],metadata='group',coef=s[1],stderr=s[2],pval=s[4])
  }))
  result$qval <- p.adjust(result$pval,'BH')
  write.table(result,file.path(output,'all_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  saveRDS(models,file.path(output,'fits/models.rds'))
}
exercise <- function(name,backend) {
  out <- file.path(scratch,name)
  fit_maaslin_batch(input,out,getwd(),backend=backend)
  summary <- read.delim(file.path(out,'summary.tsv'),colClasses='character',na.strings=character())
  f <- read.delim(gzfile(file.path(out,'features.tsv.gz')),check.names=FALSE)
  t <- read.delim(gzfile(file.path(out,'targets.tsv.gz')),check.names=FALSE)
  stopifnot(nrow(summary)==2,nrow(f)==24,nrow(t)==20,
    identical(summary$anchor,c('0.0001','')),all(summary$mismatched_features=='0'),
    all(summary$family_n=='12'),all(summary$target_rows=='10'),
    !dir.exists(file.path(input,contexts$context_id[1],'native_fit')))
  reference <- fit_maaslin_context(a,m,colnames(a),file.path(scratch,paste0(name,'_reference')),backend=backend)
  first <- f[f$context_id==contexts$context_id[1],setdiff(names(f),'context_id')]
  stopifnot(isTRUE(all.equal(first,reference,check.attributes=FALSE,tolerance=1e-12)))
}
exercise('mock',mock)
cat('[PASS] Batched full-family statistics equal retained wrapper; blank/dose identity and cleanup preserved\n')
if (requireNamespace('Maaslin2',quietly=TRUE)) {
  exercise('actual',NULL)
  cat('[PASS] Actual pinned MaAsLin2 batch retention agrees with original wrapper\n')
} else {
  if(Sys.getenv('REQUIRE_MAASLIN_BATCH_BACKEND')=='1') stop('Actual pinned backend required')
  cat('[SKIP] Actual package unavailable locally; mandatory on cluster before arrays\n')
}
unlink(scratch,recursive=TRUE)
