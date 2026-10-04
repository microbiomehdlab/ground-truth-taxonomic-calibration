#!/usr/bin/env Rscript
repo <- normalizePath('.')
source(file.path(repo,'analysis_v2/scripts/fit_clinical_spike.R'))
root <- tempfile(); dir.create(root)
set.seed(312)
n <- 36
m <- data.frame(observation_id=paste0('person',1:n),biological_sample_id=paste0('person',1:n),
                group=rep(0:1,each=n/2),spike_state=rep(c('original','spiked'),each=n/2),
                age=sample(30:75,n,TRUE),sex=sample(c('Female','Male'),n,TRUE))
x <- data.frame(observation_id=m$observation_id,
                target=2^(rnorm(n,10,1)+2*m$group)*1e-8,
                background=2^rnorm(n,10,1)*1e-8,constant=0)
write.table(m,file.path(root,'metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(x,file.path(root,'abundance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
backend <- function(input_data,input_metadata,output,...) {
  dir.create(file.path(output,'fits'),recursive=TRUE)
  models <- lapply(input_data,function(y) lm(y~group+age+sex,data=input_metadata))
  result <- do.call(rbind,lapply(names(models),function(name) {
    s <- coef(summary(models[[name]]))['group',]
    data.frame(feature=name,metadata='group',coef=s[1],stderr=s[2],pval=s[4],qval=s[4])
  }))
  saveRDS(models,file.path(output,'fits/models.rds'))
  write.table(result,file.path(output,'all_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
r <- fit_clinical_spike(root,repo,backend)
stopifnot(nrow(r)==3,r$positive_discovery[1],!r$estimable[3],
          r$p_for_BH[3]==1,all(r$residual_df==32),
          all(r$control_n==18),all(r$adenoma_spiked_n==18),
          identical(r$wrapper_q,p.adjust(r$p_for_BH,'BH')))
cat('[PASS] Clinical covariate LM, complete family BH, constant feature and group summaries\n')
# Optional cluster test uses the actual pinned backend, never an automatic fallback.
if (Sys.getenv('REQUIRE_CLINICAL_SPIKE_BACKEND')=='1') {
  actual <- tempfile(); dir.create(actual)
  file.copy(file.path(root,c('metadata.tsv','abundance.tsv')),actual)
  r2 <- fit_clinical_spike(actual,repo)
  stopifnot(max(abs(r2$beta[r2$estimable]-r$beta[r$estimable]))<1e-7,
            max(abs(r2$wrapper_q-r$wrapper_q))<1e-7)
  cat('[PASS] Actual pinned MaAsLin2 clinical spike backend\n')
}
