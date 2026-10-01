source('analysis_v2/lib/maaslin_contract.R')
source('analysis_v2/lib/maaslin_context.R')
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),'try-error'))
x <- matrix(c(.01,.02,.015,.05,.045,.06,rep(0,6)),nrow=6,
            dimnames=list(paste0('S',1:6),c('target','zero')))
meta <- data.frame(group=c(0,0,0,1,1,1),row.names=rownames(x))
mock <- function(input_data,input_metadata,output,...) {
  dir.create(output)
  models <- lapply(input_data,function(y) lm(y~group,data=input_metadata))
  rows <- do.call(rbind,lapply(seq_along(models),function(i) {
    s <- coef(summary(models[[i]]))['group',]
    data.frame(feature=names(input_data)[i],metadata='group',coef=s[1],stderr=s[2],pval=s[4],qval=s[4])
  }))
  write.table(rows,file.path(output,'all_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  saveRDS(models,file.path(output,'models.rds'))
}
root <- tempfile('maaslin_context_'); dir.create(root)
r <- fit_maaslin_context(x,meta,colnames(x),file.path(root,'good'),backend=mock)
stopifnot(r$estimable[1],!r$estimable[2],r$p_for_BH[2]==1,
          abs(r$wrapper_q[1]-min(1,2*r$raw_p[1]))<1e-12)
fails(fit_maaslin_context(x,meta,colnames(x),file.path(root,'good'),backend=mock))
warning_mock <- function(...) warning('fit diagnostic')
fails(fit_maaslin_context(x,meta,colnames(x),file.path(root,'warning'),backend=warning_mock))
stopifnot(!file.exists(file.path(root,'warning','SUCCESS')))
bad <- meta; bad$group[1]<-NA
fails(fit_maaslin_context(x,bad,colnames(x),file.path(root,'bad'),backend=mock))
stopifnot(!dir.exists(file.path(root,'bad')))
if (identical(Sys.getenv('REQUIRE_MAASLIN_BACKEND'),'1')) {
  stopifnot(requireNamespace('Maaslin2',quietly=TRUE),as.character(packageVersion('Maaslin2'))=='1.18.0')
  actual <- fit_maaslin_context(x,meta,colnames(x),file.path(root,'actual'))
  stopifnot(abs(actual$beta[1]-r$beta[1])<1e-8,abs(actual$raw_p[1]-r$raw_p[1])<1e-8)
  set.seed(41001)
  ids <- rep(paste0('person',1:12),each=2)
  g <- rep(0:1,12)
  y <- 10+rep(rnorm(12,0,1.5),each=2)+.7*g+rnorm(24,0,.2)
  pairs <- matrix(1e-8*(2^y-1),ncol=1,
                  dimnames=list(paste0('obs',1:24),'paired_target'))
  pm <- data.frame(group=g,biological_sample_id=ids,
                   spike_state=ifelse(g==0,'original','spiked'),row.names=rownames(pairs))
  paired <- fit_maaslin_context(pairs,pm,colnames(pairs),file.path(root,'actual_paired'),paired=TRUE)
  direct <- lmerTest::lmer(y~group+(1|biological_sample_id),data=pm)
  ds <- coef(summary(direct))['group',]
  stopifnot(abs(paired$beta-ds['Estimate'])<1e-7,
            abs(paired$raw_p-ds['Pr(>|t|)'])<1e-7)
  cat('[PASS] Pinned MaAsLin2 backend agrees with direct LM fixture\n')
  cat('[PASS] Paired MaAsLin2 backend agrees with direct mixed-model fixture\n')
} else cat('[NOTE] Actual MaAsLin2 backend NOT tested; mock backend only\n')
cat('[PASS] Context wrapper contract fixtures\n')
