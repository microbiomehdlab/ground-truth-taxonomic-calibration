source('analysis_v2/lib/maaslin_contract.R')
source('analysis_v2/lib/maaslin_context.R')
source('analysis_v2/scripts/audit_da1_context.R')
write_tsv <- function(x,p) write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
fixture_backend <- function(input_data,input_metadata,output,...) {
  dir.create(output)
  models <- lapply(input_data,function(y) lm(y~group+age+sex,data=input_metadata))
  rows <- do.call(rbind,lapply(seq_along(models),function(i) {
    s <- coef(summary(models[[i]]))['group',]
    data.frame(feature=names(input_data)[i],metadata='group',coef=s[1],stderr=s[2],pval=s[4],qval=s[4])
  }))
  write_tsv(rows,file.path(output,'all_results.tsv'))
  dir.create(file.path(output,'fits')); saveRDS(models,file.path(output,'fits/models.rds'))
}
fit_fixture <- function(task) {
  m <- read.delim(file.path(task,'metadata.tsv'),row.names=1)
  m$sex <- factor(m$sex,levels=c('Female','Male'))
  x <- as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE))
  actual <- identical(Sys.getenv('REQUIRE_MAASLIN_BACKEND'),'1')
  fit_maaslin_context(x,m,colnames(x),file.path(task,'fit'),covariates=c('age','sex'),
                     backend=if(actual) NULL else fixture_backend)
  writeLines(if(actual) capture.output(sessionInfo()) else
      'TEST FIXTURE ONLY: Maaslin2_1.18.0',file.path(task,'sessionInfo.txt'))
}
args <- commandArgs(TRUE)
if (length(args)==2 && args[1]=='--fixtures') {
  for (task in list.dirs(args[2],recursive=FALSE)) fit_fixture(task)
} else {
  root <- tempfile('da1-export-'); dir.create(root)
  set.seed(712)
  m <- data.frame(observation_id=paste0('o',1:24),biological_sample_id=paste0('p',1:24),
      group=rep(0:1,each=12),spike_state='original',age=sample(30:75,24),sex=rep(c('Female','Male'),12))
  v <- 6 + 1.5*m$group + .02*m$age + rnorm(24)
  x <- data.frame(observation_id=m$observation_id,signal=1e-8*(2^v-1),
      lower=1e-8*(2^(12-v/2)-1),constant=0,perfect=1e-8*(2^(2+2*m$group)-1))
  write_tsv(m,file.path(root,'metadata.tsv')); write_tsv(x,file.path(root,'abundance.tsv'))
  write_tsv(data.frame(context_id='fixture',analysis='DA1',paired=0,covariates='age,sex',
      population='community',nominal_total_dose=0,observations_n=24,family_n=4),file.path(root,'context.tsv'))
  fit_fixture(root)
  output <- file.path(root,'audit.tsv'); audit_da1_context(root,output)
  r <- read.delim(output)
  stopifnot(all(r$residual_df==20),r$beta[1]>0,r$beta[2]<0,
      r$status[3]=='NON_ESTIMABLE_CONSTANT',r$status[4]=='NON_ESTIMABLE_PERFECT_FIXED_FIT',
      all(is.na(r$ci95_low[3:4])),all(r$control_n==12),all(r$disease_n==12))
  m$sex <- factor(m$sex,levels=c('Female','Male'))
  fit <- lm(v~group+age+sex,m)
  stopifnot(max(abs(c(r$ci95_low[1],r$ci95_high[1])-confint(fit)['group',]))<1e-9)
  bad <- read.delim(file.path(root,'fit/context_results.tsv'))
  original <- bad
  for (field in c('beta','stderr','raw_p','wrapper_q','p_for_BH')) {
    bad <- original; bad[[field]][1] <- bad[[field]][1]+.01
    write_tsv(bad,file.path(root,'fit/context_results.tsv'))
    stopifnot(inherits(try(audit_da1_context(root,output),silent=TRUE),'try-error'))
  }
  write_tsv(original,file.path(root,'fit/context_results.tsv'))
  m$age[1] <- m$age[1]+10; write_tsv(m,file.path(root,'metadata.tsv'))
  stopifnot(inherits(try(audit_da1_context(root,output),silent=TRUE),'try-error'))
  cat('[PASS] DA1 audit: exact adjusted CI, both directions, constants/perfect fits, tampering\n')
}
