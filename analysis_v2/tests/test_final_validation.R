source('analysis_v2/lib/maaslin_contract.R')
source('analysis_v2/lib/maaslin_context.R')
source('analysis_v2/lib/unpaired_null.R')
set.seed(123)
a <- matrix(runif(200,1e-5,.01),20,dimnames=list(paste0('person',1:20),paste0('f',1:10)))
a[,9] <- 0; a[,10] <- rep(c(.0001,.001),each=10)
group <- rep(0:1,each=10); y <- log2(1+a/1e-8)
x <- unpaired_null_results(y,group)
for (i in 1:8) {
  m <- summary(lm(y[,i]~group))$coefficients
  stopifnot(isTRUE(all.equal(x$beta[i],unname(m[2,1]),tolerance=1e-10)),
    isTRUE(all.equal(x$raw_p[i],unname(m[2,4]),tolerance=1e-10)))
}
stopifnot(identical(x$status[9:10],c('NON_ESTIMABLE_CONSTANT','NON_ESTIMABLE_PERFECT_FIXED_FIT')))
if (requireNamespace('Maaslin2',quietly=TRUE)) {
  task <- tempfile('validation-backend-')
  actual <- fit_maaslin_context(a,data.frame(group=group,row.names=rownames(a)),colnames(a),task)
  stopifnot(identical(actual$status,x$status),
    isTRUE(all.equal(actual$raw_p,x$raw_p,tolerance=1e-7)),
    isTRUE(all.equal(actual$beta,x$beta,tolerance=1e-7)))
  unlink(task,recursive=TRUE)
  cat('[PASS] Actual pinned MaAsLin2/fast group-only LM equivalence\n')
} else {
  if (Sys.getenv('REQUIRE_VALIDATION_BACKEND')=='1') stop('Actual pinned backend required')
  cat('[SKIP] Pinned backend unavailable locally; enforced by cluster submission\n')
}
source('analysis_v2/lib/paired_difference_context.R')
task <- tempfile('sensitivity-test-'); dir.create(task)
meta <- data.frame(biological_sample_id=rep(paste0('person',1:10),2),
  group=group,spike_state=rep(c('original','spiked'),each=10),row.names=rownames(a))
write.table(data.frame(observation_id=rownames(a),a,check.names=FALSE),file.path(task,'abundance.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(observation_id=rownames(a),meta),file.path(task,'metadata.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
code <- system2(file.path(R.home('bin'),'Rscript'),c('analysis_v2/scripts/run_final_validation.R',shQuote(task),'sensitivity',shQuote(getwd())))
saved <- read.delim(file.path(task,'sensitivity_results.tsv'))
stopifnot(code==0,nrow(saved)==30,length(unique(saved$pseudocount))==3)
primary <- paired_difference_results(a,meta,colnames(a))$results
stopifnot(isTRUE(all.equal(saved$raw_p[saved$pseudocount==1e-8],primary$raw_p)))
unlink(task,recursive=TRUE)
task <- tempfile('synthetic-test-'); dir.create(task)
write.table(data.frame(n_pairs=10,family_n=20,scenario='skewed_zero_mean'),
  file.path(task,'settings.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
code <- system2(file.path(R.home('bin'),'Rscript'),c('analysis_v2/scripts/run_final_validation.R',shQuote(task),'synthetic',shQuote(getwd())))
draws <- read.delim(file.path(task,'null_draws.tsv'))
summary <- read.delim(file.path(task,'summary.tsv'))
stopifnot(code==0,nrow(draws)==1000,summary$repetitions==1000,
  summary$any_discovery_rate==mean(draws$discoveries>0))
unlink(task,recursive=TRUE)
cat('[PASS] Group-only LM formulas, degeneracy and file-based pseudocount sensitivity\n')
