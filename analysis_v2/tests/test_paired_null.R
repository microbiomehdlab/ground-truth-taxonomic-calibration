source('analysis_v2/lib/paired_null.R')
set.seed(42)
d <- matrix(rnorm(60),nrow=10)
tolerance <- rep(1e-12,ncol(d))
x <- paired_null_statistics(d,tolerance)
p <- apply(d,2,function(v)t.test(v)$p.value)
stopifnot(x['discoveries']==sum(p.adjust(p,'BH')<=0.05),x['nominal_rejections']==sum(p<=0.05))
a <- paired_signflip_diagnostic(d,tolerance,20)
b <- paired_signflip_diagnostic(d,tolerance,20)
stopifnot(identical(a,b),all(a$signs %in% c(-1,1)),nrow(a$draws)==20)
zero <- paired_null_statistics(matrix(0,10,3),rep(1e-12,3))
constant <- paired_null_statistics(matrix(1,10,3),rep(1e-12,3))
stopifnot(all(zero==0),all(constant==0))
task <- tempfile('null-test-'); dir.create(task); dir.create(file.path(task,'fit'))
colnames(d) <- paste0('feature',1:6); rownames(d) <- paste0('person',1:10)
write.table(data.frame(analysis='DA2',inference_method='paired_differences',context_id='fixture'),
  file.path(task,'context.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(feature=colnames(d),numerical_tolerance=tolerance,n_pairs=10),
  file.path(task,'fit/context_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(biological_sample_id=rownames(d),d,check.names=FALSE),
  file.path(task,'fit/paired_differences.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
out <- file.path(task,'output')
status <- system2(file.path(R.home('bin'),'Rscript'),c('analysis_v2/scripts/run_da2_signflip.R',
  shQuote(task),shQuote(out),shQuote(getwd())))
summary <- read.delim(file.path(out,'summary.tsv'))
stopifnot(status==0,summary$repetitions==1000,summary$unique_sign_vectors<=1024,
  summary$unique_sign_vectors>1,summary$family_n==6)
unlink(task,recursive=TRUE)
cat('[PASS] Sign-flip reproducibility, t-test/BH equivalence and degenerate differences\n')
