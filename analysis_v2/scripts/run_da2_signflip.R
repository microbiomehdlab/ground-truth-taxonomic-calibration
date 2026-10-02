#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop('Usage: source-context output repository')
source_dir <- args[1]; out <- args[2]; repo <- args[3]
source(file.path(repo,'analysis_v2/lib/paired_null.R'))
if (file.exists(out)) stop('Fresh output required')
context <- read.delim(file.path(source_dir,'context.tsv'),stringsAsFactors=FALSE)
if (nrow(context)!=1 || context$analysis!='DA2' || context$inference_method!='paired_differences')
  stop('DA2 paired-difference source required')
results <- read.delim(file.path(source_dir,'fit/context_results.tsv'),check.names=FALSE)
data <- read.delim(file.path(source_dir,'fit/paired_differences.tsv'),check.names=FALSE,row.names=1)
d <- as.matrix(data); storage.mode(d) <- 'double'
if (!identical(colnames(d),results$feature) || anyDuplicated(rownames(d)) ||
    any(results$n_pairs!=nrow(d))) stop('Difference identity mismatch')
# All contexts use this seed; matched populations with identical sorted person
# identities therefore receive identical signs across profilers.
x <- paired_signflip_diagnostic(d,results$numerical_tolerance)
dir.create(out,recursive=TRUE)
write.table(x$draws,file.path(out,'null_draws.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(biological_sample_id=rownames(d),x$signs,check.names=FALSE),
  file.path(out,'signs.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
count <- sum(x$draws$discoveries>0); interval <- binom.test(count,nrow(x$draws))$conf.int
write.table(data.frame(context_id=context$context_id,n_pairs=nrow(d),family_n=ncol(d),
  repetitions=nrow(x$draws),any_discovery_rate=count/nrow(x$draws),
  mc_ci_low=interval[1],mc_ci_high=interval[2],
  mean_discoveries=mean(x$draws$discoveries),
  unique_sign_vectors=nrow(unique(t(x$signs))),
  diagnostic='CONDITIONAL_SIGNFLIP_SYMMETRY_ASSUMED_NOT_BIOLOGICAL_NULL'),
  file.path(out,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('[PASS] Conditional DA2 sign-flip diagnostic:',context$context_id,'\n')
