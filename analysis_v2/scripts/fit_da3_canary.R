args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: assembled-task repository')
task<-args[1];repo<-args[2]
source(file.path(repo,'analysis_v2/lib/exact_da3_candidate.R'))
source(file.path(repo,'analysis_v2/lib/monte_carlo_da3_candidate.R'))
source(file.path(repo,'analysis_v2/lib/unpaired_null.R'))
a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE));storage.mode(a)<-'double'
m<-read.delim(file.path(task,'metadata.tsv'),stringsAsFactors=FALSE)
c<-read.delim(file.path(task,'context.tsv'),stringsAsFactors=FALSE)
if(anyDuplicated(m$biological_sample_id)||!identical(m$observation_id,rownames(a)))stop('Invalid biological identities')
y<-log2(1+a/1e-8);g<-as.integer(m$group)
start<-proc.time()[['elapsed']]
fit<-if(c$n==5)exact_da3_candidate(y,g) else monte_carlo_da3_candidate(y,g,repetitions=9999L,seed=20261002L)
elapsed<-proc.time()[['elapsed']]-start
old<-unpaired_null_results(y,g)
write.table(fit,file.path(task,'candidate_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(old,file.path(task,'parametric_comparison.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(data.frame(context_id=c$context_id,n=c$n,arm=c$arm,inference_seconds=elapsed,
    permutations=if(c$n==5)252 else 9999,method=if(c$n==5)'EXACT' else 'MONTE_CARLO',
    positive_fwer_discoveries=sum(fit$beta>0 & fit$max_statistic_fwer_p<=.05),
    interpretation='ENGINEERING_TIMING_NOT_PRODUCTION_AUTHORIZATION'),
    file.path(task,'fit_timing.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
