args<-commandArgs(trailingOnly=TRUE);if(length(args)!=2)stop('Usage: task repo')
task<-args[1];repo<-args[2]
source(file.path(repo,'analysis_v2/lib/exact_da3_candidate.R'));source(file.path(repo,'analysis_v2/lib/monte_carlo_da3_candidate.R'));source(file.path(repo,'analysis_v2/lib/unpaired_null.R'))
a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE));storage.mode(a)<-'double'
m<-read.delim(file.path(task,'metadata.tsv'),stringsAsFactors=FALSE);c<-read.delim(file.path(task,'context.tsv'),stringsAsFactors=FALSE)
if(!identical(m$observation_id,rownames(a))||anyDuplicated(m$observation_id))stop('Identity mismatch')
y<-log2(1+a/1e-8);start<-proc.time()[['elapsed']]
f<-if(c$n==5)exact_da3_candidate(y,m$group)else monte_carlo_da3_candidate(y,m$group,9999L,as.integer(c$permutation_seed))
old<-unpaired_null_results(y,m$group)
f$parametric_p<-old$raw_p;f$parametric_bh_q<-old$wrapper_q
write.table(f,file.path(task,'results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(data.frame(context_id=c$context_id,cohort=c$cohort,background=c$background,profiler=c$profiler,n=c$n,arm=c$arm,anchor=c$anchor,allocation_id=c$allocation_id,family_n=ncol(a),positive_fwer_discoveries=sum(f$beta>0 & f$max_statistic_fwer_p<=.05),positive_parametric_discoveries=sum(old$estimable & !is.na(old$beta) & old$beta>0 & old$wrapper_q<=.05),inference_seconds=proc.time()[['elapsed']]-start),file.path(task,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
