source('analysis_v2/scripts/audit_paired_pilot.R')
escaped_path <- tempfile()
messages <- c('very large eigenvalue\n - Rescale variables?',
              'tab\there\rreturn and \\literal and "quote"', '')
write_audit_tsv(data.frame(context_id=rep('pilot_test',3),convergence_messages=messages),escaped_path)
escaped <- read.delim(escaped_path,quote='',check.names=FALSE,colClasses='character')
stopifnot(nrow(escaped)==3,all(escaped$context_id=='pilot_test'),
          identical(escaped$convergence_messages,encodeString(messages,quote='')),
          length(readLines(escaped_path))==4,
          all(count.fields(escaped_path,sep='\t',quote='')==2))
x <- paired_metrics(c(0,.01,.02),c(0,.01,.02))
stopifnot(x$exactly_unchanged_pairs==3,x$sd_log_difference==0,x$deterministic_difference_at_tolerance)
x <- paired_metrics(c(.01,.02,.03),c(.02,.025,.04))
stopifnot(x$exactly_unchanged_pairs==0,x$sd_log_difference>0,!x$deterministic_difference_at_tolerance)
stopifnot(inherits(try(paired_metrics(c(NA,.1),c(.1,.2)),silent=TRUE),'try-error'))
root <- tempfile(); task <- file.path(root,'attempts','pilot_test'); out <- tempfile()
dir.create(file.path(task,'fit'),recursive=TRUE); dir.create(out)
write <- function(x,name) write.table(x,file.path(task,name),sep='\t',quote=FALSE,row.names=FALSE)
writeLines('{}',file.path(task,'FAILED.json'))
write(data.frame(context_id='pilot_test',analysis='DA2',paired=1,cohort='feng',
                 profiler='metaphlan4',population='community'),'context.tsv')
write(data.frame(observation_id=c('a0','b1','a1','b0'),biological_sample_id=c('a','b','a','b'),
                 group=c(0,1,1,0),spike_state=c('original','spiked','spiked','original')),'metadata.tsv')
write(data.frame(observation_id=c('b0','a1','b1','a0'),
                 'Fusobacterium nucleatum'=c(.02,.02,.04,.01),check.names=FALSE),'abundance.tsv')
write(data.frame(model_feature='feature_000001',feature='Fusobacterium nucleatum'),'fit/feature_map.tsv')
writeLines('Feature feature_000001 : simpleWarning test',file.path(task,'model.err'))
audit_paired(root,out,getwd())
result <- read.delim(file.path(out,'paired_variation.tsv'))
stopifnot(result$pairs==2,result$target,result$backend_warning_prefix,
          abs(result$max_absolute_native_change-.02)<1e-12)
cat('[PASS] Paired-variation metrics; no model refits\n')
