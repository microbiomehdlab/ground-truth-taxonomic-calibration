source('analysis_v2/scripts/refit_paired_diagnostics.R')
blank_path <- tempfile()
write_audit_tsv(data.frame(optimizer_code='0',convergence_messages='',diagnostic_error='',singular=FALSE),blank_path)
blank <- read_diagnostic_table(blank_path)
stopifnot(identical(blank$diagnostic_error,''),identical(blank$convergence_messages,''),!blank$singular)
v <- data.frame(context_id='pilot_test',attempt='attempt_test',model_feature=sprintf('feature_%06d',1:4),
                feature=c('bad','warning','clean','Fusobacterium nucleatum'),
                backend_warning_prefix=c(TRUE,TRUE,FALSE,TRUE),deterministic_difference_at_tolerance=c(TRUE,FALSE,FALSE,FALSE))
m <- data.frame(context_id='pilot_test',model_feature=v$model_feature,mapping_status='EXACT_NAME',
                diagnostic_error='',optimizer_code=c('-4','0','0','0'),
                convergence_messages=c('failure','Hessian','','gradient'),singular=FALSE)
s <- select_diagnostics(v,m)
stopifnot(nrow(s)==4,'Fnuc' %in% s$selection_reasons,
          identical(select_diagnostics(v[4:1,],m[4:1,]),s))
if(requireNamespace('lmerTest',quietly=TRUE)) {
    set.seed(93)
    n <- 12; id <- rep(sprintf('person_%02d',seq_len(n)),each=2); group <- rep(0:1,n)
    y <- 5+rep(rnorm(n,0,1),each=2)+group*.5+rnorm(2*n,0,.15)
    data <- data.frame(y=y,group=group,biological_sample_id=factor(id))
    model <- lmerTest::lmer(y~group+(1|biological_sample_id),data=data)
    root <- tempfile(); audit <- file.path(root,'audit'); results <- file.path(root,'results')
    out <- file.path(root,'output'); task <- file.path(results,'attempts','attempt_test')
    dir.create(audit,recursive=TRUE); dir.create(out); dir.create(file.path(task,'fit/maaslin_native/fits'),recursive=TRUE)
    v <- v[4,,drop=FALSE]; m <- m[4,,drop=FALSE]
    v$backend_warning_prefix <- FALSE; m$convergence_messages <- ''
    write_audit_tsv(v,file.path(audit,'paired_variation.tsv'))
    write_audit_tsv(m,file.path(audit,'model_diagnostics.tsv'))
    meta <- data.frame(observation_id=paste(id,group,sep='_'),biological_sample_id=id,group=group)
    write.table(meta,file.path(task,'metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
    abundance <- data.frame(observation_id=meta$observation_id,'Fusobacterium nucleatum'=(2^y-1)*1e-8,check.names=FALSE)
    write.table(abundance,file.path(task,'abundance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
    saveRDS(setNames(list(model),v$model_feature),file.path(task,'fit/maaslin_native/fits/models.rds'))
    refit_diagnostics(audit,results,out)
    compared <- read.delim(file.path(out,'refit_comparisons.tsv'),quote='',colClasses=c(error='character'))
    if(any(nzchar(compared$error))) print(compared[,c('variant','error')])
    stopifnot(nrow(compared)==5,all(compared$error==''),
              max(abs(compared$beta-lme4::fixef(model)['group']))<1e-5,
              max(abs(compared$stderr-compared$stderr[1]))<1e-4)
    cat('[PASS] Actual mixed-model refit and original-unit rescaling fixture\n')
} else {
    if(Sys.getenv('REQUIRE_REFIT_BACKEND')=='1') stop('Required lmerTest backend unavailable')
    cat('[NOTE] Actual mixed-model backend unavailable locally; container test required\n')
}
cat('[PASS] Deterministic diagnostic selection; no outcome-based selection\n')
