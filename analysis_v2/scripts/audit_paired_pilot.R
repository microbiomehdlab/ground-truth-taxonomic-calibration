#!/usr/bin/env Rscript
# Diagnostic only: no refitting, warning suppression or DA p-value changes.
paired_metrics <- function(a,b,p=1e-8) {
    stopifnot(length(a)==length(b),length(a)>1,all(is.finite(c(a,b))),
              all(c(a,b)>=0),all(c(a,b)<=1))
    x <- log2(1+a/p); y <- log2(1+b/p); delta <- y-x
    tolerance <- 100*.Machine$double.eps*max(1,max(abs(c(x,y))))
    data.frame(pairs=length(a),original_positive=sum(a>0),spiked_positive=sum(b>0),
        exactly_unchanged_pairs=sum(a==b),max_absolute_native_change=max(abs(b-a)),
        original_log_sd=sd(x),spiked_log_sd=sd(y),mean_log_difference=mean(delta),
        sd_log_difference=sd(delta),range_log_difference=diff(range(delta)),
        deterministic_difference_at_tolerance=diff(range(delta))<=tolerance,
        difference_tolerance=tolerance)
}

audit_paired <- function(root,out,repo) {
    attempts <- list.dirs(file.path(root,'attempts'),recursive=FALSE,full.names=TRUE)
    attempts <- attempts[file.exists(file.path(attempts,'FAILED.json'))]
    if (!length(attempts)) stop('No retained failed attempts')
    aliases <- read.csv(file.path(repo,'examples/spike_taxon_aliases.csv'),stringsAsFactors=FALSE)
    features <- summaries <- model_rows <- list()
    for (task in sort(attempts)) {
        read <- function(name,...) read.delim(file.path(task,name),check.names=FALSE,
                                             stringsAsFactors=FALSE,...)
        context <- read('context.tsv')
        if (nrow(context)!=1 || context$analysis!='DA2' || context$paired!=1)
            stop('Unexpected non-paired failed context: ',task)
        meta <- read('metadata.tsv',row.names=1)
        abundance <- as.matrix(read('abundance.tsv',row.names=1))
        storage.mode(abundance) <- 'double'
        if (!setequal(rownames(meta),rownames(abundance)) || anyDuplicated(rownames(meta)) ||
            anyDuplicated(rownames(abundance)) || anyDuplicated(colnames(abundance))) stop('Invalid identities')
        meta <- meta[rownames(abundance),,drop=FALSE]
        pairs <- table(meta$biological_sample_id,factor(meta$spike_state,levels=c('original','spiked')))
        if (any(pairs!=1) || anyNA(meta$spike_state) ||
            !all(meta$spike_state %in% c('original','spiked')) ||
            any(meta$group!=as.integer(meta$spike_state=='spiked'))) stop('Invalid pairing')
        ids <- sort(unique(meta$biological_sample_id))
        original <- match(ids,meta$biological_sample_id[meta$spike_state=='original'])
        spiked <- match(ids,meta$biological_sample_id[meta$spike_state=='spiked'])
        a <- abundance[meta$spike_state=='original',,drop=FALSE][original,,drop=FALSE]
        b <- abundance[meta$spike_state=='spiked',,drop=FALSE][spiked,,drop=FALSE]
        map <- read('fit/feature_map.tsv')
        if (!setequal(map$feature,colnames(abundance)) || anyDuplicated(map$feature) ||
            anyDuplicated(map$model_feature)) stop('Invalid feature mapping')
        err_path <- file.path(task,'model.err')
        err <- if(file.exists(err_path)) readLines(err_path,warn=FALSE) else character()
        # Prefixes identify features whose warnings were printed by the backend.
        # warnings.txt itself does not identify taxa: do not assign by position.
        prefixed <- grep('^Feature feature_[0-9]+',err,value=TRUE)
        warned <- unique(sub('^Feature (feature_[0-9]+).*','\\1',prefixed))
        models_path <- file.path(task,'fit/maaslin_native/fits/models.rds')
        models <- if(file.exists(models_path)) readRDS(models_path) else list()
        if (!is.list(models)) stop('Unexpected saved-model layout')
        model_names <- names(models)
        for (i in seq_along(models)) {
            model <- models[[i]]
            name <- if(length(model_names)>=i) model_names[i] else ''
            known <- !is.na(name) && name %in% map$model_feature
            diagnostic <- list(singular=NA,optimizer_code='',convergence_messages='',
                               residual_sd=NA_real_,random_intercept_sd=NA_real_,diagnostic_error='')
            if (inherits(model,'merMod')) {
                diagnostic <- tryCatch({
                    vc <- as.data.frame(lme4::VarCorr(model))
                    list(singular=lme4::isSingular(model),
                         optimizer_code=paste(model@optinfo$conv$opt,collapse=';'),
                         convergence_messages=paste(model@optinfo$conv$lme4$messages,collapse=';'),
                         residual_sd=sigma(model),
                         random_intercept_sd=vc$sdcor[vc$grp=='biological_sample_id' & vc$var1=='(Intercept)' & is.na(vc$var2)][1],
                         diagnostic_error='')
                },error=function(e) { diagnostic$diagnostic_error<-conditionMessage(e); diagnostic })
            } else diagnostic$diagnostic_error <- 'NOT_MERMOD'
            model_rows[[length(model_rows)+1]] <- data.frame(context_id=context$context_id,
                attempt=basename(task),model_index=i,model_name=if(is.na(name)) '' else name,
                model_feature=if(known) name else '',mapping_status=if(known) 'EXACT_NAME' else 'UNMAPPED_NOT_POSITIONALLY_INFERRED',
                model_class=paste(class(model),collapse=';'),diagnostic,stringsAsFactors=FALSE)
        }
        target_names <- aliases$alias[aliases$tool==context$profiler]
        local <- list()
        for (feature in colnames(abundance)) {
            model_feature <- map$model_feature[match(feature,map$feature)]
            local[[length(local)+1]] <- data.frame(context_id=context$context_id,attempt=basename(task),
                cohort=context$cohort,profiler=context$profiler,population=context$population,
                feature=feature,model_feature=model_feature,target=feature %in% target_names,
                backend_warning_prefix=model_feature %in% warned,paired_metrics(a[,feature],b[,feature]),
                stringsAsFactors=FALSE)
        }
        local <- do.call(rbind,local)
        features[[length(features)+1]] <- local
        summaries[[length(summaries)+1]] <- data.frame(context_id=context$context_id,attempt=basename(task),
            cohort=context$cohort,profiler=context$profiler,population=context$population,
            pairs=length(ids),family_n=nrow(local),target_n=sum(local$target),
            warning_prefixed_features=sum(local$backend_warning_prefix),
            warned_targets=sum(local$target & local$backend_warning_prefix),
            deterministic_difference_features=sum(local$deterministic_difference_at_tolerance),
            saved_models=length(models),named_models_mapped=sum(model_names %in% map$model_feature),
            stringsAsFactors=FALSE)
    }
    save <- function(rows,name) {
        if(length(rows)) write.table(do.call(rbind,rows),file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
    }
    save(features,'paired_variation.tsv'); save(summaries,'context_summary.tsv')
    save(model_rows,'model_diagnostics.tsv')
    writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
    print(do.call(rbind,summaries),row.names=FALSE)
}

if (sys.nframe()==0) {
    args <- commandArgs(trailingOnly=TRUE)
    if(length(args)!=3) stop('Usage: audit_paired_pilot.R results fresh-audit-directory repository')
    audit_paired(args[1],args[2],args[3])
}
