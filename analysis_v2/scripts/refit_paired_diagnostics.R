#!/usr/bin/env Rscript
source('analysis_v2/scripts/audit_paired_pilot.R')

select_diagnostics <- function(variation,models) {
    selected <- list()
    for (context in sort(unique(variation$context_id))) {
        v <- variation[variation$context_id==context,,drop=FALSE]
        m <- models[models$context_id==context & models$mapping_status=='EXACT_NAME',,drop=FALSE]
        if(anyDuplicated(m$model_feature)) stop('Duplicate saved-model mapping')
        v <- v[order(v$model_feature),,drop=FALSE]
        j <- match(v$model_feature,m$model_feature)
        available <- !is.na(j)
        code <- m$optimizer_code[j]; message <- m$convergence_messages[j]
        no_error <- available & !is.na(m$diagnostic_error[j]) & m$diagnostic_error[j]==''
        categories <- list(optimizer_failure=which(no_error & !is.na(code) & code!='' & code!='0'),
            convergence_warning=which(no_error & code=='0' & !is.na(message) & nzchar(message)),
            clean_comparison=which(no_error & code=='0' & message=='' & !v$backend_warning_prefix & !m$singular[j]),
            deterministic_difference=which(no_error & v$deterministic_difference_at_tolerance),
            Fnuc=which(available & v$feature=='Fusobacterium nucleatum'))
        chosen <- list()
        for (category in names(categories)) {
            candidates <- categories[[category]]
            if (!length(candidates)) next
            i <- candidates[1]; key <- v$model_feature[i]
            if(is.null(chosen[[key]])) chosen[[key]] <- list(i=i,reasons=category)
            else chosen[[key]]$reasons <- c(chosen[[key]]$reasons,category)
        }
        for (key in names(chosen)) {
            row <- v[chosen[[key]]$i,,drop=FALSE]
            row$selection_reasons <- paste(chosen[[key]]$reasons,collapse=';')
            selected[[length(selected)+1]] <- row
        }
    }
    if(!length(selected)) stop('No selected features')
    do.call(rbind,selected)
}

refit_diagnostics <- function(audit,results,out) {
    if(!requireNamespace('lmerTest',quietly=TRUE)) stop('lmerTest required')
    read <- function(path) read.delim(path,quote='',check.names=FALSE,stringsAsFactors=FALSE)
    selection <- select_diagnostics(read(file.path(audit,'paired_variation.tsv')),
                                    read(file.path(audit,'model_diagnostics.tsv')))
    write_audit_tsv(selection,file.path(out,'selection.tsv'))
    collected <- list(); paired <- list()
    dir.create(file.path(out,'models'))
    for (attempt in unique(selection$attempt)) {
        task <- file.path(results,'attempts',attempt)
        meta <- read.delim(file.path(task,'metadata.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
        abundance <- as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE))
        meta <- meta[rownames(abundance),,drop=FALSE]
        models <- readRDS(file.path(task,'fit/maaslin_native/fits/models.rds'))
        rows <- selection[selection$attempt==attempt,,drop=FALSE]
        for (i in seq_len(nrow(rows))) {
            row <- rows[i,,drop=FALSE]; original <- models[[row$model_feature]]
            if(!inherits(original,'merMod')) stop('Selected model is not merMod')
            y <- log2(1+abundance[,row$feature]/1e-8)
            data <- data.frame(y=y,group=as.numeric(meta$group),biological_sample_id=factor(meta$biological_sample_id))
            # Ensure original model response and pair/group identities match.
            frame <- model.frame(original)
            key <- paste(data$biological_sample_id,data$group,sep='|')
            original_key <- paste(frame$biological_sample_id,frame$group,sep='|')
            order <- match(key,original_key)
            if(anyDuplicated(key) || anyNA(order) || length(order)!=nrow(frame) ||
               max(abs(model.response(frame)[order]-y))>1e-10) stop('Saved model differs from reconstructed input')
            if(!setequal(names(lme4::fixef(original)),c('(Intercept)','group'))) stop('Unexpected fixed effects')
            warning_free <- function(expr) {
                warnings <- messages <- character(); error <- ''
                value <- tryCatch(withCallingHandlers(expr,
                    warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
                    message=function(m) {messages <<- c(messages,conditionMessage(m));invokeRestart('muffleMessage')}),
                    error=function(e) {error <<- conditionMessage(e);NULL})
                list(value=value,warnings=warnings,messages=messages,error=error)
            }
            for (variant in c('saved_original','raw_nloptwrap','scaled_nloptwrap','scaled_bobyqa','scaled_Nelder_Mead')) {
                scaled <- startsWith(variant,'scaled_'); scale <- if(scaled) sd(y) else 1
                mean <- if(scaled) mean(y) else 0
                if(!is.finite(scale) || scale<=0) scale <- 1
                data$y <- (y-mean)/scale
                optimizer <- sub('^(raw|scaled)_','',variant)
                captured <- if(variant=='saved_original') list(value=original,warnings=character(),messages=character(),error='')
                    else warning_free(lmerTest::lmer(y~group+(1|biological_sample_id),data=data,
                        REML=lme4::isREML(original),control=lme4::lmerControl(optimizer=optimizer,
                            optCtrl=if(optimizer=='nloptwrap') list(maxeval=100000) else list(maxfun=100000))))
                model <- captured$value
                stats <- list(beta=NA_real_,stderr=NA_real_,df=NA_real_,raw_p=NA_real_,singular=NA,
                    optimizer_code='',convergence_messages='',residual_sd=NA_real_,random_intercept_sd=NA_real_)
                if(!is.null(model)) {
                    extracted <- warning_free({
                        co <- coef(summary(model))['group',]
                        vc <- as.data.frame(lme4::VarCorr(model))
                        list(beta=unname(co['Estimate'])*scale,stderr=unname(co['Std. Error'])*scale,
                            df=unname(co['df']),raw_p=unname(co['Pr(>|t|)']),singular=lme4::isSingular(model),
                            optimizer_code=paste(model@optinfo$conv$opt,collapse=';'),
                            convergence_messages=paste(model@optinfo$conv$lme4$messages,collapse=';'),
                            residual_sd=sigma(model)*scale,
                            random_intercept_sd=vc$sdcor[vc$grp=='biological_sample_id' & is.na(vc$var2)][1]*scale)
                    })
                    if(!is.null(extracted$value)) stats <- extracted$value
                    captured$warnings <- c(captured$warnings,extracted$warnings)
                    captured$error <- paste(c(captured$error,extracted$error)[nzchar(c(captured$error,extracted$error))],collapse=';')
                    saveRDS(model,file.path(out,'models',paste(row$context_id,row$model_feature,variant,sep='_')))
                }
                collected[[length(collected)+1]] <- data.frame(context_id=row$context_id,attempt=attempt,
                    feature=row$feature,model_feature=row$model_feature,selection_reasons=row$selection_reasons,
                    variant=variant,response_scale=scale,response_center=mean,REML=lme4::isREML(original),stats,
                    captured_warnings=paste(captured$warnings,collapse=';'),captured_messages=paste(captured$messages,collapse=';'),
                    error=captured$error,stringsAsFactors=FALSE)
            }
            ids <- sort(unique(meta$biological_sample_id))
            originals <- y[meta$group==0][match(ids,meta$biological_sample_id[meta$group==0])]
            spikes <- y[meta$group==1][match(ids,meta$biological_sample_id[meta$group==1])]
            delta <- spikes-originals
            test <- warning_free(t.test(delta,mu=0))
            paired[[length(paired)+1]] <- data.frame(context_id=row$context_id,model_feature=row$model_feature,
                feature=row$feature,n_pairs=length(delta),mean_difference=mean(delta),sd_difference=sd(delta),
                standard_error=sd(delta)/sqrt(length(delta)),raw_p=if(is.null(test$value)) NA_real_ else test$value$p.value,
                df=length(delta)-1,error=test$error,warnings=paste(test$warnings,collapse=';'),
                label='DIAGNOSTIC_PAIRED_T_NOT_PRIMARY_NO_Q',stringsAsFactors=FALSE)
        }
    }
    write_audit_tsv(do.call(rbind,collected),file.path(out,'refit_comparisons.tsv'))
    write_audit_tsv(do.call(rbind,paired),file.path(out,'paired_difference_checks.tsv'))
    writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
    cat('[PASS] Diagnostic refits:',nrow(selection),'features;',length(collected),'model records; no primary results replaced\n')
}

if(sys.nframe()==0) {
    args <- commandArgs(trailingOnly=TRUE)
    if(length(args)!=3) stop('Usage: refit_paired_diagnostics.R audit-root pilot-results fresh-output')
    refit_diagnostics(args[1],args[2],args[3])
}
