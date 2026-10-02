# Engineering candidate: null-restricted HC2 residual bootstrap, HC3 studentization.
# Source clinical_hc3.R first. Not an automatically approved production method.
clinical_wild_bootstrap <- function(x, y, repetitions=9999L, seed=20261003L,
                                    multipliers=NULL) {
 if(length(repetitions)!=1L || !is.finite(repetitions) || repetitions<1 ||
    repetitions!=as.integer(repetitions)) stop('Invalid repetitions')
 if(is.null(rownames(x)) || is.null(rownames(y)) ||
    anyDuplicated(rownames(x)) || !identical(rownames(x),rownames(y)))
   stop('Clinical sample identities must match in order')
 if(is.null(colnames(y)) || anyDuplicated(colnames(y))) stop('Invalid feature identities')
 observed<-clinical_hc3(x,y)
 j<-match('group',colnames(x)); x0<-x[,-j,drop=FALSE]
 if(!ncol(x0) || qr(x0)$rank!=ncol(x0)) stop('Invalid restricted design')
 b0<-solve(crossprod(x0)); h0<-rowSums((x0%*%b0)*x0)
 if(any(1-h0<1e-10)) stop('Extreme restricted leverage')
 f0<-lm.fit(x0,y); r0<-matrix(f0$residuals,nrow=nrow(x))
 mu<-y-r0
 residual<-r0/sqrt(1-h0)
 if(is.null(multipliers)) {
   # Restore caller RNG; identical seed gives identical person weights across feature shards.
   had_rng<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
   if(had_rng) old_rng<-get('.Random.seed',envir=.GlobalEnv)
   on.exit(if(had_rng) assign('.Random.seed',old_rng,envir=.GlobalEnv)
           else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
             rm('.Random.seed',envir=.GlobalEnv),add=TRUE)
   set.seed(seed)
   multipliers<-matrix(sample(c(-1,1),nrow(x)*repetitions,replace=TRUE),nrow(x))
 }
 if(!is.matrix(multipliers) || !identical(dim(multipliers),c(nrow(x),as.integer(repetitions))) ||
    any(!is.finite(multipliers)) || any(!multipliers %in% c(-1,1)))
   stop('Expected person-by-draw Rademacher multipliers')
 bread<-solve(crossprod(x)); w<-x%*%bread
 h<-rowSums(w*x); threshold<-abs(observed$beta/observed$hc3_se)
 exceed<-integer(ncol(y))
 for(k in seq_len(repetitions)) {
   z<-mu+residual*multipliers[,k]
   fit<-lm.fit(x,z)
   se<-sqrt(colSums((matrix(fit$residuals,nrow=nrow(x))*(w[,j]/(1-h)))^2))
   tstar<-abs(matrix(fit$coefficients,nrow=ncol(x))[j,]/se)
   # Undefined draws count as exceedances: never silently reduce the denominator.
   exceed<-exceed+as.integer(!is.finite(tstar) | (tstar>=threshold))
 }
 p<-ifelse(observed$estimable,(exceed+1)/(repetitions+1),NA_real_)
 observed$wild_p<-p
 observed$wild_bh_q<-p.adjust(ifelse(observed$estimable,p,1),'BH')
 observed$bootstrap_draws<-repetitions
 observed$minimum_mc_p<-1/(repetitions+1)
 observed$bootstrap_reference<-'NULL_RESTRICTED_HC2_RADEMACHER_HC3_STUDENTIZED_CANDIDATE'
 observed
}
