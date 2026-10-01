source('analysis_v2/lib/maaslin_contract.R')
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),'try-error'))
x <- matrix(c(0,.01,.001,.02,0,0,0,0),nrow=4,
            dimnames=list(paste0('obs',1:4),c('target','zero')))
metadata <- data.frame(biological_sample_id=c('A','A','B','B'),
                       spike_state=c('original','spiked','original','spiked'),
                       row.names=rownames(x))
p <- prepare_maaslin_input(x,metadata[4:1,],colnames(x),paired=TRUE)
stopifnot(identical(rownames(p$metadata),rownames(x)),all(p$data$zero==0))
expected <- log2(x+1e-8)-log2(1e-8)
stopifnot(max(abs(as.matrix(p$data)-expected))<1e-12)
y <- log2(x[,'target']+1e-8)
group <- c(0,1,0,1)
stopifnot(abs(coef(lm(y~group))[2]-coef(lm(p$data$target~group))[2])<1e-12)
fails(prepare_maaslin_input(x,metadata[-1,],colnames(x)))
bad <- metadata; bad$spike_state[1]<-'spiked'
fails(prepare_maaslin_input(x,bad,colnames(x),paired=TRUE))
badx <- x; badx[1,1]<-NA
fails(prepare_maaslin_input(badx,metadata,colnames(x)))
badx <- x; badx[1,1]<-10
fails(prepare_maaslin_input(badx,metadata,colnames(x)))
fails(prepare_maaslin_input(x,metadata,c('missing')))
fails(prepare_maaslin_input(x,metadata,c('target','target')))
q <- full_family_bh(c('A','B','C'),c('C','A','B'),c(NA,.01,.04),c(FALSE,TRUE,TRUE))
stopifnot(max(abs(q$wrapper_q-c(.03,.06,1)))<1e-12)
stopifnot(is.na(q$raw_p[3]),q$p_for_BH[3]==1)
fails(full_family_bh(c('A','B'),c('A'),.01,TRUE))
fails(full_family_bh('A','A',1e-300,FALSE))
fails(full_family_bh('A','A',NA,TRUE))
cat('[PASS] MaAsLin input, fixed transform, pairing and full-family BH contracts\n')
