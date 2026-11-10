#!/usr/bin/env Rscript
# Analytic/synthetic checks only; no biological response file is read.
ff<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',ff)));source(file.path(here,'statistics_core.R'));set.seed(2026110611)
x<-seq_len(20);y<-x*2+3;ip<-sp_perm_index(20,999L);m<-sp_perm_metrics(x,y,ip);stopifnot(abs(sp_rho(x,y)-1)<1e-12,all(is.finite(m[,'rho'])),max(abs(m[,'direction']-1))<1e-12,sp_p(1,m[,'direction'])==1)
ag<-sp_agreement(rep(1,20),c(rep(1,18),-1,-1));stopifnot(ag['concordance']==.9,is.na(ag['kappa'])||ag['kappa']==0,is.na(ag['MCC']))
x<-c(rep(-1,10),rep(1,10));ag<-sp_agreement(x,x);stopifnot(ag['concordance']==1,ag['kappa']==1,ag['MCC']==1)
large<-rep(c(-1,1),219);stopifnot(sp_agreement(large,large)['MCC']==1)
ip<-sp_perm_index(20,999L);d<-sp_perm_metrics(x,x,ip);stopifnot(sp_p(1,d[,'direction'])<=.01)
cc<-cbind(a=seq_len(20),b=rep(1:4,5));f<-sp_partial_fit(cc[,1],cc[,1],cc);stopifnot(!is.finite(f$rho))
z<-sp_partial(rnorm(20),rnorm(20),cc,B=150L,P=199L);stopifnot(z$summary$genes==20,z$summary$residual_df>=3,is.finite(z$summary$raw_permutation_p))
empty<-sp_pair(data.table(gene=character(),source_deltaZ=numeric(),validation_deltaZ=numeric()),B=100L,P=199L);stopifnot(all(is.na(empty$results$estimate)),all(is.na(empty$results$raw_permutation_p)))
stopifnot(all(is.na(sp_adjust(c(NA_real_,NA_real_),2))))
cat('PASS: analytic rho, both-marginal sign null, kappa/MCC degeneracy, conditional partial residuals, fixed-family unavailable tests. No biological files loaded.\n')
