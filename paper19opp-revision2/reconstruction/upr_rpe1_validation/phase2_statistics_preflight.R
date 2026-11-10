#!/usr/bin/env Rscript
# Algebraic checks only; these values are not biological observations.
a<-commandArgs(TRUE);here<-if(length(a))a[1]else'paper19opp-revision2/reconstruction/upr_rpe1_validation';source(file.path(here,'phase2_statistics.R'))
g<-sprintf('algebra_%03d',1:120);x<-setNames(120:1,g);y<-x/2
src<-data.table(gene=g,source_delta_z=x,priority_score=x,stratum=rep(1:6,20));val<-data.table(gene=g,validation_delta_z=y);ranking<-data.table(gene=g,rank=1:120)
matches<-rbindlist(lapply(1:20,function(i)data.table(selected_gene=g[i],background_gene=g[(21+5*(i-1)):(25+5*(i-1))])))
gem<-rbindlist(lapply(1:3,function(k)data.table(gene=g,gem=as.character(k),target_cells=40,contrast_z=x+c(-.1,0,.1)[k])));vgem<-copy(gem);vgem[,contrast_z:=contrast_z/2]
z<-phase2_primary(src,val,ranking,matches,gem,vgem,B=200L,P=500L,seed=17L)
stopifnot(abs(z$endpoints$estimate[1]-1)<1e-12,abs(z$endpoints$estimate[2]-30)<1e-12,abs(z$endpoints$estimate[3]-20)<1e-12)
stopifnot(max(abs(z$bootstrap$wild_Top20-30))<1e-10,max(abs(z$bootstrap$wild_Top10-20))<1e-10,max(abs(z$bootstrap$wild_rho-1))<1e-12)
roles<-matrix(1:6,ncol=1);null<-phase2_set_permutation(y,z$sets,roles,k=1L);stopifnot(abs(mean(null))<1e-12)
for(i in 1:20){ii<-phase2_stratified_permutation(src$stratum);stopifnot(all(src$stratum[ii]==src$stratum))}
cat('Analytic estimates, shared-GEM invariance, six-role centering and strata checks passed.\n')
