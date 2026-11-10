#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(if(length(args))args[1]else'paper19opp-revision2');p<-file.path(project,'upr_robustness','ranking');out<-file.path(p,'alternative_definitions');src<-file.path(project,'reconstruction','cross_dataset','phase2_upr_rpe1');modules<-file.path(project,'upr_robustness','modules');.libPaths(c(file.path(dirname(project),'.Rlib'),.libPaths()));suppressPackageStartupMessages({library(data.table);library(digest)})
set.seed(2026110608);B<-2000L
man<-fread(file.path(out,'source_lock_manifest.tsv'));stopifnot(digest(file.path(out,'source_lock_manifest.tsv'),file=TRUE,algo='sha256')=='b94e9a421c1e6e1a7ceee9829bbe894a7c88bd2ee685438044da27f2680007e2',all(vapply(file.path(project,man$file),digest,'',file=TRUE,algo='sha256')==man$sha256))
d<-fread(file.path(modules,'scored','RPE1_gene_effects.tsv.gz'));gg<-fread(file.path(modules,'scored','RPE1_GEM_effects.tsv.gz'));point<-fread(file.path(out,'alternative_fixed_background_gaps.tsv'))
meta<-fread(file.path(src,'source_rank_inputs.tsv'));ids<-sort(meta$gene);cc<-sort(fread(file.path(src,'complete_covariate_gene_roster.tsv'))$gene);cohorts<-list(all438=ids,complete370=cc)
old<-list(all438=fread(file.path(src,'frozen_background_matches.tsv')),complete370=fread(file.path(src,'complete_covariate_background_matches.tsv')));new<-fread(file.path(out,'source_alternative_background_matches.tsv'))
state<-readRDS(file.path(project,'reconstruction','upr_rpe1_validation','phase2_statistical_state.rds'));signs<-state$primary$validation_wild$wild_signs;gems<-state$primary$validation_wild$gem_labels
ci<-function(v){v<-v[is.finite(v)];if(!length(v))return(c(NA_real_,NA_real_));as.numeric(quantile(v,c(.025,.975)))}
ans<-list();j<-0L;draws<-list()
for(a in unique(point$axis))for(mm in unique(point$mode)){
 dd<-d[axis==a&mode==mm&gene%in%ids];yall<-setNames(dd$delta_z,dd$gene)[ids];gr<-gg[axis==a&mode==mm&gene%in%ids];tot<-gr[,.(n=sum(n_target)),by=gene];total<-setNames(tot$n,tot$gene);R<-matrix(0,length(ids),length(gems),dimnames=list(ids,gems));R[cbind(match(gr$gene,ids),match(gr$gem,gems))]<-gr$n_target/total[gr$gene]*(gr$contrast_z-yall[gr$gene]);wildall<-sweep(R%*%signs,1,yall,'+')
 for(n in names(cohorts)){
  g<-cohorts[[n]];y<-yall[g];wild<-wildall[g,,drop=FALSE];W<-readRDS(file.path(p,paste0('gene_bootstrap_weights_',n,'.rds')))[,g,drop=FALSE]
  for(m in unique(point$method))for(k in c(10L,20L)){
   obs<-point[axis==a&mode==mm&cohort==n&method==m&K==k];if(nrow(obs)!=1L)next
   if(m=='UPR_priority')bg<-old[[n]][list==if(n=='all438')'UPR'else'CC_UPR']else{label<-paste(a,mm,n,m,sep='|');bg<-new[list==label]};if(k==10)bg<-bg[Top10==TRUE]
   if(!obs$estimable){j<-j+1L;ans[[j]]<-cbind(obs,data.table(endpoint=c('signed_mean_gap','absolute_mean_gap','positive_gap'),estimate=NA_real_,gene_low=NA_real_,gene_high=NA_real_,matched_set_low=NA_real_,matched_set_high=NA_real_,technical_low=NA_real_,technical_high=NA_real_));next}
   sa<-unique(bg$selected_gene);sb<-bg$background_gene;ia<-match(sa,g);ib<-match(sb,g);stopifnot(length(sa)==k,length(sb)==5*k,!anyNA(c(ia,ib)))
   boot<-tech<-setboot<-matrix(NA_real_,B,3L,dimnames=list(NULL,c('signed_mean_gap','absolute_mean_gap','positive_gap')))
   for(e in colnames(boot)){
    v<-if(e=='signed_mean_gap')y else if(e=='absolute_mean_gap')abs(y)else as.numeric(y>0);names(v)<-g
    boot[,e]<-rowSums(sweep(W[,ia,drop=FALSE],2,v[ia],'*'))/rowSums(W[,ia,drop=FALSE])-rowSums(sweep(W[,ib,drop=FALSE],2,v[ib],'*'))/rowSums(W[,ib,drop=FALSE])
    vv<-if(e=='signed_mean_gap')wild else if(e=='absolute_mean_gap')abs(wild)else wild>0;tech[,e]<-colMeans(vv[ia,,drop=FALSE])-colMeans(vv[ib,,drop=FALSE])
    byset<-bg[,.(contrast=v[selected_gene[1]]-mean(v[background_gene])),by=selected_gene]$contrast;setboot[,e]<-replicate(B,mean(byset[sample.int(k,k,replace=TRUE)]))
   }
   observed<-c(signed_mean_gap=mean(y[sa])-mean(y[sb]),absolute_mean_gap=mean(abs(y[sa]))-mean(abs(y[sb])),positive_gap=mean(y[sa]>0)-mean(y[sb]>0));stopifnot(abs(observed[1]-obs$signed_mean_gap)<1e-10)
   j<-j+1L;ans[[j]]<-data.table(axis=a,mode=mm,cohort=n,method=m,K=k,estimable=TRUE,endpoint=names(observed),estimate=observed,gene_low=vapply(seq_len(3),function(i)ci(boot[,i])[1],0.),gene_high=vapply(seq_len(3),function(i)ci(boot[,i])[2],0.),matched_set_low=vapply(seq_len(3),function(i)ci(setboot[,i])[1],0.),matched_set_high=vapply(seq_len(3),function(i)ci(setboot[,i])[2],0.),technical_low=vapply(seq_len(3),function(i)ci(tech[,i])[1],0.),technical_high=vapply(seq_len(3),function(i)ci(tech[,i])[2],0.))
   draws[[paste(a,mm,n,m,k,sep='|')]]<-list(gene=boot,matched_set=setboot,technical=tech)
  }
 }
}
fwrite(rbindlist(ans,fill=TRUE),file.path(out,'alternative_fixed_background_gap_intervals.tsv'),sep='\t',na='NA');saveRDS(draws,file.path(out,'alternative_matched_gap_draws.rds'),compress='gzip');writeLines(capture.output(sessionInfo()),file.path(out,'matched_interval_sessionInfo.txt'));cat('Predeclared conditional intervals for existing matched-gap estimates complete; no new tests.\n')
