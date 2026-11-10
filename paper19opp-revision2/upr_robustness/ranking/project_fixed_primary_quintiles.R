#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(if(length(args))args[1]else'paper19opp-revision2');p<-file.path(project,'upr_robustness','ranking');src<-file.path(project,'reconstruction','cross_dataset','phase2_upr_rpe1');val<-file.path(project,'reconstruction','upr_rpe1_validation');options(stringsAsFactors=FALSE)
# This executes only the already hashed metric function definitions, not the main analysis.
needed<-c('weighted_median','wm','qtstat','max_t_intervals')
for(e in parse(file.path(p,'analyze_rank_robustness.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%needed)eval(e)
r<-function(n)read.delim(file.path(p,n),check.names=FALSE);wt<-function(z,n)write.table(z,file.path(p,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
map<-r('source_quintile_membership.tsv');full<-subset(map,cohort=='all438');qall<-setNames(full$quintile,full$gene)
state<-readRDS(file.path(val,'phase2_statistical_state.rds'));wild<-state$primary$validation_wild$effects;eq<-readRDS(file.path(val,'equal100_target_cell_effect_vectors.rds'))
ex<-read.delim(file.path(src,'rpe1_UPR_exact_denominator_gene_effects.tsv'));yex<-setNames(ex$deltaZ,ex$gene)
sx<-read.delim(file.path(src,'source_UPR_exact_denominator_gene_effects.tsv'));xex<-setNames(sx$deltaZ,sx$gene)
summary<-list();trends<-list();adj<-list();orderfreq<-list();mem<-list();sens<-list();signsens<-list();j<-0L
for(n in c('complete370','strict390')){
 d<-subset(map,cohort==n);g<-d$gene;z<-setNames(d$RPE1_d,g);q<-qall[g];W<-readRDS(file.path(p,paste0('gene_bootstrap_weights_',n,'.rds')))[,g,drop=FALSE]
 mem[[n]]<-data.frame(cohort=n,gene=g,primary438_quintile=q,within_cohort_quintile=d$quintile)
 s<-do.call(rbind,lapply(1:5,function(i){zz<-z[q==i];means<-apply(W[,q==i,drop=FALSE],1,function(w)wm(zz,w));prop<-apply(W[,q==i,drop=FALSE],1,function(w)wm(as.numeric(zz>0),w));data.frame(cohort=n,source_groups='fixed_primary438',quintile=i,genes=length(zz),RPE1_mean=mean(zz),RPE1_median=median(zz),positive=sum(zz>0),positive_fraction=mean(zz>0),mean_low=quantile(means,.025,na.rm=TRUE),mean_high=quantile(means,.975,na.rm=TRUE),positive_low=quantile(prop,.025,na.rm=TRUE),positive_high=quantile(prop,.975,na.rm=TRUE))}))
 summary[[n]]<-s;point<-qtstat(z,q);boots<-t(apply(W,1,function(w)qtstat(z,q,w)));technical<-t(apply(wild[g,,drop=FALSE],2,function(v)qtstat(v,q)))
 trends[[n]]<-data.frame(cohort=n,source_groups='fixed_primary438',endpoint=names(point),estimate=point,gene_low=apply(boots,2,quantile,.025),gene_high=apply(boots,2,quantile,.975),technical_low=apply(technical,2,quantile,.025),technical_high=apply(technical,2,quantile,.975))
 gm<-sapply(1:5,function(i)apply(W[,q==i,drop=FALSE],1,function(w)wm(z[q==i],w)));tm<-sapply(1:5,function(i)colMeans(wild[g[q==i],,drop=FALSE]));em<-sapply(1:5,function(i)colMeans(eq$validation[g[q==i],,drop=FALSE]));pg<-setNames(diff(s$RPE1_mean),paste0('Q',2:5,'-Q',1:4));gd<-t(apply(gm,1,diff));td<-t(apply(tm,1,diff));colnames(gd)<-colnames(td)<-names(pg)
 adj[[n]]<-rbind(data.frame(cohort=n,source_groups='fixed_primary438',resampling='whole_gene',max_t_intervals(gd,pg)),data.frame(cohort=n,source_groups='fixed_primary438',resampling='shared_GEM_wild',max_t_intervals(td,pg)))
 orderfreq[[n]]<-data.frame(cohort=n,source_groups='fixed_primary438',observed_means_nondecreasing=all(pg>=-1e-12),resampling=c('whole_gene','shared_GEM_wild','equal100'),draws=c(nrow(gm),nrow(tm),nrow(em)),all_five_means_ordered_fraction=c(mean(apply(gm,1,function(v)all(diff(v)>=-1e-12))),mean(apply(tm,1,function(v)all(diff(v)>=-1e-12))),mean(apply(em,1,function(v)all(diff(v)>=-1e-12)))))
 for(a in c('exact_target_denominator','equal100')){yy<-if(a=='exact_target_denominator')matrix(yex[g],ncol=1)else eq$validation[g,,drop=FALSE];dd<-t(apply(yy,2,function(v)qtstat(v,q)));j<-j+1L;sens[[j]]<-data.frame(cohort=n,source_groups='fixed_primary438',sensitivity=a,endpoint=colnames(dd),median=apply(dd,2,median),low=apply(dd,2,quantile,.025),high=apply(dd,2,quantile,.975),draws=nrow(dd))}
}
wt(do.call(rbind,mem),'primary_quintile_projection_membership.tsv');wt(do.call(rbind,summary),'primary_quintile_projection_group_summaries.tsv');wt(do.call(rbind,trends),'primary_quintile_projection_trend_estimates.tsv');wt(do.call(rbind,adj),'primary_quintile_projection_adjacent_intervals.tsv');wt(do.call(rbind,orderfreq),'primary_quintile_projection_ordered_fractions.tsv');wt(do.call(rbind,sens),'primary_quintile_projection_exact_equal100.tsv')
# Sign agreement uses each sensitivity's own source and validation effects, while list memberships stay fixed.
lists<-r('fixed_ranking_membership.tsv');met<-state$source_metadata
si<-0L;for(n in c('all438','complete370','strict390')){g<-subset(map,cohort==n)$gene;roster<-list(all_genes=g);ll<-subset(lists,cohort==n);for(m in unique(ll$method))for(k in c(10L,20L)){a<-head(ll$gene[ll$method==m][order(ll$rank[ll$method==m])],k);if(all(a%in%g))roster[[paste(m,k,sep='|')]]<-a}
 for(l in names(roster)){a<-roster[[l]];ec<-mean(sign(xex[a])==sign(yex[a]));eqc<-colMeans(sign(eq$source[a,,drop=FALSE])==sign(eq$validation[a,,drop=FALSE]));si<-si+1L;signsens[[si]]<-data.frame(cohort=n,list=l,genes=length(a),exact_denominator_concordance=ec,equal100_median=median(eqc),equal100_low=quantile(eqc,.025),equal100_high=quantile(eqc,.975),draws=length(eqc))}}
wt(do.call(rbind,signsens),'sign_concordance_exact_equal100.tsv')
cat('Predeclared primary438 quintile projections and sign sensitivities completed; no formal tests or original outputs changed.\n');writeLines(capture.output(sessionInfo()),file.path(p,'projection_sessionInfo.txt'))
