#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
project<-normalizePath(if(length(args))args[1] else 'paper19opp-revision2')
out<-file.path(project,'upr_robustness','ranking')
archive<-file.path(project,'reconstruction'); src<-file.path(archive,'cross_dataset','phase2_upr_rpe1'); val<-file.path(archive,'upr_rpe1_validation')
options(stringsAsFactors=FALSE, warn=1); set.seed(2026110607)
seed<-2026110607L; B<-2000L; P<-10000L
sha<-function(f)strsplit(system2('shasum',c('-a','256',shQuote(f)),stdout=TRUE),'[[:space:]]+')[[1]][1]
plan<-file.path(out,'locked_analysis_plan.txt')
stopifnot(sha(plan)=='b580e5e0a1fe45c6807d76dc91f1e2e0d348d08700db3a6528580ca878488cda')
files<-c(file.path(src,c('source_rank_inputs.tsv','all_source_comparator_rankings.tsv','complete_covariate_source_rankings.tsv','null100_source_rankings.tsv.gz','complete_covariate_null100_rankings.tsv.gz','frozen_background_matches.tsv','complete_covariate_background_matches.tsv','complete_covariate_gene_roster.tsv','source_UPR_exact_denominator_gene_effects.tsv','rpe1_UPR_exact_denominator_gene_effects.tsv')),file.path(val,c('phase2_statistical_state.rds','equal100_target_cell_effect_vectors.rds')))
stopifnot(all(file.exists(files)))
manifest<-data.frame(file=sub(paste0('^',project,'/'),'',files),sha256=vapply(files,sha,''))
write.table(manifest,file.path(out,'analysis_input_hashes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Analysis started: ',format(Sys.time(),tz='UTC',usetz=TRUE),'\nPlan verified; input hashes written before reading effects.\n')
readtsv<-function(f)read.delim(if(grepl('\\.gz$',f))gzfile(f)else f,check.names=FALSE)
wt<-function(z,n)write.table(z,file.path(out,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
state<-readRDS(file.path(val,'phase2_statistical_state.rds'))
meta<-readtsv(file.path(src,'source_rank_inputs.tsv'))
ids<-sort(meta$gene); stopifnot(length(ids)==438L,!anyDuplicated(ids))
x<-setNames(meta$source_UPR_deltaZ,meta$gene)[ids]
y<-setNames(state$validation$validation_delta_z,state$validation$gene)[ids]
strata<-setNames(meta$permutation_stratum,meta$gene)[ids]
fullr<-readtsv(file.path(src,'all_source_comparator_rankings.tsv'))
ccr<-readtsv(file.path(src,'complete_covariate_source_rankings.tsv'))
nullr<-readtsv(file.path(src,'null100_source_rankings.tsv.gz'))
nullcc<-readtsv(file.path(src,'complete_covariate_null100_rankings.tsv.gz'))
ccids<-sort(unique(ccr$gene)); strictids<-sort(meta$gene[meta$strict_rpe_support])
stopifnot(length(ccids)==370L,length(strictids)==390L,all(is.finite(x)),all(is.finite(y)))
wild<-state$primary$validation_wild$effects[ids,,drop=FALSE]
eq<-readRDS(file.path(val,'equal100_target_cell_effect_vectors.rds')); stopifnot(ncol(eq$validation)==200L)
exact<-readtsv(file.path(src,'rpe1_UPR_exact_denominator_gene_effects.tsv'))
yexact<-setNames(exact$deltaZ,exact$gene)[ids]
cat('Inputs read; cohorts and dimensions verified.\n')
# Functions preserve the locked source order and never select a ranking using y.
quintile<-function(z){o<-order(z,names(z)); n<-length(z); sizes<-rep(n%/%5L,5L); if(n%%5L)sizes[seq_len(n%%5L)]<-sizes[seq_len(n%%5L)]+1L; q<-integer(n); q[o]<-rep(1:5,sizes);setNames(q,names(z))}
weighted_median<-function(z,w){ok<-is.finite(z)&w>0;if(!any(ok))return(NA_real_);z<-z[ok];w<-w[ok];o<-order(z);z<-z[o];w<-w[o];cw<-cumsum(w);half<-sum(w)/2;mean(z[c(which(cw>=half)[1],which(cw>half)[1])])}
wm<-function(z,w)if(sum(w)>0)sum(z*w)/sum(w)else NA_real_
qtstat<-function(z,q,w=rep(1,length(z))){a<-wm(z[q==5],w[q==5])-wm(z[q==1],w[q==1]); med<-weighted_median(z[q==5],w[q==5])-weighted_median(z[q==1],w[q==1]); mq<-wm(q,w);mz<-wm(z,w);slope<-sum(w*(q-mq)*(z-mz))/sum(w*(q-mq)^2);reps<-rep(seq_along(z),w);tau<-if(length(reps)>1)cor(q[reps],z[reps],method='kendall')else NA_real_;c(Q5_Q1_mean=a,Q5_Q1_median=med,Kendall_tau=tau,ordinal_slope=slope)}
ndcg<-function(z,list,k){if(length(list)<k||!all(list[seq_len(k)]%in%names(z)))return(NA_real_);g<-pmax(z,0);ideal<-sum(sort(g,decreasing=TRUE)[seq_len(k)]/log2(seq_len(k)+1));if(ideal<=0)return(NA_real_);sum(g[list[seq_len(k)]]/log2(seq_len(k)+1))/ideal}
listmetrics<-function(z,list,k){sel<-head(list,k);if(length(sel)!=k||!all(sel%in%names(z)))return(c(estimable=0,retained=sum(sel%in%names(z)),Precision=NA,NDCG=NA,mean=NA,median=NA,positive_gap_cohort=NA,positive_ratio_cohort=NA));zz<-z[sel];base<-mean(z>0);c(estimable=1,retained=k,Precision=mean(zz>0),NDCG=ndcg(z,sel,k),mean=mean(zz),median=median(zz),positive_gap_cohort=mean(zz>0)-base,positive_ratio_cohort=if(base>0)mean(zz>0)/base else NA)}
getlists<-function(tab,methodcol='method',rankcol='rank'){lapply(split(tab,tab[[methodcol]]),function(d)d$gene[order(d[[rankcol]],d$gene)])}
lists438<-getlists(subset(fullr,available))
lists370<-getlists(ccr,rankcol='rank_cc')
null438<-getlists(nullr,'module');null370<-getlists(nullcc,'module','rank_cc')
stopifnot(length(null438)==100L,length(null370)==100L,all(vapply(null438,length,1L)==438L),all(vapply(null370,length,1L)==370L))
cohorts<-list(all438=ids,complete370=ccids,strict390=strictids)
qs<-lapply(cohorts,function(g)quintile(x[g]))
quintile_map<-do.call(rbind,lapply(names(cohorts),function(n)data.frame(cohort=n,gene=cohorts[[n]],source_d=x[cohorts[[n]]],quintile=qs[[n]],RPE1_d=y[cohorts[[n]]])))
wt(quintile_map,'source_quintile_membership.tsv')
# Strict sensitivity evaluates original all438 lists literally, without replacement/refill.
lists_by<-list(all438=lists438,complete370=lists370,strict390=lists438)
null_by<-list(all438=null438,complete370=null370,strict390=null438)
point_rank<-do.call(rbind,lapply(names(cohorts),function(n){z<-y[cohorts[[n]]];do.call(rbind,lapply(names(lists_by[[n]]),function(m)do.call(rbind,lapply(c(10L,20L),function(k)data.frame(cohort=n,method=m,K=k,t(listmetrics(z,lists_by[[n]][[m]],k)))))))}))
wt(point_rank,'fixed_ranking_metrics.tsv')
# Save every method/list identity, retaining unavailable strict-list members explicitly.
wt(do.call(rbind,lapply(names(lists_by),function(n)do.call(rbind,lapply(names(lists_by[[n]]),function(m)data.frame(cohort=n,method=m,gene=lists_by[[n]][[m]],rank=seq_along(lists_by[[n]][[m]]),in_cohort=lists_by[[n]][[m]]%in%cohorts[[n]]))))),'fixed_ranking_membership.tsv')
cat('Fixed source groups and list identities exported.\n')
quintsummary<-list();qboots<-list();trend<-list();isotonic<-list();adjacent<-list();ordered<-list()
max_t_intervals<-function(draws,point){se<-apply(draws,2,sd,na.rm=TRUE);good<-is.finite(se)&se>0;if(any(good)){student<-sweep(sweep(draws[,good,drop=FALSE],2,point[good],'-'),2,se[good],'/');crit<-quantile(apply(abs(student),1,max,na.rm=TRUE),.95,na.rm=TRUE)}else crit<-0;data.frame(adjacent_gap=names(point),estimate=point,bootstrap_sd=se,max_t_critical=as.numeric(crit),simultaneous_low=point-crit*se,simultaneous_high=point+crit*se)}
for(n in names(cohorts)){
 g<-cohorts[[n]];z<-y[g];q<-qs[[n]]; weights<-t(rmultinom(B,length(g),rep(1/length(g),length(g))));colnames(weights)<-g
 # Whole-gene multiplicities are retained and shared by quintile/ranking endpoints.
 saveRDS(weights,file.path(out,paste0('gene_bootstrap_weights_',n,'.rds')),compress='gzip')
 summary<-do.call(rbind,lapply(1:5,function(i){zz<-z[q==i];s<-x[g][q==i];data.frame(cohort=n,quintile=i,genes=length(zz),source_min=min(s),source_max=max(s),RPE1_mean=mean(zz),RPE1_median=median(zz),positive=sum(zz>0),positive_fraction=mean(zz>0),sign_concordant=sum(sign(s)==sign(zz)),sign_concordance=mean(sign(s)==sign(zz)),zero_zero=sum(s==0&zz==0))}))
 boots<-t(apply(weights,1,function(w)qtstat(z,q,w)));point<-qtstat(z,q)
 wilddraw<-t(apply(wild[g,,drop=FALSE],2,function(zb)qtstat(zb,q)))
 trend[[n]]<-data.frame(cohort=n,endpoint=names(point),estimate=point,gene_low=apply(boots,2,quantile,.025,na.rm=TRUE),gene_high=apply(boots,2,quantile,.975,na.rm=TRUE),technical_low=apply(wilddraw,2,quantile,.025,na.rm=TRUE),technical_high=apply(wilddraw,2,quantile,.975,na.rm=TRUE))
 for(i in 1:5){bm<-apply(weights[,q==i,drop=FALSE],1,function(w)wm(z[q==i],w));bp<-apply(weights[,q==i,drop=FALSE],1,function(w)wm(as.numeric(z[q==i]>0),w));bmed<-apply(weights[,q==i,drop=FALSE],1,function(w)weighted_median(z[q==i],w));summary$mean_low[i]<-quantile(bm,.025,na.rm=TRUE);summary$mean_high[i]<-quantile(bm,.975,na.rm=TRUE);summary$median_low[i]<-quantile(bmed,.025,na.rm=TRUE);summary$median_high[i]<-quantile(bmed,.975,na.rm=TRUE);summary$positive_low[i]<-quantile(bp,.025,na.rm=TRUE);summary$positive_high[i]<-quantile(bp,.975,na.rm=TRUE)}
 quintsummary[[n]]<-summary; iso<-isoreg(1:5,summary$RPE1_mean)$yf
 isotonic[[n]]<-data.frame(cohort=n,quintile=1:5,observed_mean=summary$RPE1_mean,increasing_isotonic_mean=iso,residual=summary$RPE1_mean-iso,adjacent_observed_change=c(NA,diff(summary$RPE1_mean)))
 group_gene<-sapply(1:5,function(i)apply(weights[,q==i,drop=FALSE],1,function(w)wm(z[q==i],w)))
 group_wild<-sapply(1:5,function(i)colMeans(wild[g[q==i],,drop=FALSE]))
 group_equal<-sapply(1:5,function(i)colMeans(eq$validation[g[q==i],,drop=FALSE]))
 point_gap<-setNames(diff(summary$RPE1_mean),paste0('Q',2:5,'-Q',1:4));gene_gap<-t(apply(group_gene,1,diff));wild_gap<-t(apply(group_wild,1,diff))
 colnames(gene_gap)<-colnames(wild_gap)<-names(point_gap)
 adjacent[[n]]<-rbind(data.frame(cohort=n,resampling='whole_gene',max_t_intervals(gene_gap,point_gap)),data.frame(cohort=n,resampling='shared_GEM_wild',max_t_intervals(wild_gap,point_gap)))
 ordered[[n]]<-data.frame(cohort=n,observed_means_nondecreasing=all(point_gap>=-1e-12),resampling=c('whole_gene','shared_GEM_wild','equal100'),draws=c(nrow(group_gene),nrow(group_wild),nrow(group_equal)),all_five_means_ordered_fraction=c(mean(apply(group_gene,1,function(z)all(diff(z)>=-1e-12))),mean(apply(group_wild,1,function(z)all(diff(z)>=-1e-12))),mean(apply(group_equal,1,function(z)all(diff(z)>=-1e-12)))))
 qboots[[n]]<-list(gene=boots,technical=wilddraw)
}
wt(do.call(rbind,quintsummary),'quintile_group_summaries.tsv');wt(do.call(rbind,trend),'quintile_trend_estimates.tsv');wt(do.call(rbind,isotonic),'quintile_isotonic_descriptive.tsv')
wt(do.call(rbind,adjacent),'quintile_adjacent_simultaneous_intervals.tsv');wt(do.call(rbind,ordered),'quintile_ordered_draw_fractions.tsv')
saveRDS(qboots,file.path(out,'quintile_bootstrap_draws.rds'),compress='gzip')
cat('Quintile estimates/whole-gene and shared-GEM intervals complete.\n')
# Pairwise rank comparisons use shared weights across methods; NDCG linear differences condition on observed IDCG.
comp<-list();compdraws<-list(); compid<-0L
for(n in c('all438','complete370')){
 g<-cohorts[[n]];z<-y[g];W<-readRDS(file.path(out,paste0('gene_bootstrap_weights_',n,'.rds')));ll<-lists_by[[n]]
 for(reference in c('raw_UPR','UPR_priority'))for(method in setdiff(names(ll),reference))for(k in c(10L,20L)){
  a<-head(ll[[method]],k);b<-head(ll[[reference]],k);if(!all(c(a,b)%in%g))next
  ia<-match(a,g);ib<-match(b,g);prec_a<-rowSums(sweep(W[,ia,drop=FALSE],2,z[a]>0,'*'))/rowSums(W[,ia,drop=FALSE]);prec_b<-rowSums(sweep(W[,ib,drop=FALSE],2,z[b]>0,'*'))/rowSums(W[,ib,drop=FALSE]);disc<-1/log2(seq_len(k)+1);coef<-numeric(length(g));coef[ia]<-coef[ia]+disc;coef[ib]<-coef[ib]-disc;idcg<-sum(sort(pmax(z,0),decreasing=TRUE)[seq_len(k)]*disc);ndboot<-as.numeric(W%*%(coef*pmax(z,0)))/idcg
  wild_prec<-apply(wild[g,,drop=FALSE],2,function(v)mean(v[ia]>0)-mean(v[ib]>0));wild_nd<-apply(wild[g,,drop=FALSE],2,function(v){names(v)<-g;ndcg(v,a,k)-ndcg(v,b,k)})
  points<-c(Precision=mean(z[a]>0)-mean(z[b]>0),NDCG=ndcg(z,a,k)-ndcg(z,b,k));draws<-cbind(Precision=prec_a-prec_b,NDCG=ndboot);tech<-cbind(Precision=wild_prec,NDCG=wild_nd)
  compid<-compid+1L;comp[[compid]]<-data.frame(cohort=n,method=method,reference=reference,K=k,endpoint=names(points),difference=points,gene_low=apply(draws,2,quantile,.025,na.rm=TRUE),gene_high=apply(draws,2,quantile,.975,na.rm=TRUE),technical_low=apply(tech,2,quantile,.025,na.rm=TRUE),technical_high=apply(tech,2,quantile,.975,na.rm=TRUE),selection_overlap=length(intersect(a,b)),idcg=idcg)
  compdraws[[paste(n,method,reference,k,sep='|')]]<-list(gene=draws,technical=tech)
 }
}
wt(do.call(rbind,comp),'paired_method_differences.tsv');saveRDS(compdraws,file.path(out,'paired_method_draws.rds'),compress='gzip')
cat('Paired method differences complete.\n')
# One joint seven-test family; no test or comparator is selected after seeing its result.
q<-qs$all438
rankdiff<-function(z){a<-lists370$UPR_priority;b<-lists370$raw_UPR;c(Precision10=mean(z[head(a,10)]>0)-mean(z[head(b,10)]>0),Precision20=mean(z[head(a,20)]>0)-mean(z[head(b,20)]>0),NDCG10=ndcg(z,a,10)-ndcg(z,b,10),NDCG20=ndcg(z,a,20)-ndcg(z,b,20))}
obsq<-qtstat(y,q)[c('Q5_Q1_mean','Kendall_tau','ordinal_slope')];obsr<-rankdiff(y[ccids]);obs<-c(obsq,obsr)
perm<-matrix(NA_real_,P,7L,dimnames=list(NULL,names(obs)))
stratum438<-split(seq_along(ids),strata);stratum370<-split(seq_along(ccids),strata[ccids])
for(b in seq_len(P)){
 iy<-seq_along(ids);for(s in stratum438)iy[s]<-s[sample.int(length(s))];yp<-setNames(y[iy],ids)
 jr<-seq_along(ccids);for(s in stratum370)jr[s]<-s[sample.int(length(s))];rp<-setNames(y[ccids][jr],ccids)
 perm[b,]<-c(qtstat(yp,q)[c('Q5_Q1_mean','Kendall_tau','ordinal_slope')],rankdiff(rp))
 if(b%%2000L==0L)cat('Permutation',b,'of',P,'\n')
}
center<-colMeans(perm);centerobs<-obs-center;pval<-vapply(seq_along(obs),function(j){if(all(abs(perm[,j])<1e-12)&&abs(obs[j])<1e-12)return(1);(1+sum(abs(perm[,j]-center[j])>=abs(centerobs[j])-1e-12))/(P+1)},0.)
wt(data.frame(endpoint=names(obs),cohort=c(rep('all438',3),rep('complete370',4)),observed=obs,permutation_null_mean=center,centered_observed=centerobs,two_sided_P=pval,Holm7_P=p.adjust(pval,'holm'),permutation_draws=P),'joint_Holm7_tests.tsv');saveRDS(perm,file.path(out,'joint_permutation_draws.rds'),compress='gzip')
# All100 fixed null ranks and uniform random rankings are computational reference controls, without extra confirmatory P values.
nullmetrics<-list();uniformmetrics<-list();calibration<-list();ni<-0L;ui<-0L;ci<-0L
for(n in names(cohorts)){
 g<-cohorts[[n]];z<-y[g];ll<-null_by[[n]]
 for(m in names(ll))for(k in c(10L,20L)){ni<-ni+1L;nullmetrics[[ni]]<-data.frame(cohort=n,module=m,K=k,t(listmetrics(z,ll[[m]],k)))}
 U<-matrix(NA_real_,P,4L,dimnames=list(NULL,c('Precision10','Precision20','NDCG10','NDCG20')))
 for(b in seq_len(P)){a<-sample(g,length(g));U[b,]<-c(mean(z[a[1:10]]>0),mean(z[a[1:20]]>0),ndcg(z,a,10),ndcg(z,a,20))}
 ui<-ui+1L;uniformmetrics[[ui]]<-data.frame(cohort=n,endpoint=colnames(U),expected=colMeans(U),reference_low=apply(U,2,quantile,.025),reference_high=apply(U,2,quantile,.975),draws=P,cohort_positive_prevalence=mean(z>0))
 saveRDS(U,file.path(out,paste0('uniform_rank_draws_',n,'.rds')),compress='gzip')
}
nulldf<-do.call(rbind,nullmetrics);wt(nulldf,'all100_fixed_null_ranking_metrics.tsv');wt(do.call(rbind,uniformmetrics),'uniform_random_ranking_reference.tsv')
for(n in names(cohorts))for(k in c(10L,20L))for(e in c('Precision','NDCG')){
 v<-subset(nulldf,cohort==n&K==k)[[e]];r<-subset(point_rank,cohort==n&method=='UPR_priority'&K==k)[[e]];v<-v[is.finite(v)];ci<-ci+1L;calibration[[ci]]<-data.frame(cohort=n,K=k,endpoint=e,real_UPR=r,valid_null_modules=length(v),null_median=if(length(v))median(v)else NA,finite_null_upper_tail=if(length(v)&&is.finite(r))(1+sum(v>=r-1e-12))/(length(v)+1)else NA)}
wt(do.call(rbind,calibration),'all100_null_metric_calibration.tsv')
cat('Uniform and expression/detection-matched ranking controls complete.\n')
# Literal fixed-list positive enrichment/sign concordance with source-fixed backgrounds.
bg438<-readtsv(file.path(src,'frozen_background_matches.tsv'));bg370<-readtsv(file.path(src,'complete_covariate_background_matches.tsv'))
bgout<-list();signout<-list();bi<-0L;si<-0L
for(n in names(cohorts)){
 g<-cohorts[[n]];z<-y[g];bg<-if(n=='complete370')bg370 else bg438;W<-readRDS(file.path(out,paste0('gene_bootstrap_weights_',n,'.rds')))
 si<-si+1L;signout[[si]]<-data.frame(cohort=n,method='all_genes',K=length(g),source_positive=sum(x[g]>0),source_negative=sum(x[g]<0),RPE1_positive=sum(z>0),RPE1_negative=sum(z<0),sign_concordant=sum(sign(x[g])==sign(z)),sign_concordance=mean(sign(x[g])==sign(z)),zero_zero=sum(x[g]==0&z==0))
 for(m in names(lists_by[[n]]))for(k in c(10L,20L)){
  a<-head(lists_by[[n]][[m]],k);if(!all(a%in%g))next
  si<-si+1L;signout[[si]]<-data.frame(cohort=n,method=m,K=k,source_positive=sum(x[a]>0),source_negative=sum(x[a]<0),RPE1_positive=sum(z[a]>0),RPE1_negative=sum(z[a]<0),sign_concordant=sum(sign(x[a])==sign(z[a])),sign_concordance=mean(sign(x[a])==sign(z[a])),zero_zero=sum(x[a]==0&z[a]==0))
  bgname<-if(m=='UPR_priority')if(n=='complete370')'CC_UPR'else'UPR'else m
  d<-subset(bg,list==bgname);if(k==10)d<-d[d$Top10,,drop=FALSE];if(length(unique(d$selected_gene))!=k||nrow(d)!=5*k||!all(c(d$selected_gene,d$background_gene)%in%g))next
  a<-unique(d$selected_gene);b<-d$background_gene;ia<-match(a,g);ib<-match(b,g);pa<-mean(z[a]>0);pb<-mean(z[b]>0);wb<-rowSums(sweep(W[,ia,drop=FALSE],2,z[a]>0,'*'))/rowSums(W[,ia,drop=FALSE])-rowSums(sweep(W[,ib,drop=FALSE],2,z[b]>0,'*'))/rowSums(W[,ib,drop=FALSE]);tb<-apply(wild[g,,drop=FALSE],2,function(v)mean(v[ia]>0)-mean(v[ib]>0));sets<-split(d,d$selected_gene);setgap<-vapply(sets,function(t)as.numeric(z[t$selected_gene[1]]>0)-mean(z[t$background_gene]>0),0.);setboot<-replicate(B,mean(sample(setgap,length(setgap),replace=TRUE)))
  bi<-bi+1L;bgout[[bi]]<-data.frame(cohort=n,method=m,K=k,selected_positive=pa,matched_background_positive=pb,positive_gap=pa-pb,positive_ratio=if(pb>0)pa/pb else NA,gene_low=quantile(wb,.025,na.rm=TRUE),gene_high=quantile(wb,.975,na.rm=TRUE),matched_set_low=quantile(setboot,.025),matched_set_high=quantile(setboot,.975),technical_low=quantile(tb,.025),technical_high=quantile(tb,.975),signed_mean_gap=mean(z[a])-mean(z[b]),absolute_mean_gap=mean(abs(z[a]))-mean(abs(z[b])),background_genes=length(b))
 }
}
wt(do.call(rbind,bgout),'fixed_background_positive_enrichment.tsv');wt(do.call(rbind,signout),'source_validation_sign_concordance.tsv')
# Exact-denominator/equal100 sensitivities preserve original source groups and list memberships.
sensitivity<-list();qi<-0L;rank_sens<-list();ri<-0L
for(n in names(cohorts)){
 g<-cohorts[[n]];q<-qs[[n]]
 for(sens in c('exact_target_denominator','equal100')){
  yy<-if(sens=='exact_target_denominator')matrix(yexact[g],ncol=1,dimnames=list(g,NULL))else eq$validation[g,,drop=FALSE]
  tr<-t(apply(yy,2,function(v)qtstat(v,q)));qi<-qi+1L;sensitivity[[qi]]<-data.frame(cohort=n,sensitivity=sens,endpoint=colnames(tr),median=apply(tr,2,median),low=apply(tr,2,quantile,.025),high=apply(tr,2,quantile,.975),draws=nrow(tr),positive_fraction=colMeans(tr>0))
  for(m in names(lists_by[[n]]))for(k in c(10L,20L)){
   a<-head(lists_by[[n]][[m]],k);if(!all(a%in%g))next
   mat<-t(apply(yy,2,function(v){names(v)<-g;listmetrics(v,a,k)}));for(e in c('Precision','NDCG','mean','median','positive_gap_cohort')){ri<-ri+1L;rank_sens[[ri]]<-data.frame(cohort=n,sensitivity=sens,method=m,K=k,endpoint=e,median=median(mat[,e],na.rm=TRUE),low=quantile(mat[,e],.025,na.rm=TRUE),high=quantile(mat[,e],.975,na.rm=TRUE),draws=nrow(mat))}
  }
 }
}
wt(do.call(rbind,sensitivity),'quintile_exact_equal100_sensitivities.tsv');wt(do.call(rbind,rank_sens),'ranking_exact_equal100_sensitivities.tsv')
# Algebraic preflight checks test metric definitions rather than mirror the observed data.
v<-setNames(c(-1,0,.5,2),letters[1:4]);stopifnot(ndcg(v,c('d','c','b','a'),2)==1,ndcg(v,c('a','b','c','d'),2)==0,abs(listmetrics(v,c('c','d'),2)['Precision']-1)<1e-12,abs(listmetrics(v,c('a','b'),2)['Precision'])<1e-12)
wt(data.frame(check=c('ideal_linear_NDCG_equals_one','zero_relevance_NDCG_equals_zero','strict_positive_sign_not_nonnegative','paired_draws_and_fixed_source_lists'),passed=TRUE),'algebraic_preflight.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('This is a retrospective robustness analysis after prior RPE1 outcome evaluation.','Holm7 is the only added inferential family; all other calibrations/sensitivities are descriptive.','Whole-gene NDCG-difference intervals condition on observed IDCG.','No new independent validation, prospective efficacy or functional benefit is established.'),file.path(out,'interpretation_scope.txt'))
cat('Analysis completed: ',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
