#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(if(length(args))args[1]else'paper19opp-revision2');p<-file.path(project,'upr_robustness','ranking');out<-file.path(p,'alternative_definitions');modules<-file.path(project,'upr_robustness','modules');src<-file.path(project,'reconstruction','cross_dataset','phase2_upr_rpe1');.libPaths(c(file.path(dirname(project),'.Rlib'),.libPaths()));suppressPackageStartupMessages({library(data.table);library(digest)})
options(stringsAsFactors=FALSE,warn=1);set.seed(2026110608);B<-2000L;P<-10000L
# Verify source choices before opening any newly scored validation effect table.
manifest<-fread(file.path(out,'source_lock_manifest.tsv'));stopifnot(digest(file.path(out,'source_lock_manifest.tsv'),file=TRUE,algo='sha256')=='b94e9a421c1e6e1a7ceee9829bbe894a7c88bd2ee685438044da27f2680007e2');stopifnot(all(vapply(file.path(project,manifest$file),digest,'',file=TRUE,algo='sha256')==manifest$sha256))
inputs<-fread(file.path(out,'source_input_hashes.tsv'));stopifnot(all(vapply(file.path(project,inputs$file),digest,'',file=TRUE,algo='sha256')==inputs$sha256))
qa<-fread(file.path(modules,'scored','RPE1_Hallmark_reference_QA.tsv'));stopifnot(nrow(qa)==2L,all(qa$passed))
validation_files<-file.path(modules,'scored',c('RPE1_gene_effects.tsv.gz','RPE1_GEM_effects.tsv.gz','RPE1_Hallmark_reference_QA.tsv'))
fwrite(data.table(file=sub(paste0('^',project,'/'),'',validation_files),sha256=vapply(validation_files,digest,'',file=TRUE,algo='sha256')),file.path(out,'validation_input_hashes.tsv'),sep='\t')
writeLines(paste('First alternative-definition RPE1 effect read UTC',format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'validation_first_read.txt'))
cat('Source manifest/input hashes verified; opening released validation effects',format(Sys.time(),tz='UTC',usetz=TRUE),'\n');flush.console()
d<-fread(validation_files[1]);gg<-fread(validation_files[2]);stopifnot(all(d$study=='RPE1'),all(gg$study=='RPE1'))
# Only previously specified metric functions are evaluated from the original extension.
needed<-c('weighted_median','wm','qtstat','ndcg','listmetrics','getlists')
for(e in parse(file.path(p,'analyze_rank_robustness.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%needed)eval(e)
ranks<-fread(file.path(out,'source_alternative_rankings.tsv.gz'));nullr<-fread(file.path(out,'source_parent_specific_null_rankings.tsv.gz'));source_real<-fread(file.path(modules,'scored','K562_gene_effects.tsv.gz'));meta<-fread(file.path(src,'source_rank_inputs.tsv'));ids<-sort(meta$gene);ccids<-sort(fread(file.path(src,'complete_covariate_gene_roster.tsv'))$gene)
orig438<-getlists(fread(file.path(src,'all_source_comparator_rankings.tsv'))[available==TRUE]);orig370<-getlists(fread(file.path(src,'complete_covariate_source_rankings.tsv')),rankcol='rank_cc')
cohorts<-list(all438=ids,complete370=ccids);original<-list(all438=orig438,complete370=orig370)
axes<-c('Reactome_UPR_full','Reactome_UPR_no_Hallmark','GO_ER_UPR_full','GO_ER_UPR_no_Hallmark');modes<-c('exact_Lminus_target','fixed_L6780')
oldstate<-readRDS(file.path(project,'reconstruction','upr_rpe1_validation','phase2_statistical_state.rds'));signs<-oldstate$primary$validation_wild$wild_signs;gemlabels<-oldstate$primary$validation_wild$gem_labels
qmap<-fread(file.path(p,'source_quintile_membership.tsv'))[cohort=='all438'];qall<-setNames(qmap$quintile,qmap$gene)
newmatches<-fread(file.path(out,'source_alternative_background_matches.tsv'));oldmatches<-list(all438=fread(file.path(src,'frozen_background_matches.tsv')),complete370=fread(file.path(src,'complete_covariate_background_matches.tsv')))
point<-list();diffs<-list();nullmetrics<-list();cal<-list();uniform<-list();bggaps<-list();quintiles<-list();pi<-di<-ni<-ci<-ui<-bi<-qi<-0L
wild_effects<-function(gr,y,g){gr<-gr[gene%in%g];chk<-gr[,.(point=weighted.mean(contrast_z,n_target),n=sum(n_target)),by=gene];stopifnot(setequal(chk$gene,g),max(abs(setNames(chk$point,chk$gene)[g]-y))<1e-8);total<-setNames(chk$n,chk$gene);R<-matrix(0,length(g),length(gemlabels),dimnames=list(g,gemlabels));stopifnot(all(gr$gem%in%gemlabels));v<-gr$n_target/total[gr$gene]*(gr$contrast_z-y[gr$gene]);R[cbind(match(gr$gene,g),match(gr$gem,gemlabels))]<-v;sweep(R%*%signs,1,y,'+')}
wt<-function(z,n)fwrite(z,file.path(out,n),sep='\t',na='NA')
for(a in axes)for(mode_name in modes){
 dd<-d[axis==a&mode==mode_name&gene%in%ids];stopifnot(nrow(dd)==438L,all(dd$estimable),all(is.finite(dd$delta_z)));yy<-setNames(dd$delta_z,dd$gene)[ids]
 gs<-gg[axis==a&mode==mode_name&gene%in%ids];wy<-wild_effects(gs,yy,ids)
 sx<-source_real[axis==a&mode==mode_name&gene%in%ids];xx<-setNames(sx$delta_z,sx$gene)
 for(n in names(cohorts)){
  g<-cohorts[[n]];y<-yy[g];W<-readRDS(file.path(p,paste0('gene_bootstrap_weights_',n,'.rds')))[,g,drop=FALSE];wild<-wy[g,,drop=FALSE]
  ll<-original[[n]];ll<-c(ll,getlists(ranks[axis==a&mode==mode_name&cohort==n],methodcol='kind'))
  ll<-lapply(ll,as.character);stopifnot(all(vapply(ll,function(v)all(head(v,20)%in%g),TRUE)))
  for(m in names(ll))for(k in c(10L,20L)){
   sel<-head(ll[[m]],k);me<-listmetrics(y,sel,k);ia<-match(sel,g);pb<-rowSums(sweep(W[,ia,drop=FALSE],2,y[sel]>0,'*'))/rowSums(W[,ia,drop=FALSE]);tp<-colMeans(wild[ia,,drop=FALSE]>0);tn<-apply(wild,2,function(v){names(v)<-g;ndcg(v,sel,k)})
   pi<-pi+1L;point[[pi]]<-data.table(axis=a,mode=mode_name,cohort=n,method=m,K=k,estimable=me['estimable'],Precision=me['Precision'],NDCG=me['NDCG'],mean_signed=me['mean'],median_signed=me['median'],mean_absolute=mean(abs(y[sel])),cohort_positive=mean(y>0),positive_gap_cohort=me['positive_gap_cohort'],positive_ratio_cohort=me['positive_ratio_cohort'],source_RPE_same_readout_sign_concordance=mean(sign(xx[sel])==sign(y[sel])),Precision_gene_low=quantile(pb,.025,na.rm=TRUE),Precision_gene_high=quantile(pb,.975,na.rm=TRUE),Precision_technical_low=quantile(tp,.025),Precision_technical_high=quantile(tp,.975),NDCG_technical_low=quantile(tn,.025),NDCG_technical_high=quantile(tn,.975))
   if(m=='UPR_priority'){bg<-oldmatches[[n]][list==if(n=='all438')'UPR'else'CC_UPR']}else if(m%in%c('alternative_priority','alternative_raw')){label<-paste(a,mode_name,n,m,sep='|');bg<-newmatches[list==label]}else bg<-NULL
   if(!is.null(bg)){if(k==10)bg<-bg[Top10==TRUE];valid<-nrow(bg)==5*k&&length(unique(bg$selected_gene))==k&&all(c(bg$selected_gene,bg$background_gene)%in%g);bi<-bi+1L;if(valid){sb<-unique(bg$selected_gene);bb<-bg$background_gene;bggaps[[bi]]<-data.table(axis=a,mode=mode_name,cohort=n,method=m,K=k,estimable=TRUE,signed_mean_gap=mean(y[sb])-mean(y[bb]),absolute_mean_gap=mean(abs(y[sb]))-mean(abs(y[bb])),positive_gap=mean(y[sb]>0)-mean(y[bb]>0),selected_positive=mean(y[sb]>0),background_positive=mean(y[bb]>0))}else bggaps[[bi]]<-data.table(axis=a,mode=mode_name,cohort=n,method=m,K=k,estimable=FALSE,signed_mean_gap=NA_real_,absolute_mean_gap=NA_real_,positive_gap=NA_real_,selected_positive=NA_real_,background_positive=NA_real_)}
  }
  # New definition score compared with every unchanged method and its own raw-source score.
  aa<-ll$alternative_priority
  for(reference in setdiff(names(ll),'alternative_priority'))for(k in c(10L,20L)){
   sa<-head(aa,k);sb<-head(ll[[reference]],k);ia<-match(sa,g);ib<-match(sb,g);pa<-rowSums(sweep(W[,ia,drop=FALSE],2,y[sa]>0,'*'))/rowSums(W[,ia,drop=FALSE]);pb<-rowSums(sweep(W[,ib,drop=FALSE],2,y[sb]>0,'*'))/rowSums(W[,ib,drop=FALSE]);disc<-1/log2(seq_len(k)+1);coef<-numeric(length(g));coef[ia]<-coef[ia]+disc;coef[ib]<-coef[ib]-disc;ideal<-sum(sort(pmax(y,0),decreasing=TRUE)[seq_len(k)]*disc);nb<-as.numeric(W%*%(coef*pmax(y,0)))/ideal
   tp<-colMeans(wild[ia,,drop=FALSE]>0)-colMeans(wild[ib,,drop=FALSE]>0);tn<-apply(wild,2,function(v){names(v)<-g;ndcg(v,sa,k)-ndcg(v,sb,k)})
   di<-di+1L;diffs[[di]]<-data.table(axis=a,mode=mode_name,cohort=n,method='alternative_priority',reference=reference,K=k,endpoint=c('Precision','NDCG'),difference=c(mean(y[sa]>0)-mean(y[sb]>0),ndcg(y,sa,k)-ndcg(y,sb,k)),gene_low=c(quantile(pa-pb,.025,na.rm=TRUE),quantile(nb,.025,na.rm=TRUE)),gene_high=c(quantile(pa-pb,.975,na.rm=TRUE),quantile(nb,.975,na.rm=TRUE)),technical_low=c(quantile(tp,.025),quantile(tn,.025)),technical_high=c(quantile(tp,.975),quantile(tn,.975)),selection_overlap=length(intersect(sa,sb)),conditional_observed_IDCG=ideal)
  }
  # Current438-source quintile projection remains fixed for alternative RPE1 readouts.
  q<-qall[g];for(i in 1:5){z<-y[q==i];qi<-qi+1L;quintiles[[qi]]<-data.table(axis=a,mode=mode_name,cohort=n,source_groups='current_Hallmark438',quintile=i,genes=length(z),mean_signed=mean(z),median_signed=median(z),positive_fraction=mean(z>0))}
  U<-matrix(NA_real_,P,4,dimnames=list(NULL,c('Precision10','Precision20','NDCG10','NDCG20')))
  for(b in seq_len(P)){rr<-sample(g,length(g));U[b,]<-c(mean(y[rr[1:10]]>0),mean(y[rr[1:20]]>0),ndcg(y,rr,10),ndcg(y,rr,20))};ui<-ui+1L;uniform[[ui]]<-data.table(axis=a,mode=mode_name,cohort=n,endpoint=colnames(U),mean=colMeans(U),low=apply(U,2,quantile,.025),high=apply(U,2,quantile,.975),draws=P)
  if(mode_name=='fixed_L6780'){
   nn<-nullr[parent==a&cohort==n];nl<-getlists(nn,methodcol='axis')
   for(mod in names(nl))for(k in c(10L,20L)){mm<-listmetrics(y,nl[[mod]],k);ni<-ni+1L;nullmetrics[[ni]]<-data.table(parent=a,mode=mode_name,cohort=n,module=mod,K=k,Precision=mm['Precision'],NDCG=mm['NDCG'],mean_signed=mm['mean'])}
  }
 }
 cat('Completed alternative readout',a,mode_name,'\n');flush.console()
}
point<-rbindlist(point);nullmetrics<-rbindlist(nullmetrics)
for(a in axes)for(n in names(cohorts))for(method_name in c('UPR_priority','alternative_priority','alternative_raw'))for(k in c(10L,20L))for(e in c('Precision','NDCG','mean_signed')){vv<-nullmetrics[parent==a&cohort==n&K==k][[e]];rr<-point[axis==a&mode=='fixed_L6780'&cohort==n&method==method_name&K==k][[e]];ci<-ci+1L;cal[[ci]]<-data.table(axis=a,mode='fixed_L6780',cohort=n,method=method_name,K=k,endpoint=e,observed=rr,valid_parent_specific_null_modules=sum(is.finite(vv)),finite_null_upper_tail=(1+sum(vv>=rr-1e-12,na.rm=TRUE))/(1+sum(is.finite(vv))),null_median=median(vv,na.rm=TRUE))}
wt(point,'alternative_readout_fixed_list_metrics.tsv');wt(rbindlist(diffs),'alternative_priority_paired_comparisons.tsv');wt(nullmetrics,'parent_specific_null_list_metrics.tsv');wt(rbindlist(cal),'parent_specific_fixed_mode_null_calibration.tsv');wt(rbindlist(uniform),'alternative_uniform_random_reference.tsv');wt(rbindlist(bggaps,fill=TRUE),'alternative_fixed_background_gaps.tsv');wt(rbindlist(quintiles),'alternative_readout_current_quintile_projection.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'evaluation_sessionInfo.txt'));writeLines(c('All results are retrospective definition/mode sensitivity estimates; no new formal P-value family was added.','Current lists and source-frozen alternative lists were evaluated on every specified readout/mode.','Exact-primary real readouts have no exact-mode parent-null calibration; parent-specific null comparisons use fixed_L6780 only.','Whole-gene NDCG-difference intervals condition on observed same-cohort IDCG; technical intervals keep controls/SDs fixed.','Positive signs/linear-gain ranking do not establish biological benefit or a new mechanism.'),file.path(out,'evaluation_scope.txt'))
cat('Alternative ranking sensitivity evaluation complete',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
