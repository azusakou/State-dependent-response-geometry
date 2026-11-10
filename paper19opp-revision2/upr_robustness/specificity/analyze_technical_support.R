#!/usr/bin/env Rscript
# Previously specified technical/depth sensitivities; fixed readouts and panels.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);project<-normalizePath(if(length(a))a[1]else'paper19opp-revision2');stage<-if(length(a)>1)a[2]else'depth';out<-file.path(project,'upr_robustness/specificity');mod<-file.path(project,'upr_robustness/modules');source(file.path(out,'statistics_core.R'));set.seed(2026110613L)
z<-readRDS(file.path(out,'primary_signature_state.rds'));pairs<-z$pairs;axes<-z$definitions$axis;panels<-fread(file.path(out,'fixed_target_panels.tsv'));studies<-c('K562','RPE1','Adamson','Frangieh');comps<-c('RPE1','Adamson','Frangieh_Control','Frangieh_IFNg','Frangieh_Coculture')
if(stage=='depth'){
 B<-200L;values<-list();check<-list()
 for(st in studies){
  ca<-readRDS(file.path(mod,'scored',paste0(st,'_cell_cache.rds')));md<-as.data.table(ca$metadata);allowed<-if(st=='K562')md$gene%in%unique(panels$gene)else paste(md$gene,md$context,sep='|')%in%panels[validation_study==st,paste(gene,validation_context,sep='|')];keep<-which(md$is_control==FALSE&allowed);sc<-ca$real_scores_exact[keep,axes,drop=FALSE];md<-md[keep];rm(ca);gc(FALSE)
  ge<-fread(file.path(mod,'scored',paste0(st,'_gene_effects.tsv.gz')))[mode=='exact_Lminus_target'&axis%in%axes];gm<-fread(file.path(mod,'scored',paste0(st,'_GEM_effects.tsv.gz')))[mode=='exact_Lminus_target'&axis%in%axes];key<-paste(md$gene,md$context,md$gem,sep='|');v<-sc
  for(ax in axes){g<-gm[axis==ax];j<-match(key,paste(g$gene,g$context,g$gem,sep='|'));stopifnot(!anyNA(j));v[,ax]<-(sc[,ax]-g$matched_control_mean[j])/g$control_sd[j]
   ck<-data.table(gene=md$gene,context=md$context,x=v[,ax])[,.(new_delta=mean(x)),by=.(gene,context)];gg<-merge(ck,ge[axis==ax,.(gene,context,delta_z)],by=c('gene','context'));err<-abs(gg$new_delta-gg$delta_z);stopifnot(all(!is.finite(gg$delta_z)|err<1e-8));check[[length(check)+1L]]<-data.table(study=st,axis=ax,max_error=if(any(is.finite(err)))max(err,na.rm=TRUE)else NA_real_)
  }
  values[[st]]<-list(metadata=md,values=v);rm(sc,v,md,ge,gm);gc(FALSE)
 }
 fwrite(rbindlist(check),file.path(out,'depth_cell_residual_reconstruction.tsv'),sep='\t');ans<-list();draws<-list();i<-0L
 for(cp in comps)for(depth in if(cp=='RPE1')100L else c(30L,15L)){
  pp<-panels[comparison==cp];st<-unique(pp$validation_study);ctx<-unique(pp$validation_context);sv<-values$K562;vv<-values[[st]];sx<-vector('list',nrow(pp));vy<-sx
  for(k in seq_len(nrow(pp))){g<-pp$gene[k];si<-which(sv$metadata$gene==g);vi<-which(vv$metadata$gene==g&vv$metadata$context==ctx);sx[[k]]<-if(length(si)>=depth)replicate(B,sample(si,depth,replace=FALSE))else NULL;vy[[k]]<-if(length(vi)>=depth)replicate(B,sample(vi,depth,replace=FALSE))else NULL}
  for(ax in axes){i<-i+1L;xs<-ys<-matrix(NA_real_,nrow(pp),B,dimnames=list(pp$gene,NULL));for(k in seq_len(nrow(pp))){if(!is.null(sx[[k]]))xs[k,]<-colMeans(matrix(sv$values[sx[[k]],ax],depth,B));if(!is.null(vy[[k]]))ys[k,]<-colMeans(matrix(vv$values[vy[[k]],ax],depth,B))}
   dd<-rbindlist(lapply(seq_len(B),function(b){ok<-is.finite(xs[,b])&is.finite(ys[,b]);ag<-sp_agreement(xs[ok,b],ys[ok,b]);data.table(draw=b,rho=sp_rho(xs[ok,b],ys[ok,b]),C=if(any(ok))mean(sign(xs[ok,b])*ys[ok,b])else NA_real_,direction=ag['concordance'],genes=sum(ok))}));dd[,`:=`(comparison=cp,axis=ax,target_depth=depth)];draws[[i]]<-dd
   ans[[i]]<-rbindlist(lapply(c('rho','C','direction'),function(metric){v<-dd[[metric]];f<-v[is.finite(v)];data.table(comparison=cp,axis=ax,target_depth=depth,metric=metric,fixed_target_genes=nrow(pp),min_estimable_genes=min(dd$genes),median=if(length(f))median(f)else NA_real_,resampling_low=if(length(f))quantile(f,.025)else NA_real_,resampling_high=if(length(f))quantile(f,.975)else NA_real_,finite_draws=length(f))}));cat(cp,depth,ax,'depth sensitivity complete\n');flush.console()
  }
 }
 fwrite(rbindlist(ans),file.path(out,'equal_depth_signature_summaries.tsv'),sep='\t');fwrite(rbindlist(draws),file.path(out,'equal_depth_signature_draws.tsv.gz'),sep='\t')
}else if(stage=='wild'){
 B<-2000L;gs<-lapply(studies,function(st)fread(file.path(mod,'scored',paste0(st,'_GEM_effects.tsv.gz')))[mode=='exact_Lminus_target'&axis%in%axes]);names(gs)<-studies
 signs<-lapply(c('K562','RPE1'),function(st){kk<-sort(unique(gs[[st]][genuine_GEM_metadata==TRUE,as.character(gem)]));matrix(sample(c(-1,1),length(kk)*B,replace=TRUE),length(kk),B,dimnames=list(kk,NULL))});names(signs)<-c('K562','RPE1')
 wild<-function(gm,genes,point,sign){
  gm<-gm[gene%in%genes&is.finite(contrast_z)&n_target>0];n<-gm[,.(n=sum(n_target)),by=gene];pp<-gm[,.(point=weighted.mean(contrast_z,n_target)),by=gene];pv<-setNames(pp$point,pp$gene)[genes];stopifnot(max(abs(pv-point[genes]))<1e-8)
  total<-setNames(n$n,n$gene);gm[,residual:=n_target/total[gene]*(contrast_z-point[gene])];R<-matrix(0,length(genes),nrow(sign),dimnames=list(genes,rownames(sign)));R[cbind(match(gm$gene,genes),match(as.character(gm$gem),rownames(sign)))]<-gm$residual;sweep(R%*%sign,1,point[genes],'+')
 }
 ans<-list();i<-0L
 for(cp in comps)for(ax in axes){
  pp<-pairs[comparison==cp&axis==ax&estimable_pair==TRUE];gg<-pp$gene;x<-setNames(pp$source_deltaZ,gg);y<-setNames(pp$validation_deltaZ,gg);xx<-wild(gs$K562[axis==ax],gg,x,signs$K562);yy<-if(cp=='RPE1')wild(gs$RPE1[axis==ax],gg,y,signs$RPE1)else matrix(y,nrow=length(y),ncol=B,dimnames=list(gg,NULL))
  dr<-vapply(seq_len(B),function(b){ag<-sp_agreement(xx[,b],yy[,b]);c(rho=sp_rho(xx[,b],yy[,b]),C=mean(sign(xx[,b])*yy[,b]),direction=ag['concordance'],kappa=ag['kappa'],MCC=ag['MCC'])},numeric(5))
  i<-i+1L;ans[[i]]<-rbindlist(lapply(seq_len(nrow(dr)),function(j){q<-sp_ci(dr[j,]);data.table(comparison=cp,axis=ax,metric=sub('\\..*','',rownames(dr)[j]),technical_low=q[1],technical_high=q[2],finite_draws=sum(is.finite(dr[j,])),uncertainty_scope=if(cp=='RPE1')'K562 and RPE1 shared genuine-GEM weights; controlSD fixed'else'K562 shared-GEM only; external technical uncertainty not represented')}))
 }
 fwrite(rbindlist(ans),file.path(out,'shared_GEM_technical_signature_intervals.tsv'),sep='\t')
}else stop('Use depth or wild stage')
writeLines(capture.output(sessionInfo()),file.path(out,paste0(stage,'_session_info.txt')))
