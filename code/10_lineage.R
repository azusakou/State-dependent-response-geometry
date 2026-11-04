#!/usr/bin/env Rscript
# Estimate stress associations and neighbourhood sensitivities.
# Usage: Rscript run_lineage.R [input_directory] [output_directory]
suppressPackageStartupMessages({library(data.table);library(BiocNeighbors)})
args<-commandArgs(trailingOnly=TRUE)
script_dir<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])))
input<-"inputs"
out<-"results/lineage_sensitivity"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache<-readRDS(file.path(input,'lineage_inputs.rds'))
all_inputs<-cache$data;desired_features<-cache$features
stress_axes<-c('hallmark__HYPOXIA','hallmark__REACTIVE_OXYGEN_SPECIES_PATHWAY','hallmark__UNFOLDED_PROTEIN_RESPONSE')
base_cov<-c('cell_cycle_state','stress_state','chromatin_regulator_proxy','metabolic_state','mitochondrial_fraction')
set.seed(20261103)
fit_one<-function(d,y,extended=TRUE,nonlinear=FALSE) {
  dd<-copy(d);dd[,yy:=as.numeric(y)]
  vars<-base_cov
  if(extended) vars<-c(vars,'log_baseline_n','log_descendant_n','log_detected_genes')
  rhs<-paste(vars,collapse='+')
  if(uniqueN(dd$future_condition)>1)rhs<-paste(rhs,'factor(future_condition)',sep='+')
  if(extended && uniqueN(dd$dominant_library)>1)rhs<-paste(rhs,'factor(dominant_library)',sep='+')
  if(nonlinear)rhs<-sub('stress_state','splines::ns(stress_state,df=3)',rhs,fixed=TRUE)
  fit<-lm(as.formula(paste('yy~',rhs)),dd)
  if(nonlinear) return(fit)
  x<-model.matrix(fit); b<-coef(fit); keep<-is.finite(b);x<-x[,keep,drop=FALSE];b<-b[keep]
  bread<-solve(crossprod(x));score<-rowsum(x*residuals(fit),dd$clone_id,reorder=FALSE);g<-nrow(score);n<-nrow(x);p<-ncol(x)
  v<-(g/(g-1))*((n-1)/(n-p))*bread%*%crossprod(score)%*%bread
  se<-sqrt(diag(v));tcrit<-qt(.975,g-1)
  j<-match('stress_state',names(b))
  data.table(beta=b[j],se=se[j],lower=b[j]-tcrit*se[j],upper=b[j]+tcrit*se[j],p=2*pt(-abs(b[j]/se[j]),g-1),units=nrow(dd),clones=uniqueN(dd$clone_id),residual_df=df.residual(fit))
}
neighbours<-function(d,features,npc=15,k=NULL) {
  ans<-vector('list',nrow(d));ks<-integer(nrow(d))
  for(co in unique(d$future_condition)) {
    rr<-which(d$future_condition==co)
    x<-scale(as.matrix(d[rr,..features]));x[!is.finite(x)]<-0
    if(is.finite(npc))x<-prcomp(x,rank.=min(npc,ncol(x),nrow(x)-1),center=FALSE,scale.=FALSE)$x
    kk<-if(is.null(k)) min(60,max(15,floor(sqrt(nrow(x)))),nrow(x)-1) else min(k,nrow(x)-1)
    z<-BiocNeighbors::findKNN(x,k=kk,BNPARAM=KmknnParam())$index
    for(i in seq_along(rr))ans[[rr[i]]]<-rr[z[i,]]
    ks[rr]<-kk
  }
  list(index=ans,k=ks)
}
geometry<-function(y,nn,shrink=0,threshold=.05) {
  rbindlist(lapply(seq_len(nrow(y)),function(i) {
    yy<-y[nn$index[[i]],,drop=FALSE];cv<-cov(yy)
    if(shrink>0)cv<-(1-shrink)*cv+shrink*diag(diag(cv))
    eig<-pmax(eigen(cv,symmetric=TRUE,only.values=TRUE)$values,0); pp<-eig/max(sum(eig),1e-12);loc<-colMeans(yy)
    data.table(rank=if(sum(eig)>0)exp(-sum(pp[pp>0]*log(pp[pp>0])))else 0,entropy=-sum(loc[loc>0]*log(loc[loc>0])),reliable_branches=sum(loc>=threshold))
  }))
}
coef_results<-list();support_results<-list();curves<-list();nonlinearity<-list();rarity<-list();tech<-list()
for(nm in names(all_inputs)) {
  d<-all_inputs[[nm]];feat<-setdiff(intersect(desired_features,names(d)),stress_axes);branches<-grep('^p_branch_',names(d),value=TRUE);yy<-as.matrix(d[,..branches])
  grid<-unique(rbindlist(list(data.table(setting='original_geometry',npc=15,k=NA_integer_,repr='all',shrink=0),data.table(setting=paste0('k',c(15,30,60)),npc=15,k=c(15,30,60),repr='all',shrink=0),data.table(setting=paste0('pc',c(5,10)),npc=c(5,10),k=NA_integer_,repr='all',shrink=0),data.table(setting=c('pathway_only','regulon_only','unreduced_axes'),npc=c(15,15,Inf),k=NA_integer_,repr=c('pathway','regulon','all'),shrink=0),data.table(setting=c('shrink_0.01','shrink_0.05','shrink_0.20'),npc=15,k=NA_integer_,repr='all',shrink=c(.01,.05,.20)))))
  for(j in seq_len(nrow(grid))) {
    gg<-grid[j]; f<-if(gg$repr=='pathway')feat[grepl('^(progeny|hallmark)__',feat)]else if(gg$repr=='regulon')feat[grepl('^collectri__',feat)]else feat
    nn<-neighbours(d,f,npc=gg$npc,k=if(is.na(gg$k))NULL else gg$k);g<-geometry(yy,nn,shrink=gg$shrink)
    for(extended in c(FALSE,TRUE)) {
      z<-fit_one(d,g$rank,extended);z[,`:=`(dataset=nm,setting=gg$setting,adjustment=if(extended)'states_counts_detection_library'else'original_states',npc=gg$npc,k_min=min(nn$k),k_max=max(nn$k),features=length(f),shrink=gg$shrink,median_rank=median(g$rank))]
      coef_results[[length(coef_results)+1]]<-z
    }
    if(gg$setting=='original_geometry') {
      maxdelta<-max(abs(g$rank-d$local_effective_rank_stress_state_loo));cat(nm,'original agreement',maxdelta,'\n');stopifnot(maxdelta<1e-9)
      for(metric in c('entropy','reliable_branches')) {
        z<-fit_one(d,g[[metric]],TRUE);z[,`:=`(dataset=nm,metric=metric,threshold=.05)];support_results[[length(support_results)+1]]<-z
      }
      for(t in c(.01,.10)) {
        gt<-geometry(yy,nn,threshold=t);z<-fit_one(d,gt$reliable_branches,TRUE);z[,`:=`(dataset=nm,metric='reliable_branches',threshold=t)];support_results[[length(support_results)+1]]<-z
      }
      linear<-fit_one(d,g$rank,TRUE,TRUE)
      # Partial effect contrasts with cluster sandwich uncertainty, holding
      # the empirical distribution of nuisance predictors fixed (g-computation).
      x<-model.matrix(linear);b<-coef(linear);keep<-is.finite(b);x<-x[,keep,drop=FALSE];b<-b[keep];bread<-solve(crossprod(x));ss<-rowsum(x*residuals(linear),d$clone_id,reorder=FALSE);ng<-nrow(ss);v<-(ng/(ng-1))*((nrow(x)-1)/(nrow(x)-ncol(x)))*bread%*%crossprod(ss)%*%bread
      ref<-copy(d);ref[,stress_state:=median(d$stress_state)];xr<-colMeans(model.matrix(delete.response(terms(linear)),ref)[,keep,drop=FALSE]);tt<-delete.response(terms(linear))
      cc<-rbindlist(lapply(seq(quantile(d$stress_state,.05),quantile(d$stress_state,.95),length.out=31),function(s){nd<-copy(d);nd[,stress_state:=s];a<-colMeans(model.matrix(tt,nd)[,keep,drop=FALSE])-xr;est<-sum(a*b);se<-sqrt(as.numeric(t(a)%*%v%*%a));data.table(stress=s,adjusted_rank_contrast=est,lower=est-qt(.975,ng-1)*se,upper=est+qt(.975,ng-1)*se)}));cc[,dataset:=nm];curves[[nm]]<-cc
      # Compare AIC descriptively only: neighbouring outcomes share information.
      dd<-copy(d);dd[,yy:=g$rank];rhs<-paste(c(base_cov,'log_baseline_n','log_descendant_n','log_detected_genes'),collapse='+');if(uniqueN(d$future_condition)>1)rhs<-paste(rhs,'factor(future_condition)',sep='+');rhs<-paste(rhs,'factor(dominant_library)',sep='+');linfit<-lm(as.formula(paste('yy~',rhs)),dd)
      nonlinearity[[nm]]<-data.table(dataset=nm,linear_aic=AIC(linfit),spline_aic=AIC(linear),difference_spline_minus_linear=AIC(linear)-AIC(linfit),contrast_95_vs_5=cc$adjusted_rank_contrast[31]-cc$adjusted_rank_contrast[1])
      # Available confound correlations; no causal interpretation.
      tech[[nm]]<-rbindlist(lapply(c('log_baseline_n','log_descendant_n','log_detected_genes','mitochondrial_fraction'),function(vv)data.table(dataset=nm,covariate=vv,spearman=cor(d$stress_state,d[[vv]],method='spearman'))))
    }
  }
  for(m in c(2,5,10)) {
    ds<-d[descendant_cells>=m];y<-as.matrix(ds[,..branches]);cnt<-round(y*ds$descendant_cells)
    nn<-neighbours(ds,feat)
    # Same retained clone set before and after rarefaction distinguishes selection.
    pre<-fit_one(ds,geometry(y,nn)$rank,TRUE);pre[,`:=`(dataset=nm,m=m,replicate=0L,method='retained_before_downsampling')];rarity[[length(rarity)+1]]<-pre
    set.seed(20261103+match(nm,names(all_inputs))*100+m)
    for(bi in seq_len(100)) {
      subcounts<-t(vapply(seq_len(nrow(cnt)),function(i){c<-cnt[i,];draw<-integer(length(c));left<-m;population<-sum(c);for(j in seq_len(length(c)-1)){draw[j]<-rhyper(1,c[j],population-c[j],left);left<-left-draw[j];population<-population-c[j]};draw[length(c)]<-left;draw},numeric(ncol(cnt))))
      z<-fit_one(ds,geometry(subcounts/m,nn)$rank,TRUE);z[,`:=`(dataset=nm,m=m,replicate=bi,method='without_replacement_equal_descendants')];rarity[[length(rarity)+1]]<-z
    }
    cat(nm,'rarefaction m',m,'done, retained',nrow(ds),'\n')
  }
}
fwrite(rbindlist(coef_results,fill=TRUE),file.path(out,'lineage_sensitivity.tsv'),sep='\t')
fwrite(rbindlist(support_results),file.path(out,'lineage_alternative_metrics.tsv'),sep='\t')
fwrite(rbindlist(curves),file.path(out,'lineage_adjusted_nonlinear_effect.tsv'),sep='\t')
fwrite(rbindlist(nonlinearity),file.path(out,'lineage_nonlinearity_diagnostic.tsv'),sep='\t')
fwrite(rbindlist(rarity),file.path(out,'lineage_rarefaction.tsv'),sep='\t')
fwrite(rbindlist(tech),file.path(out,'lineage_stress_quality_correlations.tsv'),sep='\t')
rar_dt<-rbindlist(rarity)
rar_summary<-rar_dt[replicate>0,.(median_beta=median(beta),minimum_beta=min(beta),maximum_beta=max(beta),negative_runs=sum(beta<0),median_lower=median(lower),median_upper=median(upper)),by=.(dataset,m)]
fwrite(rar_summary,file.path(out,'lineage_rarefaction_summary.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('All lineage sensitivity analyses completed.\n')

# Additional individual clone metrics and full-pipeline cluster bootstrap.
a<-list()
for(nm in names(all_inputs)){
 d<-all_inputs[[nm]];branches<-grep('^p_branch_',names(d),value=TRUE)
 d[,observed_branches:=rowSums(as.matrix(.SD)>0),.SDcols=branches]
 for(metric in c('transition_entropy','observed_branches')){
  r<-fit_one(d,d[[metric]],TRUE);r[,`:=`(dataset=nm,metric=metric)];a[[length(a)+1]]<-r
 }
}
fwrite(rbindlist(a),file.path(out,'lineage_individual_clone_metrics.tsv'),sep='\t')

answer<-list();curve_boot<-list()
for(nm in names(all_inputs)) {
 d<-all_inputs[[nm]];ids<-unique(d$clone_id);feat<-setdiff(intersect(desired_features,names(d)),stress_axes);branches<-grep('^p_branch_',names(d),value=TRUE)
 set.seed(202611031+match(nm,names(all_inputs)))
 for(b in seq_len(200)) {
   map<-data.table(clone_id=sample(ids,length(ids),replace=TRUE),bootstrap_clone_id=seq_along(ids))
   z<-merge(map,d,by='clone_id',allow.cartesian=TRUE);z[,clone_id:=bootstrap_clone_id]
   nn<-neighbours(z,feat);y<-as.matrix(z[,..branches]);g<-geometry(y,nn)
   a<-fit_one(z,g$rank,TRUE);a[,`:=`(dataset=nm,replicate=b)]
   spfit<-fit_one(z,g$rank,TRUE,TRUE)
   grid_stress<-seq(quantile(d$stress_state,.05),quantile(d$stress_state,.95),length.out=31)
   nd<-z[rep(1L,32L)];nd[,stress_state:=c(grid_stress,median(d$stress_state))]
   effect<-predict(spfit,newdata=nd);effect<-effect[1:31]-effect[32]
   curve_boot[[length(curve_boot)+1]]<-data.table(dataset=nm,replicate=b,stress=grid_stress,contrast=effect)
   answer[[length(answer)+1]]<-a
   if(b%%25==0)cat(nm,b,'complete\n')
 }
}
a<-rbindlist(answer);fwrite(a,file.path(out,'lineage_full_pipeline_bootstrap.tsv'),sep='\t')
s<-a[,.(replicates=.N,median_beta=median(beta),lower=quantile(beta,.025),upper=quantile(beta,.975),negative_runs=sum(beta<0)),by=dataset]
fwrite(s,file.path(out,'lineage_full_pipeline_bootstrap_summary.tsv'),sep='\t');print(s)
cb<-rbindlist(curve_boot)
fwrite(cb,file.path(out,'lineage_curve_bootstrap.tsv'),sep='\t')
cs<-cb[,.(bootstrap_lower=quantile(contrast,.025),bootstrap_upper=quantile(contrast,.975)),by=.(dataset,stress)]
fwrite(cs,file.path(out,'lineage_curve_bootstrap_summary.tsv'),sep='\t')
