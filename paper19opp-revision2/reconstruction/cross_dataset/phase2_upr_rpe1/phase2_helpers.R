suppressPackageStartupMessages({library(data.table);library(Matrix);library(HDF5Array);library(hdf5r)})
read_annot<-function(f,group,name){
 node<-f[[group]][[name]]
 if(inherits(node,'H5Group')){z<-node[['codes']][];lev<-node[['categories']][];ans<-rep(NA_character_,length(z));ok<-z>=0;ans[ok]<-lev[z[ok]+1L];return(ans)}
 ans<-node[]
 if('__categories'%in%names(f[[group]])&&name%in%names(f[[group]][['__categories']]))ans<-f[[group]][['__categories']][[name]][as.integer(ans)+1L]
 ans
}
make_coefficients<-function(universe,spec){
 pos<-matrix(0,length(universe),length(spec),dimnames=list(universe,names(spec)));neg<-pos
 for(a in seq_along(spec)){
  p<-spec[[a]]$positive;n<-spec[[a]]$negative
  stopifnot(length(p)>1L,all(c(p,n)%in%universe))
  pos[match(p,universe),a]<-1/length(p)
  if(length(n))neg[match(n,universe),a]<-1/length(n)
 }
 list(positive=pos,negative=neg)
}
score_log_counts<-function(logcounts,universe,spec,coeff,assigned_gene){
 pos<-as.matrix(t(logcounts)%*%coeff$positive);neg<-as.matrix(t(logcounts)%*%coeff$negative)
 score<-pos-neg;target_index<-match(assigned_gene,universe);targetlog<-numeric(ncol(logcounts));ok<-!is.na(target_index)
 targetlog[ok]<-logcounts[cbind(target_index[ok],which(ok))]
 for(a in seq_along(spec)){
  p<-spec[[a]]$positive;n<-spec[[a]]$negative
  ip<-which(assigned_gene%in%p);im<-which(assigned_gene%in%n)
  if(length(ip))score[ip,a]<-score[ip,a]+(pos[ip,a]-targetlog[ip])/(length(p)-1L)
  if(length(im))score[im,a]<-score[im,a]-if(length(n)>1L)(neg[im,a]-targetlog[im])/(length(n)-1L)else-neg[im,a]
 }
 list(scores=score,positive=pos,negative=neg,target_log=targetlog)
}
control_base_moments<-function(base_scores,gem){
 gems<-sort(unique(gem));ng<-tabulate(match(gem,gems),length(gems))
 mu<-rowsum(base_scores,match(gem,gems),reorder=TRUE)/ng
 list(gems=gems,n=ng,mean=mu,sd=apply(base_scores,2,sd))
}
# Exact numerator correction under a FIXED normalization denominator.
# The nonlinear leave-target-out denominator is a separately declared UPR check.
target_control_moments<-function(gene,spec,positive,negative,target_log,gem,base){
 mu<-base$mean;sig<-base$sd
 affected<-which(vapply(spec,function(a)gene%in%c(a$positive,a$negative),logical(1)))
 if(length(affected)){
  j<-match(gene,colnames(target_log));stopifnot(!is.na(j),all(is.finite(target_log[,j])))
  lg<-target_log[,j]
  for(a in affected){
   p<-spec[[a]]$positive;n<-spec[[a]]$negative;q<-positive[,a]-negative[,a]
   if(gene%in%p)q<-q+(positive[,a]-lg)/(length(p)-1L)
   if(gene%in%n)q<-q-if(length(n)>1L)(negative[,a]-lg)/(length(n)-1L)else-negative[,a]
   sig[a]<-sd(q);mu[,a]<-as.numeric(rowsum(q,match(gem,base$gems),reorder=TRUE))/base$n
  }
 }
 list(mean=mu,sd=sig)
}
gene_effects_from_cache<-function(cache,spec){
 base<-control_base_moments(cache$control_positive-cache$control_negative,cache$control_cells$gem_group)
 units<-cache$unit_metadata;counts<-units$n_cells;ans<-list();gans<-list()
 genes<-cache$analysis_genes
 for(g in genes){
  ii<-which(units$gene==g);stopifnot(length(ii)>0L)
  cm<-target_control_moments(g,spec,cache$control_positive,cache$control_negative,cache$control_target_log,cache$control_cells$gem_group,base)
  gi<-match(units$gem_group[ii],base$gems);stopifnot(!anyNA(gi))
  delta<-cache$unit_score_sums[ii,,drop=FALSE]/counts[ii]-cm$mean[gi,,drop=FALSE]
  for(a in seq_along(spec)){
   raw<-sum(delta[,a]*counts[ii])/sum(counts[ii]);sigma<-cm$sd[a]
   f<-if(any(counts[ii]>=5))mean(delta[counts[ii]>=5,a]>0)else NA_real_
   ans[[length(ans)+1L]]<-data.table(gene=g,axis=names(spec)[a],delta_raw=raw,control_sd=sigma,deltaZ=raw/sigma,f_positive_GEM=f,source_cells=sum(counts[ii]),source_GEM_labels=length(ii),source_GEM_ge5=sum(counts[ii]>=5))
   gans[[length(gans)+1L]]<-data.table(gene=g,gem_group=units$gem_group[ii],axis=names(spec)[a],n_target=counts[ii],control_mean=cm$mean[gi,a],treated_mean=cache$unit_score_sums[ii,a]/counts[ii],delta_raw=delta[,a],deltaZ=delta[,a]/sigma,control_sd=sigma)
  }
 }
 list(gene=rbindlist(ans),gem=rbindlist(gans))
}
