#!/usr/bin/env Rscript
source("code/common.R")
suppressPackageStartupMessages(library(data.table))
set.seed(191104)
out<-'results/prediction_checks';dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-read_norman_input()
norm<-function(v)sqrt(sum(v*v,na.rm=TRUE));cos<-function(a,b)sum(a*b,na.rm=TRUE)/max(norm(a)*norm(b),1e-10)
classify<-function(p,mean,red){q<-p.adjust(p,'BH');fifelse(is.finite(q)&q<.05&mean>0,'synergistic',fifelse(is.finite(q)&q<.05&mean<0,'antagonistic',fifelse(red,'redundant','context-dependent')))}
compute<-function(a,feature_scale=FALSE){
 d<-a$delta;dg<-a$delta_gem
 if(feature_scale){ss<-apply(d,1,sd);ss[!is.finite(ss)|ss<1e-8]<-1;d<-d/ss;dg<-sweep(dg,1,ss,'/')}
 singles<-split(which(a$is_single&!is.na(a$single_gene)),a$single_gene[a$is_single&!is.na(a$single_gene)])
 ds<-lapply(singles,function(idx)as.numeric(d[,idx,drop=FALSE]%*%a$n_cells[idx]/sum(a$n_cells[idx])))
 dgs<-lapply(singles,function(idx){m<-matrix(NA_real_,nrow(d),dim(dg)[3]);for(g in seq_len(ncol(m))){w<-a$n_condition_gem[idx,g];m[,g]<-as.numeric(matrix(dg[,idx,g,drop=FALSE],nrow(d),length(idx))%*%w/sum(w))};m})
 rows<-list();projs<-list()
 for(j in which(a$is_pair)){
  genes<-strsplit(a$conditions[j],'_',fixed=TRUE)[[1]];if(length(genes)!=2||!all(genes%in%names(ds)))next
  da<-ds[[genes[1]]];db<-ds[[genes[2]]];pair<-d[,j];add<-da+db;unit<-add/max(norm(add),1e-10);iv<-pair-add
  v<-vapply(seq_len(dim(dg)[3]),function(g)sum((dg[,j,g]-dgs[[genes[1]]][,g]-dgs[[genes[2]]][,g])*unit,na.rm=TRUE),numeric(1));v<-v[is.finite(v)]
  se<-sd(v)/sqrt(length(v));me<-mean(v);p<-2*pt(-abs(me/se),length(v)-1);ci<-me+c(-1,1)*qt(.975,length(v)-1)*se
  red<-max(cos(pair,da),cos(pair,db))>.8&&norm(pair)<=1.1*max(norm(da),norm(db))
  rows[[length(rows)+1L]]<-data.table(combination=a$conditions[j],gene_a=genes[1],gene_b=genes[2],n_cells=a$n_cells[j],projection=sum(iv*unit),mean=me,se=se,p=p,ci_low=ci[1],ci_high=ci[2],redundant_candidate=red)
  projs[[length(projs)+1L]]<-data.table(combination=a$conditions[j],gem_group=seq_along(v),projection=v)
 }
 r<-rbindlist(rows);r[,q:=p.adjust(p,'BH')];r[,interaction_class:=classify(p,mean,redundant_candidate)]
 list(summary=r,projections=rbindlist(projs))
}
primary<-compute(x);scaled<-compute(x,TRUE)
altproj<-fread('inputs/norman_alternative_projections.tsv')[normalization=='log1p_CPM']
altmeta<-fread('inputs/norman_alternative_metadata.tsv')[normalization=='log1p_CPM']
cpm_summary<-altproj[,.(mean=mean(projection),se=sd(projection)/sqrt(.N),p=2*pt(-abs(mean(projection)/(sd(projection)/sqrt(.N))),.N-1),ci_low=mean(projection)-qt(.975,.N-1)*sd(projection)/sqrt(.N),ci_high=mean(projection)+qt(.975,.N-1)*sd(projection)/sqrt(.N)),by=combination]
cpm_summary<-merge(cpm_summary,altmeta[,.(combination,redundant_candidate,projection)],by='combination')
cpm_summary[,q:=p.adjust(p,'BH')];cpm_summary[,interaction_class:=classify(p,mean,redundant_candidate)]
cpm<-list(summary=cpm_summary,projections=altproj[,.(combination,gem_group,projection)])
orig<-fread('results/norman_control/combination_interactions.tsv')
test<-merge(primary$summary,orig,by='combination');stopifnot(max(abs(test$projection-test$interaction_projection))<1e-10,max(abs(test$mean-test$replicate_projection_mean))<1e-10,all(test$interaction_class.x==test$interaction_class.y))
fwrite(primary$projections,file.path(out,'norman_gem_projections.tsv'),sep='\t')
# Exact sign-flip test and GEM bootstrap are technical sensitivity checks only.
ss<-primary$projections[,{
 v<-projection;n<-length(v);signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)))
 null<-as.numeric(signs%*%v)/n;p<-mean(abs(null)>=abs(mean(v))-1e-12)
 bo<-replicate(2000,mean(sample(v,n,replace=TRUE)));ci<-quantile(bo,c(.025,.975))
 # Direction recomputed without the held-out GEM is not available from this vector alone.
 .(signflip_p=p,bootstrap_ci_low=ci[1],bootstrap_ci_high=ci[2],n_gem_groups=n)
},by=combination]
ss[,signflip_q:=p.adjust(signflip_p,'BH')]
ss<-merge(primary$summary,ss,by='combination');ss[,signflip_class:=fifelse(signflip_q<.05&mean>0,'synergistic',fifelse(signflip_q<.05&mean<0,'antagonistic',fifelse(redundant_candidate,'redundant','context-dependent')))]
# Leave-one-GEM-out sign stability, within the same pooled experiment.
loo<-primary$projections[,.(loo_mean_min=min(vapply(seq_along(projection),function(i)mean(projection[-i]),numeric(1))),loo_mean_max=max(vapply(seq_along(projection),function(i)mean(projection[-i]),numeric(1)))),by=combination]
ss<-merge(ss,loo,by='combination')
fwrite(ss,file.path(out,'norman_technical_uncertainty.tsv'),sep='\t')
variants<-rbindlist(list(primary$summary[,normalization:='log1p_CP10k'],cpm$summary[,normalization:='log1p_CPM'],scaled$summary[,normalization:='log1p_CP10k_gene_SD_scaled']),fill=TRUE)
fwrite(variants,file.path(out,'norman_normalization_sensitivity.tsv'),sep='\t')
fwrite(rbindlist(list(primary$projections[,normalization:='log1p_CP10k'],cpm$projections[,normalization:='log1p_CPM'],scaled$projections[,normalization:='log1p_CP10k_gene_SD_scaled'])),file.path(out,'norman_normalization_gem_projections.tsv'),sep='\t')
ct<-variants[,.N,by=.(normalization,interaction_class)];fwrite(ct,file.path(out,'norman_normalization_class_counts.tsv'),sep='\t')
comparison<-merge(variants,primary$summary[,.(combination,original_class=interaction_class)],by='combination')
ag<-comparison[,.(n=.N,same_class=sum(interaction_class==original_class),original_positive_retained=sum(original_class=='synergistic'&interaction_class=='synergistic'),original_positive_count=sum(original_class=='synergistic')),by=normalization]
fwrite(ag,file.path(out,'norman_normalization_agreement.tsv'),sep='\t')
# Threshold sensitivity for conditional t-test evidence; no biological replication is introduced.
th<-rbindlist(lapply(c(.01,.025,.05,.10),function(q){primary$summary[,.(n=.N,positive=sum(p.adjust(p,'BH')<q&mean>0),negative=sum(p.adjust(p,'BH')<q&mean<0),q_threshold=q)]}))
fwrite(th,file.path(out,'norman_fdr_sensitivity.tsv'),sep='\t')
cat('Primary projections and sensitivity statistics recomputed. Alternative CPM statistics use frozen per-GEM projections.\n')

suppressPackageStartupMessages(library(data.table));set.seed(191105)
out<-'results/prediction_checks';a<-read_norman_input();d<-a$delta;dg<-a$delta_gem;norm<-function(x)sqrt(sum(x*x,na.rm=TRUE))
singles<-split(which(a$is_single&!is.na(a$single_gene)),a$single_gene[a$is_single&!is.na(a$single_gene)])
dgs<-lapply(singles,function(idx){m<-matrix(NA_real_,nrow(d),dim(dg)[3]);for(g in seq_len(ncol(m))){w<-a$n_condition_gem[idx,g];m[,g]<-as.numeric(matrix(dg[,idx,g,drop=FALSE],nrow(d),length(idx))%*%w/sum(w))};m})
rows<-list()
for(j in which(a$is_pair)){
 genes<-strsplit(a$conditions[j],'_',fixed=TRUE)[[1]];if(length(genes)!=2||!all(genes%in%names(dgs)))next
 sg<-dgs[[genes[1]]]+dgs[[genes[2]]];proj<-vapply(seq_len(ncol(sg)),function(g){u<-rowMeans(sg[,-g,drop=FALSE],na.rm=TRUE);u<-u/max(norm(u),1e-10);sum((dg[,j,g]-sg[,g])*u,na.rm=TRUE)},numeric(1));v<-proj[is.finite(proj)];me<-mean(v);se<-sd(v)/sqrt(length(v));p<-2*pt(-abs(me/se),length(v)-1)
 for(g in seq_along(proj))rows[[length(rows)+1L]]<-data.table(combination=a$conditions[j],gem_group=g,projection=proj[g],mean=me,se=se,p=p)
}
r<-rbindlist(rows);z<-unique(r[,.(combination,mean,se,p)]);z[,q:=p.adjust(p,'BH')];z[,positive:=q<.05&mean>0]
fwrite(r[,.(combination,gem_group,projection)],file.path(out,'norman_crossfit_gem_projections.tsv'),sep='\t');fwrite(z,file.path(out,'norman_crossfit_direction_sensitivity.tsv'),sep='\t')
orig<-fread(file.path(out,'norman_technical_uncertainty.tsv'));print(merge(orig[interaction_class=='synergistic',.(combination,original_mean=mean,original_q=q)],z,by='combination'));cat('Crossfit positives:',sum(z$positive),'\n')
