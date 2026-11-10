#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);stopifnot(length(a)>=2L);cache<-readRDS(a[1]);out<-a[2]
cells<-as.data.table(cache$cells);upr<-colnames(cache$raw_upr);stopifnot(length(upr)==91L,nrow(cache$raw_upr)==nrow(cells),nrow(cache$raw_target)==nrow(cells),identical(colnames(cache$raw_target),cache$target_genes))
res<-list()
for(g in cache$target_genes){
 member<-g%in%upr;normalizer_member<-g%in%cache$universe
 target<-if(normalizer_member)cache$raw_target[,g]else rep(0,nrow(cells));stopifnot(all(is.finite(target)),all(target>=0),all(cells$common_library-target>0))
 denom<-cells$common_library-target
 cols<-if(member)setdiff(upr,g)else upr
 # Recompute each individual log term from counts; no aggregate-score rescaling.
 score<-rowMeans(log1p(1e4*cache$raw_upr[,cols,drop=FALSE]/denom))
 for(co in unique(cells$context)){
  tr<-which(cells$gene==g&cells$context==co&!cells$is_control);ct<-which(cells$context==co&cells$is_control)
  if(!length(tr)||length(ct)<2)next
  res[[paste(g,co)]]<-data.table(study=cache$study,context=co,gene=g,treated_cells=length(tr),control_cells=length(ct),target_in_UPR=member,target_in_normalization=normalizer_member,members_remaining=length(cols),treated_mean=mean(score[tr]),control_mean=mean(score[ct]),control_sd=sd(score[ct]),effect_z=(mean(score[tr])-mean(score[ct]))/sd(score[ct]),min_target_excluded_library=min(denom[c(tr,ct)]),max_target_library_fraction=max(target[c(tr,ct)]/cells$common_library[c(tr,ct)]))
 }
}
r<-rbindlist(res);stopifnot(all(is.finite(r$effect_z)));fwrite(r,out,sep='\t');cat(cache$study,'exact denominator-corrected moments',nrow(r),'rows\n')
