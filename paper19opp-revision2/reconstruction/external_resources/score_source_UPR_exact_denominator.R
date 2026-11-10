#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Matrix)})
r<-'paper19opp-revision2/reconstruction';cachefile<-file.path(r,'cross_dataset/phase2_upr_rpe1/original91_source_raw_cache.rds');stopifnot(file.exists(cachefile));ca<-readRDS(cachefile);cells<-as.data.table(ca$cells)
upr<-colnames(ca$raw_upr);stopifnot(length(upr)==91L,nrow(ca$raw_upr)==nrow(cells),identical(rownames(ca$raw_upr),cells$cell_id))
ctrow<-match(cells$cell_id[which(cells$is_control)],ca$control_cells$cell_id);stopifnot(!anyNA(ctrow));records<-list()
for(g in ca$fixed_external_genes){
 use<-which(cells$is_control|cells$gene==g);z<-cells[use];ct<-which(z$is_control);tr<-which(!z$is_control&z$gene==g);stopifnot(length(tr)>0,length(ct)>1)
 member<-g%in%upr;normalizer_member<-g%in%ca$universe
 target<-numeric(length(use));if(normalizer_member){j<-match(g,colnames(ca$raw_control_targets));stopifnot(!is.na(j));target[ct]<-ca$raw_control_targets[match(z$cell_id[ct],ca$control_cells$cell_id),j];target[tr]<-z$target_raw_count[tr]}
 denominator<-z$common_library-target;stopifnot(all(is.finite(target)),all(target>=0),all(denominator>0))
 selected<-if(member)setdiff(upr,g)else upr;raw<-as.matrix(ca$raw_upr[use,selected,drop=FALSE]);score<-rowMeans(log1p(1e4*raw/denominator))
 records[[g]]<-data.table(study='Replogle2022',context='K562',gene=g,treated_cells=length(tr),control_cells=length(ct),target_in_UPR=member,target_in_normalization=normalizer_member,members_remaining=length(selected),treated_mean=mean(score[tr]),control_mean=mean(score[ct]),control_sd=sd(score[ct]),effect_z=(mean(score[tr])-mean(score[ct]))/sd(score[ct]),min_target_excluded_library=min(denominator),max_target_library_fraction=max(target/z$common_library))
}
x<-rbindlist(records);stopifnot(all(is.finite(x$effect_z)));fwrite(x,file.path(r,'external_resources/replogle_UPR_exact_target_excluded_moments.tsv'),sep='\t');cat('Exact original91/7848 source composition moments:',nrow(x),'fixed-cohort-union genes; no RPE1 expression accessed.\n')
