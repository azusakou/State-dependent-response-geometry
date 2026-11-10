#!/usr/bin/env Rscript
# Implements the predeclared real-UPR denominator sensitivity; no re-ranking.
a<-commandArgs(TRUE);stopifnot(length(a)>=2L);mode<-a[1];o<-a[2]
suppressPackageStartupMessages({library(data.table);library(Matrix)})
universe<-fread(file.path(o,'phase2_gene_universe.tsv'))$gene;upr<-fread(file.path(o,'phase2_UPR_members.tsv'))$gene;genes<-fread(file.path(o,'phase2_eligible_genes.tsv'))$gene
stopifnot(length(universe)==6780L,length(upr)==90L,length(genes)==438L)
if(mode=='source'){
 ca<-readRDS(file.path(o,'original91_source_raw_cache.rds'));cells<-as.data.table(ca$cells);cells[,`:=`(lib=phase2_library,own_count=target_raw_count)];raw<-ca$raw_upr;ctcells<-ca$control_cells;ctraw<-ca$raw_control_targets;label<-'K562'
}else if(mode=='rpe1'){
 ca<-readRDS(file.path(o,'RPE1_phase2_cache.rds'));cells<-as.data.table(ca$cell_metadata);cells[,`:=`(lib=common_library,own_count=assigned_target_count)];raw<-ca$cell_UPR_raw_counts;ctcells<-ca$control_cells;ctraw<-ca$control_target_counts;label<-'RPE1'
}else stop('Unknown mode')
stopifnot(nrow(raw)==nrow(cells),all(upr%in%colnames(raw)))
gene_out<-list();gem_out<-list()
for(g in genes){
 ii<-which(cells$is_control|cells$gene==g);z<-cells[ii];ct<-which(z$is_control);tr<-which(!z$is_control&z$gene==g);stopifnot(length(tr)>0,length(ct)>1)
 in_norm<-g%in%universe;member<-g%in%upr;own<-numeric(nrow(z))
 if(in_norm){j<-match(g,colnames(ctraw));stopifnot(!is.na(j));cr<-match(z$cell_id[ct],ctcells$cell_id);stopifnot(!anyNA(cr));own[ct]<-ctraw[cr,j];own[tr]<-z$own_count[tr]}
 denom<-z$lib-own;stopifnot(all(is.finite(own)),all(own>=0),all(denom>0))
 members<-if(member)setdiff(upr,g)else upr
 counts<-as.matrix(raw[ii,members,drop=FALSE]);score<-rowMeans(log1p(1e4*counts/denom));sigma<-sd(score[ct]);stopifnot(is.finite(sigma),sigma>0)
 tm<-data.table(gem_group=z$gem_group[tr],score=score[tr])[,.(n_target=.N,treated_mean=mean(score)),by=gem_group]
 cm<-data.table(gem_group=z$gem_group[ct],score=score[ct])[,.(n_control=.N,control_mean=mean(score)),by=gem_group]
 mm<-merge(tm,cm,by='gem_group',all.x=TRUE);stopifnot(!anyNA(mm$control_mean));mm[,`:=`(gene=g,axis='UPR',context=label,control_sd=sigma,target_in_normalization=in_norm,target_in_UPR=member,members_remaining=length(members))];mm[,`:=`(delta_raw=treated_mean-control_mean,deltaZ=(treated_mean-control_mean)/sigma)]
 d<-weighted.mean(mm$delta_raw,mm$n_target);f<-mean(mm[n_target>=5]$delta_raw>0)
 gene_out[[g]]<-data.table(gene=g,axis='UPR',context=label,delta_raw=d,control_sd=sigma,deltaZ=d/sigma,f_positive_GEM=f,source_cells=sum(mm$n_target),source_GEM_labels=nrow(mm),source_GEM_ge5=sum(mm$n_target>=5),target_in_normalization=in_norm,target_in_UPR=member,members_remaining=length(members),control_cells=length(ct))
 gem_out[[g]]<-mm
}
fwrite(rbindlist(gene_out),file.path(o,paste0(mode,'_UPR_exact_denominator_gene_effects.tsv')),sep='\t');fwrite(rbindlist(gem_out),file.path(o,paste0(mode,'_UPR_exact_denominator_GEM_effects.tsv.gz')),sep='\t');cat(mode,'exact realUPR sensitivity complete;438 locked genes; no ranking changed.\n')
