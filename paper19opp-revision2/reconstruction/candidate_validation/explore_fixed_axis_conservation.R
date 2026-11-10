#!/usr/bin/env Rscript
# Exploratory follow-up of all fixed axes after the candidate-priority test.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);root<-if(length(a))normalizePath(a[1])else normalizePath('.');out<-file.path(root,'paper19opp-revision2/reconstruction/candidate_validation');cross<-file.path(root,'paper19opp-revision2/reconstruction/cross_dataset');set.seed(202611061L);B<-2000L;P<-10000L
membership<-fread(file.path(out,'final_axis_membership.tsv'));axes<-unique(membership$program);cols<-paste0('score_',axes);stopifnot(length(axes)==12)
rank<-fread(file.path(out,'harmonized_discovery_ranking.tsv'));x<-as.data.table(readRDS(file.path(cross,'harmonized_score_part_1.rds')));x[,gene:=as.character(gene)];ctrl<-x[gene=='non-targeting'];mu<-vapply(cols,function(v)mean(ctrl[[v]]),numeric(1));sdv<-vapply(cols,function(v)sd(ctrl[[v]]),numeric(1));source<-x[gene!='non-targeting',c(list(discovery_cells=.N),lapply(.SD,mean)),by=gene,.SDcols=cols];source<-source[gene%in%rank[metadata_eligible==TRUE,gene]];D<-sweep(sweep(as.matrix(source[,..cols]),2,mu,'-'),2,sdv,'/');source[,c(paste0('effect_',axes)):=as.data.table(D)];fwrite(source[,c('gene','discovery_cells',paste0('effect_',axes)),with=FALSE],file.path(out,'discovery_12axis_effects.tsv'),sep='\t');generic<-colMeans(D)
q2<-function(v)as.numeric(quantile(v[is.finite(v)],c(.025,.975)))
corrows<-list();vectorrows<-list();effectrows<-list();generows<-list()
run<-function(cells,study,minimum,candidate_file){
 gset<-fread(file.path(out,candidate_file))$gene
 for(ctx in unique(cells$context)){
  z<-cells[context==ctx];ct<-z[is_control==TRUE];cmean<-vapply(cols,function(v)mean(ct[[v]]),numeric(1));csd<-vapply(cols,function(v)sd(ct[[v]]),numeric(1));stopifnot(all(csd>0))
  ym<-z[is_control==FALSE&gene%in%gset,c(list(validation_cells=.N),lapply(.SD,mean)),by=gene,.SDcols=cols];ym<-ym[validation_cells>=minimum];mm<-merge(source,ym,by='gene',suffixes=c('_source','_external'));n<-nrow(mm);stopifnot(n>=10)
  X<-as.matrix(mm[,paste0('effect_',axes),with=FALSE]);Y<-sweep(sweep(as.matrix(mm[,paste0(cols,'_external'),with=FALSE]),2,cmean,'-'),2,csd,'/')
  ids<-replicate(B,sample.int(n,n,replace=TRUE),simplify=FALSE);perms<-replicate(P,sample.int(n,n,replace=FALSE),simplify=FALSE)
  for(j in seq_along(axes)){
   rho<-cor(X[,j],Y[,j],method='spearman');bb<-vapply(ids,function(ii)suppressWarnings(cor(X[ii,j],Y[ii,j],method='spearman')),numeric(1));qq<-q2(bb);rr<-rank(X[,j]);tt<-rank(Y[,j]);pn<-vapply(perms,function(ii)cor(rr,tt[ii]),numeric(1));pv<-(1+sum(abs(pn)>=abs(rho)-1e-12))/(1+length(pn));corrows[[length(corrows)+1L]]<<-data.table(study=study,context=ctx,axis=axes[j],genes=n,spearman=rho,conditional_gene_bootstrap_low=qq[1],conditional_gene_bootstrap_high=qq[2],two_sided_gene_permutation_p=pv)
  }
  for(gr in c('three_programs','nine_hallmarks','all_twelve')){
   jj<-switch(gr,three_programs=1:3,nine_hallmarks=4:12,all_twelve=1:12);xx<-X[,jj,drop=FALSE];yy<-Y[,jj,drop=FALSE];gg<-matrix(generic[jj],nrow=n,ncol=length(jj),byrow=TRUE)
   cs<-rowSums(xx*yy)/sqrt(rowSums(xx^2)*rowSums(yy^2));cg<-rowSums(gg*yy)/sqrt(rowSums(gg^2)*rowSums(yy^2));cs[!is.finite(cs)]<-NA_real_;cg[!is.finite(cg)]<-NA_real_
   eidentity<-rowMeans((xx-yy)^2);ezero<-rowMeans(yy^2);egeneric<-rowMeans((gg-yy)^2)
   for(ref in c('zero','generic_mean')){
    er<-if(ref=='zero')ezero else egeneric;delta<-sqrt(mean(eidentity))-sqrt(mean(er));vb<-vapply(ids,function(ii)sqrt(mean(eidentity[ii]))-sqrt(mean(er[ii])),numeric(1));qq<-q2(vb);vectorrows[[length(vectorrows)+1L]]<<-data.table(study=study,context=ctx,axis_group=gr,genes=n,reference=ref,identity_rmse=sqrt(mean(eidentity)),reference_rmse=sqrt(mean(er)),identity_minus_reference_rmse=delta,conditional_gene_bootstrap_low=qq[1],conditional_gene_bootstrap_high=qq[2],mean_identity_cosine=mean(cs,na.rm=TRUE),median_identity_cosine=median(cs,na.rm=TRUE),mean_generic_cosine=mean(cg,na.rm=TRUE),mean_cosine_advantage=mean(cs-cg,na.rm=TRUE))
   }
   generows[[length(generows)+1L]]<<-data.table(study=study,context=ctx,gene=mm$gene,axis_group=gr,identity_cosine=cs,generic_cosine=cg,identity_MSE=eidentity,zero_MSE=ezero,generic_MSE=egeneric,discovery_cells=mm$discovery_cells,validation_cells=mm$validation_cells)
  }
  ef<-rbindlist(lapply(seq_along(axes),function(j)data.table(study=study,context=ctx,gene=mm$gene,axis=axes[j],discovery_effect=X[,j],external_effect=Y[,j],discovery_cells=mm$discovery_cells,validation_cells=mm$validation_cells)));effectrows[[length(effectrows)+1L]]<<-ef
 }
}
run(fread(file.path(out,'adamson_cell_scores_harmonized.tsv.gz')),'Adamson2016',100L,'harmonized_adamson_candidates.tsv')
run(fread(file.path(out,'frangieh_cell_scores_harmonized.tsv.gz')),'Frangieh2021',30L,'harmonized_frangieh_candidates.tsv')
c<-rbindlist(corrows);stopifnot(nrow(c)==48);c[,BH_q_all_48:=p.adjust(two_sided_gene_permutation_p,'BH')];fwrite(c,file.path(out,'exploratory_axis_conservation.tsv'),sep='\t');fwrite(rbindlist(vectorrows),file.path(out,'exploratory_vector_prediction.tsv'),sep='\t');fwrite(rbindlist(effectrows),file.path(out,'exploratory_all_gene_axis_effects.tsv.gz'),sep='\t');fwrite(rbindlist(generows),file.path(out,'exploratory_gene_vector_metrics.tsv.gz'),sep='\t');fwrite(data.table(axis=axes,discovery_generic_mean_effect=generic,discovery_genes=nrow(source)),file.path(out,'discovery_generic_mean_12axis.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(out,'axis_followup_sessionInfo.txt'));print(c[order(BH_q_all_48)])
