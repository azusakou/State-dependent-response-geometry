#!/usr/bin/env Rscript
# Own-target exclusion is applied identically to treated and matched controls.
suppressPackageStartupMessages({library(data.table);library(hdf5r)})
root<-'paper19opp-revision2/reconstruction';out<-file.path(root,'external_resources');started<-Sys.time()
axis<-'score_hallmark_UNFOLDED_PROTEIN_RESPONSE';upr<-fread(file.path(root,'cross_dataset/final_axis_membership.tsv'))[program=='hallmark_UNFOLDED_PROTEIN_RESPONSE'&in_common==TRUE,gene];stopifnot(length(upr)==91L)
annot<-function(f,g,n){o<-f[[g]][[n]];if(inherits(o,'H5Group')){v<-o[['codes']][];lev<-o[['categories']][];return(lev[v+1L])};v<-o[];if('__categories'%in%f[[g]]$ls()$name&&n%in%f[[g]][['__categories']]$ls()$name)v<-f[[g]][['__categories']][[n]][][v+1L];v}
records<-list();gemrecords<-list();checks<-list()
process<-function(d,study,gset,own,libfield,make_gem=FALSE){
 d<-copy(d);d[,S:=get(axis)];for(co in unique(d$context)){
  z<-d[context==co];origrows<-match(z$cell_id,d$cell_id);ctl<-which(z$is_control==TRUE);stopifnot(length(ctl)>1)
  origmu<-mean(z$S[ctl]);origsd<-sd(z$S[ctl]);base<-z[is_control==FALSE&gene%in%gset,.(treated_cells=.N,treated_mean=mean(S)),by=gene]
  base[,`:=`(study=study,context=co,target_in_upr=gene%in%upr,members_remaining=91L,control_cells=length(ctl),control_mean=origmu,control_sd=origsd,original_treated_mean=treated_mean,original_control_mean=origmu,original_control_sd=origsd,treated_own_logcp10k=NA_real_,control_own_logcp10k=NA_real_)]
  for(g in intersect(base$gene,upr)){
   stopifnot(g%in%names(own));v<-own[[g]][origrows];adj<-(91*z$S-v)/90;tr<-which(z$gene==g&!z$is_control)
   base[gene==g,`:=`(members_remaining=90L,treated_mean=mean(adj[tr]),control_mean=mean(adj[ctl]),control_sd=sd(adj[ctl]),treated_own_logcp10k=mean(v[tr]),control_own_logcp10k=mean(v[ctl]))]
  }
  base[,`:=`(effect_z=(treated_mean-control_mean)/control_sd,original_effect_z=(original_treated_mean-original_control_mean)/original_control_sd)]
  records[[paste(study,co)]]<<-base
  if(make_gem){
   cc<-z[is_control==TRUE,.(control_cells=.N,control_mean=mean(S),control_sd=sd(S)),by=.(gem=gem_group)]
   tt<-z[is_control==FALSE&gene%in%base$gene,.(treated_cells=.N,treated_mean=mean(S)),by=.(gene,gem=gem_group)]
   gg<-merge(tt,cc,by='gem');gg<-merge(gg,base[,.(gene,global_control_sd=control_sd,members_remaining)],by='gene')
   gg[,effect_z:=(treated_mean-control_mean)/global_control_sd]
   for(g in intersect(base$gene,upr)){
    vv<-own[[g]][origrows];adj<-(91*z$S-vv)/90
    tr<-data.table(gene=z$gene,gem=z$gem_group,ctl=z$is_control,score=adj)
    cm<-tr[ctl==TRUE,.(control_cells=.N,control_mean=mean(score),control_sd=sd(score)),by=gem]
    tm<-tr[ctl==FALSE&gene==g,.(treated_cells=.N,treated_mean=mean(score)),by=.(gene,gem)]
    qq<-merge(tm,cm,by='gem');qq[,`:=`(members_remaining=90L,global_control_sd=base[gene==g,control_sd])];qq[,effect_z:=(treated_mean-control_mean)/global_control_sd]
    gg<-rbindlist(list(gg[gene!=g],qq),use.names=TRUE,fill=TRUE)
   }
   gemrecords[[study]]<<-gg
  }
 }
 saveRDS(list(cell_id=d$cell_id,own_logcp10k=own),file.path(out,paste0(tolower(study),'_UPR_own_expression.rds')),compress='gzip')
}
# Replogle: dense count rows; all normalization values come from the frozen scored cells.
x<-as.data.table(readRDS(file.path(root,'cross_dataset/harmonized_score_part_1.rds')));x[,`:=`(gene=as.character(gene),context='K562')];gset<-fread(file.path(root,'candidate_validation/harmonized_discovery_ranking.tsv'))[metadata_eligible==TRUE,gene]
f<-H5File$new('/Volumes/wang/paper6yhn/raw/replogle_k562_essential/K562_essential_raw_singlecell_01.h5ad','r');gen<-as.character(annot(f,'var','gene_name'));ri<-match(x$cell_id,as.character(annot(f,'obs','cell_barcode')));stopifnot(!anyNA(ri))
needed<-intersect(gset,upr);own<-list();rawtarget<-matrix(0,length(needed),length(annot(f,'obs','cell_barcode')));rowix<-lapply(needed,function(g)which(gen==g));stopifnot(all(lengths(rowix)>0))
# Stream contiguous full-gene cell blocks once; avoid repeated strided whole-file reads.
for(st in seq.int(1,ncol(rawtarget),by=5000L)){en<-min(st+4999L,ncol(rawtarget));block<-f[['X']][,st:en];for(j in seq_along(needed))rawtarget[j,st:en]<-colSums(block[rowix[[j]],,drop=FALSE]);if(en%%50000L==0L){cat('Replogle raw block through',en,'elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'seconds\n');flush.console()}}
for(j in seq_along(needed))own[[needed[j]]]<-log1p(1e4*rawtarget[j,ri]/x$common_library);rm(rawtarget,block);f$close_all();process(x,'Replogle2022',gset,own,'common_library',TRUE);cat('Replogle own exclusions',length(own),'elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'seconds\n');flush.console()
# Adamson: CSR raw counts, with the exact 63(mod) control and frozen common library.
x<-fread(file.path(root,'candidate_validation/adamson_cell_scores_harmonized.tsv.gz'));gset<-fread(file.path(root,'candidate_validation/harmonized_adamson_candidates.tsv'))$gene;needed<-intersect(gset,upr)
f<-H5File$new('/Volumes/wang/paper4rfv/data/processed/h5ad/core_unpacked/adamson_2016_upr__GSM2406681_10X010.h5ad','r');gen<-as.character(annot(f,'var','_index'));ids<-as.character(annot(f,'obs','_index'));ri<-match(x$cell_id,ids);stopifnot(!anyNA(ri),identical(ri,seq_along(ids)));ptr<-f[['X']][['indptr']][];wanted<-match(needed,gen);stopifnot(!anyNA(wanted));ownmat<-matrix(0,nrow(x),length(needed));colnames(ownmat)<-needed
for(st in seq.int(1,nrow(x),by=2000L)){en<-min(st+1999L,nrow(x));pp<-ptr[st:(en+1L)];ix<-seq.int(pp[1]+1,tail(pp,1));jj<-f[['X']][['indices']][ix]+1L;vv<-f[['X']][['data']][ix];rr<-rep(seq.int(st,en),diff(pp));cc<-match(jj,wanted);ok<-!is.na(cc);ownmat[cbind(rr[ok],cc[ok])]<-vv[ok]}
f$close_all();own<-lapply(seq_along(needed),function(j)log1p(1e4*ownmat[,j]/x$library_size));names(own)<-needed
# Treated own-expression reconstruction also validates orientation/normalization.
for(g in needed){ii<-which(x$gene==g);stopifnot(max(abs(own[[g]][ii]-x$target_logcp10k[ii]))<1e-9)}
process(x,'Adamson2016',gset,own,'library_size');cat('Adamson own exclusions',length(own),'elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'seconds\n');flush.console()
# Frangieh: CSC gene columns; restrict to the frozen gene cohort and context controls.
x<-fread(file.path(root,'candidate_validation/frangieh_cell_scores_harmonized.tsv.gz'));gset<-fread(file.path(root,'candidate_validation/harmonized_frangieh_candidates.tsv'))$gene;needed<-intersect(gset,upr)
f<-H5File$new('/Volumes/wang/paper13/raw/perturbation/SCP1064/FrangiehIzar2021_RNA.h5ad','r');gen<-as.character(annot(f,'var','gene_symbol'));ids<-as.character(annot(f,'obs','cell_name'));ri<-match(x$cell_id,ids);stopifnot(!anyNA(ri));ptr<-f[['X']][['indptr']][];own<-list()
for(g in needed){j<-match(g,gen);stopifnot(!is.na(j));raw<-numeric(length(ids));lo<-ptr[j]+1;hi<-ptr[j+1L];if(hi>=lo)raw[f[['X']][['indices']][lo:hi]+1L]<-f[['X']][['data']][lo:hi];own[[g]]<-log1p(1e4*raw[ri]/x$library_size_common);ii<-which(x$gene==g);stopifnot(max(abs(own[[g]][ii]-x$target_logcp10k[ii]))<1e-9)}
f$close_all();process(x,'Frangieh2021',gset,own,'library_size_common');cat('Frangieh own exclusions',length(own),'elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'seconds\n')
r<-rbindlist(records);stopifnot(all(is.finite(r$effect_z)),max(abs(r[target_in_upr==FALSE]$effect_z-r[target_in_upr==FALSE]$original_effect_z))<1e-12)
fwrite(r,file.path(out,'UPR_own_target_excluded_effects.tsv'),sep='\t');fwrite(rbindlist(gemrecords),file.path(out,'Replogle_UPR_own_target_excluded_GEM_effects.tsv'),sep='\t')
writeLines(c('The frozen91-member UPR mean is replaced, only when the perturbed target is a member, by(91*original_score-own_logCP10k)/90 for BOTH treated cells and every matched control cell. Control mean AND SD are target-specific after exclusion. All nonmember target effects are unchanged. The7,848-gene library denominator is unchanged; raw duplicate feature counts are aggregated before log. Cohorts are unchanged and no ranking is updated. RPE1 expression values are never accessed by this script.',paste('Elapsed seconds',as.numeric(difftime(Sys.time(),started,units='secs')))),file.path(out,'UPR_own_target_exclusion_provenance.txt'))
