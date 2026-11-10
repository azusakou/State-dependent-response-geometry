#!/usr/bin/env Rscript
# Paired specification diagnostics and normalization-matched negative controls.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);project<-normalizePath(if(length(a))a[1]else'paper19opp-revision2');stage<-if(length(a)>1)a[2]else'paired'
out<-file.path(project,'upr_robustness/specificity');mod<-file.path(project,'upr_robustness/modules');source(file.path(out,'statistics_core.R'));set.seed(2026110612L);B<-4999L
z<-readRDS(file.path(out,'primary_signature_state.rds'));defs<-z$definitions;d<-z$pairs;panels<-fread(file.path(out,'fixed_target_panels.tsv'));comps<-unique(panels$comparison)
if(stage=='paired'){
 paired<-function(a,b){
  namesa<-names(a);a<-a[,.(gene,x=source_deltaZ,y=validation_deltaZ)];b<-b[,.(gene,x2=source_deltaZ,y2=validation_deltaZ)];q<-merge(a,b,by='gene');q<-q[is.finite(x)&is.finite(y)&is.finite(x2)&is.finite(y2)];n<-nrow(q)
  if(n<3L)return(data.table(genes=n,rho1=NA_real_,rho2=NA_real_,rho_difference=NA_real_,rho_difference_low=NA_real_,rho_difference_high=NA_real_,C1=NA_real_,C2=NA_real_,C_difference=NA_real_,C_difference_low=NA_real_,C_difference_high=NA_real_))
  r1<-sp_rho(q$x,q$y);r2<-sp_rho(q$x2,q$y2);c1<-mean(sign(q$x)*q$y);c2<-mean(sign(q$x2)*q$y2)
  boot<-vapply(seq_len(B),function(j){ii<-sample.int(n,n,replace=TRUE);c(rho=sp_rho(q$x[ii],q$y[ii])-sp_rho(q$x2[ii],q$y2[ii]),C=mean(sign(q$x[ii])*q$y[ii])-mean(sign(q$x2[ii])*q$y2[ii]))},numeric(2));cr<-sp_ci(boot['rho',]);cc<-sp_ci(boot['C',])
  data.table(genes=n,rho1=r1,rho2=r2,rho_difference=r1-r2,rho_difference_low=cr[1],rho_difference_high=cr[2],C1=c1,C2=c2,C_difference=c1-c2,C_difference_low=cc[1],C_difference_high=cc[2])
 }
 ans<-list();i<-0L
 for(cp in comps)for(pa in unique(defs$parent_module)){
  full<-defs[parent_module==pa&variant=='full',axis]
  for(alt in defs[parent_module==pa&variant!='full',axis]){i<-i+1L;rr<-paired(d[comparison==cp&axis==full&estimable_pair==TRUE],d[comparison==cp&axis==alt&estimable_pair==TRUE]);rr[,`:=`(comparison=cp,axis1=full,axis2=alt,interpretation='paired specification diagnostic; no additional significance test')];ans[[i]]<-rr}
 }
 fwrite(rbindlist(ans,fill=TRUE),file.path(out,'paired_full_disjoint_unique_specification_differences.tsv'),sep='\t')
 ans<-list();i<-0L
 for(ax in defs$axis)for(cp in c('Frangieh_IFNg','Frangieh_Coculture')){
  i<-i+1L;rr<-paired(d[comparison=='Frangieh_Control'&axis==ax&estimable_pair==TRUE],d[comparison==cp&axis==ax&estimable_pair==TRUE]);rr[,`:=`(axis=ax,comparison1='Frangieh_Control',comparison2=cp,interpretation='direct paired context diagnostic; no additional significance test')];ans[[i]]<-rr
 }
 fwrite(rbindlist(ans,fill=TRUE),file.path(out,'paired_Frangieh_context_differences.tsv'),sep='\t')
}else if(stage=='nulls'){
 files<-file.path(mod,'scored',paste0(c('K562','RPE1','Adamson','Frangieh'),'_gene_effects.tsv.gz'));dd<-lapply(seq_along(files),function(i){x<-fread(files[i],select=c('study','context','gene','axis','mode','delta_z','estimable'));x[,study:=c('K562','RPE1','Adamson','Frangieh')[i]];x});x<-rbindlist(dd,use.names=TRUE,fill=TRUE)[mode=='fixed_L6780']
 x[,comparison:=fifelse(study=='RPE1','RPE1',fifelse(study=='Adamson','Adamson',fifelse(study=='Frangieh'&context=='Control','Frangieh_Control',fifelse(study=='Frangieh'&context%in%c('IFNγ','IFNg','IFN-gamma'),'Frangieh_IFNg',fifelse(study=='Frangieh'&context%in%c('Co-culture','Coculture','co-culture'),'Frangieh_Coculture',NA_character_)))))]
 nm<-unique(fread(file.path(mod,'null100_per_readout_memberships.tsv.gz'))[,.(family,module)]);stopifnot(nrow(nm)==800L)
 availability<-copy(defs);availability[,size_expression_detection_matched_nulls:=vapply(axis,function(ax)nm[family==ax,.N],integer(1))];availability[,interpretation:=fifelse(size_expression_detection_matched_nulls>0,'fixedL real/null matched assay available','no size-matched null family; specificity unresolved')];fwrite(availability,file.path(out,'readout_specificity_availability.tsv'),sep='\t')
 scores<-list();i<-0L
 for(cp in comps){
  gg<-panels[comparison==cp,gene];s<-x[study=='K562'&gene%in%gg,.(gene,axis,s=delta_z)];v<-x[comparison==cp&study!='K562'&gene%in%gg,.(gene,axis,v=delta_z)];q<-merge(s,v,by=c('gene','axis'),all=TRUE)
  for(ax in c(unique(nm$module),unique(nm$family),'Hallmark_UPR90')){
   qq<-q[axis==ax&is.finite(s)&is.finite(v)];i<-i+1L;ag<-sp_agreement(qq$s,qq$v);scores[[i]]<-data.table(comparison=cp,axis=ax,genes=nrow(qq),rho=sp_rho(qq$s,qq$v),C=if(nrow(qq))mean(sign(qq$s)*qq$v)else NA_real_,direction=ag['concordance'],kappa=ag['kappa'],MCC=ag['MCC'])
  }
 }
 scores<-rbindlist(scores);fwrite(scores,file.path(out,'all_fixed_library_real_null_readout_metrics.tsv.gz'),sep='\t');cal<-list();i<-0L
 for(cp in comps)for(fa in unique(nm$family))for(metric in c('rho','C')){
  nn<-scores[comparison==cp&axis%in%nm[family==fa,module],get(metric)];nn<-nn[is.finite(nn)];obs<-scores[comparison==cp&axis==fa,get(metric)];i<-i+1L;cal[[i]]<-data.table(comparison=cp,axis=fa,metric=metric,observed_fixed_library=obs,valid_null_modules=length(nn),upper_tail_fraction=if(length(nn)&&is.finite(obs))(1+sum(nn>=obs-1e-12))/(1+length(nn))else NA_real_,null_median=if(length(nn))median(nn)else NA_real_,null_low=if(length(nn))quantile(nn,.025)else NA_real_,null_high=if(length(nn))quantile(nn,.975)else NA_real_,interpretation='matched fixedL real/null assay; descriptive finite-control calibration')
 }
 fwrite(rbindlist(cal),file.path(out,'fixed_library_readout_specificity_calibration.tsv'),sep='\t')
}else stop('Use paired or nulls stage')
writeLines(capture.output(sessionInfo()),file.path(out,paste0(stage,'_session_info.txt')))
