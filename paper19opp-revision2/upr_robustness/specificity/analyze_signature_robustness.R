#!/usr/bin/env Rscript
# No gene set or target population is chosen from response outcomes.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);project<-normalizePath(if(length(a))a[1]else'paper19opp-revision2');stage<-if(length(a)>1)a[2]else'primary'
out<-file.path(project,'upr_robustness/specificity');mod<-file.path(project,'upr_robustness/modules');source(file.path(out,'statistics_core.R'))
stopifnot(digest(file.path(out,'inference_protocol.txt'),file=TRUE,algo='sha256')=='c093d504c9bba1300130fe20a06a5cd643edce717c6330253662071728d9b23a')
set.seed(2026110611L);B<-4999L;P<-19999L
defs<-fread(file.path(mod,'locked_axis_definitions.tsv'));stopifnot(nrow(defs)==16L)
defs[,subfamily_id:=fifelse(parent_module%in%c('Reactome_UPR','GO_ER_UPR'),'independent_UPR_definition',fifelse(variant=='full','functional_full',fifelse(variant=='Hallmark_disjoint','functional_Hallmark_disjoint','functional_unique')))]
panel<-fread(file.path(out,'fixed_target_panels.tsv'));comps<-c('RPE1','Adamson','Frangieh_Control','Frangieh_IFNg','Frangieh_Coculture')
stopifnot(nrow(panel)==523L,setequal(panel$comparison,comps))
write_session<-function(){writeLines(capture.output(sessionInfo()),file.path(out,paste0(stage,'_session_info.txt')))}

if(stage=='primary'){
 files<-file.path(mod,'scored',paste0(c('K562','RPE1','Adamson','Frangieh'),'_gene_effects.tsv.gz'));stopifnot(all(file.exists(files)))
 # Record all exact source paths/checksums before loading these new response values.
 hh<-data.table(file=sub(paste0('^',project,'/'),'',files),sha256=vapply(files,function(f)digest(f,file=TRUE,algo='sha256'),character(1)));fwrite(hh,file.path(out,'response_input_hashes.tsv'),sep='\t')
 dd<-lapply(seq_along(files),function(i){z<-fread(files[i]);z[,study:=c('K562','RPE1','Adamson','Frangieh')[i]];z});d<-rbindlist(dd,use.names=TRUE,fill=TRUE)
 stopifnot(all(c('gene','axis','mode','delta_z','n_target','control_sd','control_expression','estimable')%in%names(d)))
 d<-d[mode=='exact_Lminus_target'&axis%in%c(defs$axis,'Hallmark_UPR90','core_apoptosis','core_proliferation','core_stress','E2F')]
 d[,comparison:=fifelse(study=='RPE1','RPE1',fifelse(study=='Adamson','Adamson',fifelse(study=='Frangieh'&context=='Control','Frangieh_Control',fifelse(study=='Frangieh'&context%in%c('IFNγ','IFNg','IFN-gamma'),'Frangieh_IFNg',fifelse(study=='Frangieh'&context%in%c('Co-culture','Coculture','co-culture'),'Frangieh_Coculture',NA_character_)))))]
 if(any(d$study=='Frangieh'&is.na(d$comparison)))stop('Unknown Frangieh context metadata')
 source_d<-d[study=='K562'];validation_d<-d[study!='K562'];stopifnot(!anyDuplicated(source_d[,.(gene,axis)]),!anyDuplicated(validation_d[,.(gene,axis,comparison)]))
 cols<-intersect(c('gene','axis','delta_z','n_target','n_control','control_sd','control_expression','control_detection','measured_members','remaining_members','estimable'),names(d))
 sd<-source_d[,..cols];setnames(sd,setdiff(cols,c('gene','axis')),paste0('source_',setdiff(cols,c('gene','axis'))));setnames(sd,'source_delta_z','source_deltaZ')
 vd<-validation_d[,c('comparison',cols),with=FALSE];setnames(vd,setdiff(cols,c('gene','axis')),paste0('validation_',setdiff(cols,c('gene','axis'))));setnames(vd,'validation_delta_z','validation_deltaZ')
 pairs<-panel[,.(axis=defs$axis),by=.(comparison,validation_study,validation_context,gene)];pairs<-merge(pairs,sd,by=c('gene','axis'),all.x=TRUE,sort=FALSE);pairs<-merge(pairs,vd,by=c('comparison','gene','axis'),all.x=TRUE,sort=FALSE)
 pairs<-merge(pairs,defs[,.(axis,parent_module,variant,subfamily_id,full_canonical_members,measured_members)],by='axis',all.x=TRUE,sort=FALSE)
 for(st in c('source','validation')){
  base<-if(st=='source')source_d else validation_d;cc<-base[axis%in%c('core_apoptosis','core_proliferation','core_stress','E2F')]
  cc<-if(st=='source')dcast(cc,gene~axis,value.var='delta_z')else dcast(cc,gene+comparison~axis,value.var='delta_z')
  cc[,generic_magnitude:=sqrt(core_apoptosis^2+core_proliferation^2+core_stress^2)];cc<-cc[,c(if(st=='source')'gene'else c('gene','comparison'),'generic_magnitude','E2F'),with=FALSE];setnames(cc,c('generic_magnitude','E2F'),paste0(st,c('_generic_magnitude','_E2F')));pairs<-merge(pairs,cc,by=if(st=='source')'gene'else c('gene','comparison'),all.x=TRUE,sort=FALSE)
 }
 pairs[,estimable_pair:=is.finite(source_deltaZ)&is.finite(validation_deltaZ)&source_estimable&validation_estimable]
 pairs[,`:=`(source_log_count=log(source_n_target),validation_log_count=log(validation_n_target))];setorder(pairs,comparison,axis,gene);fwrite(pairs,file.path(out,'all_fixed_panel_signature_effects.tsv.gz'),sep='\t')
 rm(d,dd,sd,vd,source_d,validation_d);gc(FALSE)
 allres<-allboot<-allperm<-allloo<-list();idx<-0L;cache<-new.env(parent=emptyenv())
 for(cp in comps)for(ax in defs$axis){
  idx<-idx+1L;q<-pairs[comparison==cp&axis==ax&estimable_pair==TRUE];setorder(q,gene);key<-digest(q$gene)
  if(!exists(key,envir=cache,inherits=FALSE))assign(key,list(perm=sp_perm_index(nrow(q),P),boot=if(nrow(q))replicate(B,sample.int(nrow(q),nrow(q),replace=TRUE))else matrix(integer(),0,B)),envir=cache)
  ci<-get(key,envir=cache);fit<-sp_pair(q,B,P,ci$perm,ci$boot)
  rr<-fit$results;rr[,`:=`(comparison=cp,axis=ax,subfamily_id=defs[axis==ax,subfamily_id],fixed_target_genes=panel[comparison==cp,.N],measured_readout_members=defs[axis==ax,measured_members])];allres[[idx]]<-rr
  bb<-fit$bootstrap;bb[,`:=`(comparison=cp,axis=ax)];allboot[[idx]]<-bb
  pp<-fit$permutation;pp[,`:=`(comparison=cp,axis=ax)];allperm[[idx]]<-pp
  ll<-fit$leave_one_out;ll[,`:=`(comparison=cp,axis=ax)];allloo[[idx]]<-ll
  cat(format(Sys.time(),tz='UTC',usetz=TRUE),cp,ax,'n=',nrow(q),'rho=',rr[metric=='rho',estimate],'\n');flush.console()
 }
 res<-rbindlist(allres,fill=TRUE);res[metric=='rho',BH80:=sp_adjust(raw_permutation_p,80L)];res[metric=='rho',BH20_subfamily:=sp_adjust(raw_permutation_p,20L),by=subfamily_id];res[metric%in%c('direction','C'),BH160_secondary:=sp_adjust(raw_permutation_p,160L)]
 fwrite(res,file.path(out,'all80_association_direction_consistency.tsv'),sep='\t');fwrite(res[metric=='rho'],file.path(out,'primary80_rank_associations.tsv'),sep='\t');fwrite(res[metric%in%c('direction','C')],file.path(out,'secondary160_direction_consistency.tsv'),sep='\t');fwrite(res[metric%in%c('kappa','MCC')],file.path(out,'chance_adjusted_direction_agreement.tsv'),sep='\t')
 fwrite(rbindlist(allboot),file.path(out,'gene_bootstrap_draws.tsv.gz'),sep='\t');fwrite(rbindlist(allperm),file.path(out,'gene_permutation_draws.tsv.gz'),sep='\t');fwrite(rbindlist(allloo,fill=TRUE),file.path(out,'leave_one_target_out.tsv.gz'),sep='\t')
 meta<-rbindlist(lapply(defs$axis,function(ax){m<-sp_meta(pairs[axis==ax&estimable_pair==TRUE],B);m$summary[,axis:=ax];m$summary}),fill=TRUE);fwrite(meta,file.path(out,'study_balanced_signed_effect_summary.tsv'),sep='\t')
 saveRDS(list(pairs=pairs,definitions=defs,protocol_hash=digest(file.path(out,'inference_protocol.txt'),file=TRUE,algo='sha256'),B=B,P=P),file.path(out,'primary_signature_state.rds'),compress='xz')
 print(res[metric=='rho',.(comparison,axis,genes,estimate,gene_bootstrap_low,gene_bootstrap_high,raw_permutation_p,BH80)]);write_session()
}else if(stage=='controls'){
 z<-readRDS(file.path(out,'primary_signature_state.rds'));pairs<-z$pairs;models<-list(cell_support=c('source_log_count','validation_log_count'),control_expression=c('source_control_expression','validation_control_expression'),generic_magnitude=c('source_generic_magnitude','validation_generic_magnitude'),E2F=c('source_E2F','validation_E2F'));ans<-list();ii<-0L
 for(cp in comps)for(ax in defs$axis)for(mm in names(models)){
  q<-pairs[comparison==cp&axis==ax&estimable_pair==TRUE];vars<-models[[mm]];covs<-as.matrix(q[,..vars]);fit<-sp_partial(q$source_deltaZ,q$validation_deltaZ,covs,B,P);rr<-fit$summary;rr[,`:=`(comparison=cp,axis=ax,control_model=mm,fixed_target_genes=panel[comparison==cp,.N])];ii<-ii+1L;ans[[ii]]<-rr
  cat(cp,ax,mm,'n=',rr$genes,'partial=',rr$estimate,'\n');flush.console()
 }
 rr<-rbindlist(ans,fill=TRUE);rr[,BH320:=sp_adjust(raw_permutation_p,320L)];fwrite(rr,file.path(out,'all320_minimal_covariate_sensitivities.tsv'),sep='\t');write_session()
}else stop('Unknown stage; use primary or controls')
