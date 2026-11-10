#!/usr/bin/env Rscript
project<-normalizePath(commandArgs(TRUE)[1]);out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');source(file.path(out,'scoring_helpers.R'));verify_locks(out);dest<-file.path(out,'canonical_Hallmark94')
locks<-fread(file.path(dest,'technical_addendum_lock.tsv'));stopifnot(all(vapply(file.path(dest,locks$file),digest,character(1),file=TRUE,algo='sha256')==locks$sha256))
execution<-fread(file.path(dest,'execution_lock_manifest.tsv'));stopifnot(all(vapply(file.path(dest,execution$file),digest,character(1),file=TRUE,algo='sha256')==execution$sha256))
members<-fread(file.path(dest,'fixed_canonical_members.tsv'))$assay_gene;stopifnot(length(members)==94L)
panels<-fread(file.path(project,'paper19opp-revision2/upr_robustness/specificity/fixed_target_panels.tsv'));QA<-list()
sd_columns<-function(x){n<-nrow(x);sqrt(pmax((colSums(x*x)-colSums(x)^2/n)/(n-1),0))}
for(study_name in c('K562','RPE1','Adamson','Frangieh')){
 started<-Sys.time();cat('Canonical94 scoring',study_name,format(started,tz='UTC',usetz=TRUE),'\n');flush.console()
 c<-readRDS(file.path(out,'scored',paste0(study_name,'_cell_cache.rds')));additional<-readRDS(file.path(dest,paste0(study_name,'_additional_raw_counts.rds')));md<-c$metadata;ctl<-c$control_cells;ci<-c$control_rows;units<-c$unit_metadata
 stopifnot(identical(additional$cell_id,md$cell_id),identical(additional$source_row,md$source_row))
 old90<-c$real_spec$Hallmark_UPR90$positive;stopifnot(length(old90)==90L,setequal(c(old90,colnames(additional$raw_counts)),members))
 raw<-cbind(as.matrix(c$raw_real_union[,old90,drop=FALSE]),additional$raw_counts);raw<-raw[,members,drop=FALSE]
 spec<-list(Hallmark_UPR94_canonical=list(positive=members,negative=character()),Hallmark_UPR90_reference_check=list(positive=old90,negative=character()));coeff<-make_coeff(members,spec)
 L<-md$common_library;Lminus<-L;subtract<-md$gene%in%c$universe&!md$is_control;Lminus[subtract]<-Lminus[subtract]-md$assigned_target_count[subtract];stopifnot(all(Lminus>0))
 cell_fixed<-score_logs(t(log1p(raw/L*1e4)),members,spec,coeff,md$gene)
 cell_exact<-score_logs(t(log1p(raw/Lminus*1e4)),members,spec,coeff,md$gene)
 e_cell_fixed<-max(abs(cell_fixed[,2]-c$real_scores_fixed[,'Hallmark_UPR90']));e_cell_exact<-max(abs(cell_exact[,2]-c$real_scores_exact[,'Hallmark_UPR90']));stopifnot(e_cell_fixed<1e-10,e_cell_exact<1e-10)
 key<-paste(units$gene,units$context,units$gem,sep='|');unit_id<-match(paste(md$gene,md$context,md$gem,sep='|'),key);ti<-which(!md$is_control)
 sumf<-rowsum(cell_fixed[ti,,drop=FALSE],unit_id[ti],reorder=TRUE);sume<-rowsum(cell_exact[ti,,drop=FALSE],unit_id[ti],reorder=TRUE);stopifnot(identical(as.integer(rownames(sumf)),seq_len(nrow(units))))
 all_genes<-all_gems<-all_controls<-list()
 for(ctx in unique(md$context)){
  cc<-which(ctl$context==ctx);control_raw<-raw[ci[cc],,drop=FALSE];Lct<-ctl$common_library[cc];target_set<-if(study_name=='K562')c$targets else panels[validation_study==study_name&validation_context==ctx,gene]
  for(g in target_set){
   ui<-which(units$gene==g&units$context==ctx);stopifnot(length(ui)>0L);uu<-units[ui];w<-uu$n_target/sum(uu$n_target);target_raw<-c$control_target_raw[cc,g];Lctminus<-Lct;if(g%in%c$universe){stopifnot(all(is.finite(target_raw)));Lctminus<-Lctminus-target_raw};stopifnot(all(Lctminus>0))
   for(mode_name in c('fixed_L6780','exact_Lminus_target')){
    exact<-mode_name=='exact_Lminus_target';denom<-if(exact)Lctminus else Lct
    ctrl_score<-score_logs(t(log1p(control_raw/denom*1e4)),members,spec,coeff,rep(g,length(cc)))
    sigma<-sd_columns(ctrl_score);globalmean<-colMeans(ctrl_score)
    cm<-do.call(rbind,lapply(uu$gem,function(gm)if(c$anchor_type=='matched_GEM')colMeans(ctrl_score[ctl$gem[cc]==gm,,drop=FALSE])else globalmean))
    target_mean<-(if(exact)sume else sumf)[ui,,drop=FALSE]/uu$n_target;contrast<-sweep(target_mean-cm,2L,sigma,'/');delta<-colSums(contrast*w);supported<-uu$n_target>=5L;fp<-if(any(supported))colMeans(contrast[supported,,drop=FALSE]>0)else rep(NA_real_,2L)
    remain<-c(94L-as.integer(g%in%members),90L-as.integer(g%in%old90));valid<-is.finite(delta)&is.finite(sigma)&sigma>0&remain>0
    all_genes[[length(all_genes)+1L]]<-data.table(study=study_name,context=ctx,gene=g,axis=names(spec),mode=mode_name,delta_z=ifelse(valid,delta,NA_real_),mean_target=colSums(target_mean*w),mean_control=colSums(cm*w),global_control_mean=globalmean,control_sd=sigma,n_target=sum(uu$n_target),n_control=length(cc),n_gems=nrow(uu),supported_groups_ge5=sum(supported),f_positive=fp,measured_members=c(94L,90L),remaining_members=remain,target_member=c(g%in%members,g%in%old90),target_in_universe=g%in%c$universe,estimable=valid,anchor_type=c$anchor_type)
    for(k in seq_along(spec)){
     all_gems[[length(all_gems)+1L]]<-data.table(study=study_name,context=ctx,gene=g,axis=names(spec)[k],mode=mode_name,gem=uu$gem,n_target=uu$n_target,n_control=if(c$anchor_type=='matched_GEM')vapply(uu$gem,function(gm)sum(ctl$gem[cc]==gm),integer(1))else length(cc),contrast_z=contrast[,k],matched_control_mean=cm[,k],control_sd=sigma[k],genuine_GEM_metadata=study_name%in%c('K562','RPE1'),technical_group_type=if(study_name%in%c('K562','RPE1'))'GEM'else if(study_name=='Adamson')'capture_lane'else'pooled_context')
    }
    all_controls[[length(all_controls)+1L]]<-data.table(study=study_name,context=ctx,gene=g,axis=names(spec)[1],mode=mode_name,control_mean=globalmean[1],control_sd=sigma[1],n_control=length(cc))
   }
  }
 }
 gene_effects<-rbindlist(all_genes);gems<-rbindlist(all_gems);old_effects<-fread(file.path(out,'scored',paste0(study_name,'_gene_effects.tsv.gz')))[axis=='Hallmark_UPR90',.(study,context,gene,mode,old_delta_z=delta_z,old_n_target=n_target)]
 check<-merge(gene_effects[axis=='Hallmark_UPR90_reference_check'],old_effects,by=c('study','context','gene','mode'));stopifnot(nrow(check)==nrow(old_effects),all(check$n_target==check$old_n_target));e_gene<-max(abs(check$delta_z-check$old_delta_z),na.rm=TRUE);stopifnot(e_gene<1e-8)
 gene94<-gene_effects[axis=='Hallmark_UPR94_canonical'];gem94<-gems[axis=='Hallmark_UPR94_canonical']
 fwrite(gene94,file.path(dest,paste0(study_name,'_gene_effects.tsv.gz')),sep='\t');fwrite(gem94,file.path(dest,paste0(study_name,'_GEM_effects.tsv.gz')),sep='\t');fwrite(rbindlist(all_controls),file.path(dest,paste0(study_name,'_control_moments.tsv.gz')),sep='\t')
 saveRDS(list(study=study_name,axis='Hallmark_UPR94_canonical',members=members,metadata=md,real_scores_exact=matrix(cell_exact[,1],ncol=1,dimnames=list(NULL,'Hallmark_UPR94_canonical')),real_scores_fixed=matrix(cell_fixed[,1],ncol=1,dimnames=list(NULL,'Hallmark_UPR94_canonical')),control_rows=ci,unit_metadata=units,original_90_primary_unchanged=TRUE),file.path(dest,paste0(study_name,'_cell_scores.rds')),compress='gzip')
 QA[[study_name]]<-data.table(study=study_name,cells=nrow(md),controls=length(ci),canonical_members=length(members),additional_aliases=4L,old90_cell_error_fixed=e_cell_fixed,old90_cell_error_exact=e_cell_exact,old90_gene_effect_error=e_gene,same_target_support=all(check$n_target==check$old_n_target),passed=TRUE,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
 cat('Canonical94 complete',study_name,'reference90QAmax',e_gene,'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs')),2),'s\n');flush.console();rm(c,raw,cell_fixed,cell_exact,additional);gc(FALSE)
}
fwrite(rbindlist(QA),file.path(dest,'reference90_reconstruction_QA.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(dest,'scoring_sessionInfo.txt'))
