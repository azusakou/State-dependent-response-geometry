#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
a<-commandArgs(TRUE);project<-normalizePath(if(length(a))a[1]else'paper19opp-revision2');out<-file.path(project,'upr_robustness/specificity');mod<-file.path(project,'upr_robustness/modules');source(file.path(out,'statistics_core.R'));set.seed(2026110630L);B<-4999L;P<-19999L
studies<-c('K562','RPE1','Adamson','Frangieh');panels<-fread(file.path(out,'fixed_target_panels.tsv'));defs<-fread(file.path(mod,'locked_axis_definitions.tsv'));axes<-c(defs$axis,'Hallmark_UPR90');stopifnot(length(axes)==17L)
available<-fread(file.path(out,'sequencing_depth_schema_availability.tsv'));stopifnot(all(available$recorded_UMI_available),all(available$common_library_available))
files<-c(file.path(mod,'scored',paste0(studies,'_cell_cache.rds')),file.path(mod,'scored',paste0(studies,'_gene_effects.tsv.gz')),file.path(out,c('sequencing_depth_addendum.txt','analyze_sequencing_depth.R','statistics_core.R','fixed_target_panels.tsv','sequencing_depth_schema_availability.tsv')))
# No depth mean has been aggregated when this manifest is written.
hashes<-data.table(file=sub(paste0('^',project,'/'),'',files),sha256=vapply(files,function(f)digest(f,file=TRUE,algo='sha256'),character(1)));fwrite(hashes,file.path(out,'sequencing_depth_preaggregation_hashes.tsv'),sep='\t');writeLines(as.character(Sys.time()),file.path(out,'sequencing_depth_preaggregation_time.txt'))
summarize<-function(v){ok<-is.finite(v)&v>=0;c(mean=if(any(ok))mean(log1p(v[ok]))else NA_real_,valid=sum(ok),missing=sum(!is.finite(v)),invalid=sum(is.finite(v)&v<0))}
allcov<-list()
for(st in studies){
 ca<-readRDS(file.path(mod,'scored',paste0(st,'_cell_cache.rds')));md<-as.data.table(ca$metadata);rm(ca);gc(FALSE);stopifnot(all(c('UMI_count','common_library','gene','context','is_control')%in%names(md)))
 allowed<-if(st=='K562')md$gene%in%unique(panels$gene)else paste(md$gene,md$context,sep='|')%in%panels[validation_study==st,paste(gene,validation_context,sep='|')]
 md<-md[is_control==FALSE&allowed];d<-md[,{
  u<-summarize(UMI_count);l<-summarize(common_library);.(n_target_cells=.N,mean_log1p_recorded_UMI=u['mean'],n_valid_UMI=u['valid'],n_missing_UMI=u['missing'],n_invalid_UMI=u['invalid'],mean_log1p_common_library=l['mean'],n_valid_common_library=l['valid'],n_missing_common_library=l['missing'],n_invalid_common_library=l['invalid'])
 },by=.(gene,context)];d[,study:=st];allcov[[st]]<-d;cat('Depth covariates aggregated:',st,nrow(d),'gene-context rows\n');flush.console()
}
cov<-rbindlist(allcov,use.names=TRUE);setorder(cov,study,context,gene);fwrite(cov,file.path(out,'sequencing_depth_covariates.tsv'),sep='\t')
dd<-lapply(studies,function(st){d<-fread(file.path(mod,'scored',paste0(st,'_gene_effects.tsv.gz')),select=c('context','gene','axis','mode','delta_z','n_target'));d[,study:=st];d[(axis%in%defs$axis&mode=='exact_Lminus_target')|(axis=='Hallmark_UPR90'&mode=='fixed_L6780')]});effects<-rbindlist(dd)
effects[,comparison:=fifelse(study=='RPE1','RPE1',fifelse(study=='Adamson','Adamson',fifelse(study=='Frangieh'&context=='Control','Frangieh_Control',fifelse(study=='Frangieh'&context=='IFNγ','Frangieh_IFNg',fifelse(study=='Frangieh'&context=='Co-culture','Frangieh_Coculture',NA_character_)))))]
ss<-effects[study=='K562',.(gene,axis,source_deltaZ=delta_z,source_n_target=n_target,response_mode=mode)];vv<-effects[study!='K562',.(gene,comparison,axis,validation_deltaZ=delta_z,validation_n_target=n_target)];stopifnot(!anyDuplicated(ss[,.(gene,axis)]),!anyDuplicated(vv[,.(gene,comparison,axis)]))
pairs<-panels[,.(axis=axes),by=.(comparison,validation_study,validation_context,gene)];pairs<-merge(pairs,ss,by=c('gene','axis'),all.x=TRUE);pairs<-merge(pairs,vv,by=c('gene','axis','comparison'),all.x=TRUE)
sc<-cov[study=='K562',.(gene,source_mean_log1p_recorded_UMI=mean_log1p_recorded_UMI,source_mean_log1p_common_library=mean_log1p_common_library,source_depth_cells=n_target_cells)];vc<-merge(panels,cov,by.x=c('gene','validation_study','validation_context'),by.y=c('gene','study','context'),all.x=TRUE)
vc<-vc[,.(gene,comparison,validation_mean_log1p_recorded_UMI=mean_log1p_recorded_UMI,validation_mean_log1p_common_library=mean_log1p_common_library,validation_depth_cells=n_target_cells)];pairs<-merge(pairs,sc,by='gene',all.x=TRUE);pairs<-merge(pairs,vc,by=c('gene','comparison'),all.x=TRUE)
stopifnot(nrow(pairs)==17L*523L,all(pairs$source_n_target==pairs$source_depth_cells),all(pairs$validation_n_target==pairs$validation_depth_cells));fwrite(pairs,file.path(out,'sequencing_depth_model_inputs.tsv.gz'),sep='\t')
models<-list(recorded_UMI=c('source_mean_log1p_recorded_UMI','validation_mean_log1p_recorded_UMI'),common_library=c('source_mean_log1p_common_library','validation_mean_log1p_common_library'));rows<-list();i<-0L
for(cp in unique(panels$comparison))for(ax in axes)for(model in names(models)){
 d<-pairs[comparison==cp&axis==ax];fields<-models[[model]];fit<-sp_partial(d$source_deltaZ,d$validation_deltaZ,as.matrix(d[,..fields]),B,P);r<-fit$summary;r[,`:=`(comparison=cp,axis=ax,response_mode=unique(d$response_mode),depth_model=model,fixed_target_genes=nrow(d),family='captured_RNA_depth170')];i<-i+1L;rows[[i]]<-r;cat(cp,ax,model,'n=',r$genes,'partial=',r$estimate,'\n');flush.console()
}
res<-rbindlist(rows,fill=TRUE);stopifnot(nrow(res)==170L);res[,BH170:=sp_adjust(raw_permutation_p,170L)];fwrite(res,file.path(out,'sequencing_depth170_conditional_associations.tsv'),sep='\t')
schema<-data.table(column=names(cov),meaning=c('Assigned perturbation target','Existing context label','Number of treated cells in unchanged panel','Mean log1p recorded UMI count over finite nonnegative values','Number contributing to recorded-UMI mean','Missing/nonfinite recorded-UMI values','Finite negative recorded-UMI values','Mean log1p exact frozen6780-gene library count','Number contributing to common-library mean','Missing/nonfinite common-library values','Finite negative common-library values','Existing source study'))
fwrite(schema,file.path(out,'sequencing_depth_covariate_schema.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(out,'sequencing_depth_session_info.txt'));cat('Completed all170 conditional depth diagnostics; primary families unchanged.\n')
