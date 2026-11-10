#!/usr/bin/env Rscript
# Run only after the source manifest is frozen. Validation is loaded exactly once.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);if(length(a)<4L)stop('Usage: Rscript run_phase2_statistics.R source_directory RPE1_cache.rds output_directory pre_RPE1_manifest.tsv')
srcdir<-normalizePath(a[1]);rpefile<-a[2];out<-a[3];manifest<-normalizePath(a[4]);dir.create(out,recursive=TRUE,showWarnings=FALSE)
ff<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',ff)))
source(file.path(here,'phase2_statistics.R'));source(file.path(srcdir,'phase2_helpers.R'))
verify_manifest<-function(f){
 h<-fread(f);pc<-intersect(c('file','path','relative_path'),names(h));hc<-names(h)[tolower(names(h))%in%c('sha256','sha_256')];stopifnot(length(pc)>=1L,length(hc)==1L)
 studyroot<-dirname(dirname(dirname(dirname(srcdir))))
 for(i in seq_len(nrow(h))){p<-as.character(h[[pc[1]]][i]);if(!grepl('^[A-Za-z0-9]{64}$',h[[hc]][i]))next
  candidates<-if(startsWith(p,'/'))p else c(file.path(dirname(f),p),file.path(srcdir,p),file.path(studyroot,p),p);existing<-candidates[file.exists(candidates)];if(!length(existing))stop(paste('Missing locked file',p));got<-digest(existing[1],file=TRUE,algo='sha256');if(!identical(got,as.character(h[[hc]][i])))stop(paste('Hash mismatch',p))
 }
}
verify_manifest(manifest)
# Reuse the declared scoring helper but materialize GEM effects only for real UPR.
fast_effects<-function(cache,spec){
 base<-control_base_moments(cache$control_positive-cache$control_negative,cache$control_cells$gem_group)
 units<-as.data.table(cache$unit_metadata);counts<-units$n_cells;ans<-gans<-vector('list',length(cache$analysis_genes));iiout<-0L
 for(g in cache$analysis_genes){
  iiout<-iiout+1L;ii<-which(units$gene==g);cm<-target_control_moments(g,spec,cache$control_positive,cache$control_negative,cache$control_target_log,cache$control_cells$gem_group,base);gi<-match(units$gem_group[ii],base$gems);stopifnot(!anyNA(gi));diff<-cache$unit_score_sums[ii,,drop=FALSE]/counts[ii]-cm$mean[gi,,drop=FALSE];raw<-colSums(diff*counts[ii])/sum(counts[ii]);f<-colMeans(diff[counts[ii]>=5,,drop=FALSE]>0)
  ans[[iiout]]<-data.table(gene=g,axis=names(spec),delta_raw=raw,control_sd=cm$sd,deltaZ=raw/cm$sd,f_positive_GEM=f,source_cells=sum(counts[ii]),source_GEM_labels=length(ii),source_GEM_ge5=sum(counts[ii]>=5))
  k<-match('UPR',names(spec));gans[[iiout]]<-data.table(gene=g,gem=as.character(units$gem_group[ii]),target_cells=counts[ii],contrast_z=diff[,k]/cm$sd[k])
 }
 list(gene=rbindlist(ans),gem=rbindlist(gans))
}
si<-fread(file.path(srcdir,'source_rank_inputs.tsv'))
s<-si[,.(gene,source_delta_z=source_UPR_deltaZ,priority_score=UPR_priority,stratum=permutation_stratum)]
allrank<-fread(file.path(srcdir,'all_source_comparator_rankings.tsv'));ranking<-allrank[method=='UPR_priority',.(gene,rank)]
matchall<-fread(file.path(srcdir,'frozen_background_matches.tsv'));matches<-matchall[list=='UPR',.(selected_gene,background_gene)]
sg<-fread(file.path(srcdir,'source_phase2_GEM_effects.tsv.gz'),select=c('gene','gem_group','axis','n_target','deltaZ'));sg<-sg[axis=='UPR',.(gene,gem=as.character(gem_group),target_cells=n_target,contrast_z=deltaZ)]
spec<-readRDS(file.path(srcdir,'phase2_axis_spec.rds'))
stopifnot(file.exists(rpefile));rpe<-readRDS(rpefile);ve<-fast_effects(rpe,spec);rm(rpe);gc(FALSE)
fwrite(ve$gene,file.path(out,'RPE1_gene_effects_all105.tsv.gz'),sep='\t');fwrite(ve$gem,file.path(out,'RPE1_real_UPR_GEM_effects.tsv.gz'),sep='\t')
v<-ve$gene[axis=='UPR',.(gene,validation_delta_z=deltaZ)]
stopifnot(nrow(s)==438L,nrow(v)==438L,setequal(s$gene,v$gene))
primary<-phase2_primary(s,v,ranking,matches,sg,ve$gem)
fwrite(primary$endpoints,file.path(out,'primary_Holm3_endpoints.tsv'),sep='\t');fwrite(primary$diagnostics,file.path(out,'primary_diagnostics.tsv'),sep='\t');fwrite(primary$bootstrap,file.path(out,'primary_conditional_bootstrap_draws.tsv.gz'),sep='\t');fwrite(primary$permutation_null,file.path(out,'primary_permutation_nulls.tsv.gz'),sep='\t')
ccr<-fread(file.path(srcdir,'complete_covariate_source_rankings.tsv'));stopifnot('rank_cc'%in%names(ccr));ccgenes<-sort(unique(ccr$gene));cmp<-ccr[,.(method,gene,rank=rank_cc)]
baselines<-phase2_baselines(primary,s,cmp,ccgenes,'raw_UPR');fwrite(baselines,file.path(out,'complete_case_baseline_utility.tsv'),sep='\t')
nr<-fread(file.path(srcdir,'null100_source_rankings.tsv.gz'));nsource<-nr[,.(module_id=module,gene,source_delta_z=source_deltaZ)];nranking<-nr[,.(module_id=module,gene,rank,priority_score=score)]
nvalidation<-ve$gene[grepl('^null_',axis),.(module_id=axis,gene,validation_delta_z=deltaZ)]
nmatches<-matchall[grepl('^null_',list),.(module_id=list,selected_gene,background_gene)]
ns<-phase2_module_specificity(primary,s,v,nsource,nvalidation,nmatches,nranking);fwrite(ns$all_modules,file.path(out,'primary_all100_module_controls.tsv'),sep='\t');fwrite(ns$calibration,file.path(out,'primary_module_specificity_calibration.tsv'),sep='\t')
# Uniform source-only complete-covariate sensitivity was separately frozen before access.
ccmatches<-fread(file.path(srcdir,'complete_covariate_background_matches.tsv'));ccranking<-ccr[method=='UPR_priority',.(gene,rank=rank_cc)];ccs<-s[gene%in%ccgenes];ccv<-v[gene%in%ccgenes]
ccprimary<-phase2_primary(ccs,ccv,ccranking,ccmatches[list=='CC_UPR',.(selected_gene,background_gene)],sg[gene%in%ccgenes],ve$gem[gene%in%ccgenes],seed=20261108L)
fwrite(ccprimary$endpoints,file.path(out,'complete_covariate_sensitivity_endpoints.tsv'),sep='\t')
ccnr<-fread(file.path(srcdir,'complete_covariate_null100_rankings.tsv.gz'));rkfield<-if('rank_cc'%in%names(ccnr))'rank_cc'else'rank';ccnsource<-ccnr[,.(module_id=module,gene,source_delta_z=source_deltaZ)];ccnranking<-ccnr[,.(module_id=module,gene,rank=get(rkfield),priority_score=score)]
ccnm<-ccmatches[grepl('^CC_null_',list),.(module_id=sub('^CC_','',list),selected_gene,background_gene)]
ccns<-phase2_module_specificity(ccprimary,ccs,ccv,ccnsource,nvalidation[gene%in%ccgenes],ccnm,ccnranking)
fwrite(ccns$all_modules,file.path(out,'complete_covariate_all100_module_controls.tsv'),sep='\t');fwrite(ccns$calibration,file.path(out,'complete_covariate_module_specificity_calibration.tsv'),sep='\t')
# Conditional original-priority estimates on the stricter metadata-supported subgroup.
strictgenes<-si[strict_rpe_support==TRUE,gene];sx<-primary$point_source[strictgenes];sy<-primary$point_validation[strictgenes];set.seed(20261109L);b<-vapply(seq_len(2000L),function(i){jj<-sample.int(length(strictgenes),length(strictgenes),replace=TRUE);phase2_spearman(sx[jj],sy[jj])},numeric(1));qb<-phase2_quantile(b)
fwrite(data.table(genes=length(strictgenes),signed_rho=phase2_spearman(sx,sy),gene_bootstrap_low=qb[1],gene_bootstrap_high=qb[2]),file.path(out,'strict_RPE_GEM_support_sensitivity.tsv'),sep='\t')
# Small result cache shares GEM weights and fixed lists with the exact-depth/composition checks.
saveRDS(list(primary=primary,source=s,validation=v,ranking=ranking,matches=matches,source_gem=sg,validation_gem=ve$gem,source_metadata=si,validation_gene_effects=ve$gene,source_directory=srcdir,validation_cache=rpefile),file.path(out,'phase2_statistical_state.rds'),compress='xz')
writeLines(capture.output(sessionInfo()),file.path(out,'statistics_sessionInfo.txt'));print(primary$endpoints);print(baselines);print(ccns$calibration)
