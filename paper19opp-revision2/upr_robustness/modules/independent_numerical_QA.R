#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(args[1]);out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');source(file.path(out,'scoring_helpers.R'));verify_locks(out)
results<-list()
for(study_name in c('K562','RPE1','Adamson','Frangieh')){
 cat('Independent raw-score checks',study_name,'\n');flush.console()
 dest<-file.path(out,'scored');c<-readRDS(file.path(dest,paste0(study_name,'_cell_cache.rds')));md<-c$metadata
 candidate<-which(!md$is_control&md$gene%in%unique(unlist(lapply(c$real_spec,function(s)c(s$positive,s$negative)))))
 # Deterministic checks include own-target member cells and controls, independent of their responses.
 take<-unique(c(head(candidate,12),head(which(!md$is_control),12),head(which(md$is_control),12)))
 counts<-as.matrix(c$raw_real_union[take,,drop=FALSE]);maxerr<-0
 for(mode in c('fixed_L6780','exact_Lminus_target')){
  direct<-matrix(NA_real_,length(take),length(c$real_spec))
  for(j in seq_along(take)){
   i<-take[j];L<-md$common_library[i];g<-md$gene[i]
   if(mode=='exact_Lminus_target'&&g%in%c$universe&&!md$is_control[i])L<-L-md$assigned_target_count[i]
   lg<-log1p(1e4*counts[j,]/L)
   for(k in seq_along(c$real_spec)){
    pos<-setdiff(c$real_spec[[k]]$positive,g);neg<-setdiff(c$real_spec[[k]]$negative,g)
    direct[j,k]<-if(length(pos))mean(lg[pos])else NA_real_;if(length(neg))direct[j,k]<-direct[j,k]-mean(lg[neg])
   }
  }
  saved<-if(mode=='fixed_L6780')c$real_scores_fixed[take,,drop=FALSE]else c$real_scores_exact[take,,drop=FALSE]
  error<-max(abs(direct-saved),na.rm=TRUE);stopifnot(error<1e-10);maxerr<-max(maxerr,error)
 }
 e<-fread(file.path(dest,paste0(study_name,'_gene_effects.tsv.gz')))[axis%in%names(c$real_spec)]
 g<-fread(file.path(dest,paste0(study_name,'_GEM_effects.tsv.gz')))
 rebuilt<-g[,.(reconstructed_delta_z=sum(contrast_z*n_target)/sum(n_target),reconstructed_n=sum(n_target)),by=.(study,context,gene,axis,mode)]
 z<-merge(e,rebuilt,by=c('study','context','gene','axis','mode'));stopifnot(nrow(z)==nrow(e),all(z$n_target==z$reconstructed_n));aggregate_error<-max(abs(z$delta_z-z$reconstructed_delta_z),na.rm=TRUE);stopifnot(aggregate_error<1e-10)
 results[[study_name]]<-data.table(study=study_name,scorable_cells=nrow(md),controls=sum(md$is_control),checked_cells=length(take),checked_real_axes=length(c$real_spec),maximum_direct_cellscore_error=maxerr,real_gene_mode_rows=nrow(e),maximum_GEM_reaggregation_error=aggregate_error,same_target_support=all(z$n_target==z$reconstructed_n),passed=TRUE)
 rm(c,e,g,z,rebuilt);gc(FALSE)
}
qa<-rbindlist(results);fwrite(qa,file.path(out,'independent_numerical_QA.tsv'),sep='\t');print(qa)
