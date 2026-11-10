#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(args[1]);study_name<-args[2]
out<-file.path(project,'paper19opp-revision2/upr_robustness/modules');source(file.path(out,'scoring_helpers.R'));verify_locks(out);dest<-file.path(out,'scored')
started<-Sys.time();cat('Summarizing',study_name,format(started,tz='UTC',usetz=TRUE),'\n');flush.console()
c<-readRDS(file.path(dest,paste0(study_name,'_cell_cache.rds')));md<-c$metadata;units<-c$unit_metadata;spec<-c$spec;real_spec<-c$real_spec;axes<-names(spec);real_axes<-names(real_spec);genes<-colnames(c$raw_real_union);universe<-c$universe;co<-make_coeff(genes,real_spec)
raw_ctl<-as.matrix(c$raw_real_union[c$control_rows,,drop=FALSE]);ctl<-c$control_cells;roster<-fread(file.path(out,'fixed_study_gene_rosters.tsv'))[study==study_name,gene]
panels<-fread(file.path(project,'paper19opp-revision2/upr_robustness/specificity/fixed_target_panels.tsv'))
all_effects<-all_gems<-all_members<-all_controls<-all_ctl_gems<-all_expr<-list();record<-0L
sd_columns<-function(x){n<-nrow(x);if(n<2)return(rep(NA_real_,ncol(x)));sqrt(pmax((colSums(x*x)-colSums(x)^2/n)/(n-1),0))}
unit_key<-paste(units$gene,units$context,units$gem,sep='|');cell_ui<-match(paste(md$gene,md$context,md$gem,sep='|'),unit_key)
tre<-which(!md$is_control)
realsums_fixed<-rowsum(c$real_scores_fixed[tre,,drop=FALSE],cell_ui[tre],reorder=TRUE)
realsums_exact<-rowsum(c$real_scores_exact[tre,,drop=FALSE],cell_ui[tre],reorder=TRUE)
stopifnot(identical(as.integer(rownames(realsums_fixed)),seq_len(nrow(units))))
for(context_name in unique(md$context)){
 cc<-which(ctl$context==context_name);stopifnot(length(cc)>1L);ct<-ctl[cc];Lctl<-ct$common_library
 rawctx<-raw_ctl[cc,,drop=FALSE];fixed_lg<-log1p(rawctx/Lctl*1e4)
 context_base<-c$control_fixed_base[cc,,drop=FALSE];base_mu<-colMeans(context_base);base_sd<-sd_columns(context_base)
 context_groups<-unique(units[context==context_name,gem]);fixed_base_bygem<-setNames(lapply(context_groups,function(gm){if(c$anchor_type=='matched_GEM')colMeans(context_base[ct$gem==gm,,drop=FALSE])else base_mu}),context_groups)
 context_roster<-if(study_name=='K562')roster else panels[validation_study==study_name&validation_context==context_name,gene]
 stopifnot(length(context_roster)>0L)
 for(gene_name in context_roster){
  ui<-which(units$gene==gene_name&units$context==context_name);if(!length(ui)){cat('Missing frozen gene-context support',gene_name,context_name,'\n');next}
  uu<-units[ui];weights<-uu$n_target/sum(uu$n_target);ii<-which(md$gene==gene_name&md$context==context_name&!md$is_control)
  target_raw<-c$control_target_raw[cc,gene_name];measured_raw<-all(is.finite(target_raw));in_universe<-gene_name%in%universe
  target_log<-if(measured_raw)log1p(1e4*target_raw/Lctl)else rep(NA_real_,length(cc))
  expression<-if(measured_raw)mean(target_log)else NA_real_;detection<-if(measured_raw)mean(target_raw>0)else NA_real_
  all_expr[[length(all_expr)+1L]]<-data.table(study=study_name,context=context_name,gene=gene_name,target_measured_raw=measured_raw,target_in_universe=in_universe,control_expression=expression,control_detection=detection,n_control=length(cc))
  # Fixed-denominator controls for all real and null readouts; own-member removal is exact per-cell.
  cf<-context_base;changed<-which(vapply(spec,function(s)gene_name%in%c(s$positive,s$negative),logical(1)))
  if(length(changed))for(k in changed){
   pos<-setdiff(spec[[k]]$positive,gene_name);neg<-setdiff(spec[[k]]$negative,gene_name)
   if(k<=21L){cf[,k]<-if(length(pos))rowMeans(fixed_lg[,pos,drop=FALSE])else NA_real_;if(length(neg))cf[,k]<-cf[,k]-rowMeans(fixed_lg[,neg,drop=FALSE])}
   else cf[,k]<-(length(spec[[k]]$positive)*context_base[,k]-target_log)/(length(spec[[k]]$positive)-1L)
  }
  mu_fixed<-base_mu;sd_fixed<-base_sd
  if(length(changed)){mu_fixed[changed]<-colMeans(cf[,changed,drop=FALSE]);sd_fixed[changed]<-sd_columns(cf[,changed,drop=FALSE])}
  mu_bygem_fixed<-lapply(uu$gem,function(gm){a<-fixed_base_bygem[[gm]];if(length(changed))a[changed]<-if(c$anchor_type=='matched_GEM')colMeans(cf[ct$gem==gm,changed,drop=FALSE])else mu_fixed[changed];a})
  muf<-do.call(rbind,mu_bygem_fixed);targetf<-c$unit_fixed_score_sums[ui,,drop=FALSE]/uu$n_target
  contrastf<-sweep(targetf-muf,2L,sd_fixed,'/');delf<-colSums(contrastf*weights);mean_targetf<-colSums(targetf*weights);mean_controlf<-colSums(muf*weights)
  # Exact target-excluded normalization for real readouts only.
  Lminus<-Lctl;if(in_universe){stopifnot(measured_raw);Lminus<-Lctl-target_raw};stopifnot(all(Lminus>0))
  exact_lg<-log1p(rawctx/Lminus*1e4)
  ce<-score_logs(t(exact_lg),genes,real_spec,co,rep(gene_name,length(cc)))
  mu_exact<-colMeans(ce);sd_exact<-sd_columns(ce)
  mue<-do.call(rbind,lapply(uu$gem,function(gm)if(c$anchor_type=='matched_GEM')colMeans(ce[ct$gem==gm,,drop=FALSE])else mu_exact))
  targete<-realsums_exact[ui,,drop=FALSE]/uu$n_target;contraste<-sweep(targete-mue,2L,sd_exact,'/');dele<-colSums(contraste*weights)
  member_raw<-as.matrix(c$raw_real_union[ii,,drop=FALSE]);Ltarget<-md$common_library[ii];if(in_universe)Ltarget<-Ltarget-md$assigned_target_count[ii];stopifnot(all(Ltarget>0));member_lg_exact<-log1p(member_raw/Ltarget*1e4);member_lg_fixed<-log1p(member_raw/md$common_library[ii]*1e4)
  for(mode in c('exact_Lminus_target','fixed_L6780')){
   exact<-mode=='exact_Lminus_target';ss<-if(exact)real_spec else spec;aa<-names(ss);nr<-length(ss)
   delta<-if(exact)dele else delf;sdv<-if(exact)sd_exact else sd_fixed;mt<-if(exact)colSums(targete*weights)else mean_targetf;mc<-if(exact)colSums(mue*weights)else mean_controlf
   ctrlmu<-if(exact)mu_exact else mu_fixed;contrast<-if(exact)contraste else contrastf
   remain<-vapply(ss,function(s)length(setdiff(s$positive,gene_name))+length(setdiff(s$negative,gene_name)),integer(1));size<-vapply(ss,function(s)length(s$positive)+length(s$negative),integer(1));own<-vapply(ss,function(s)gene_name%in%c(s$positive,s$negative),logical(1));valid<-is.finite(delta)&is.finite(sdv)&sdv>0&remain>0
   supported<-uu$n_target>=5L;fp<-if(any(supported))colMeans(contrast[supported,,drop=FALSE]>0)else rep(NA_real_,nr)
   all_effects[[length(all_effects)+1L]]<-data.table(study=study_name,context=context_name,gene=gene_name,axis=aa,mode=mode,delta_z=ifelse(valid,delta,NA_real_),mean_target=mt,mean_control=mc,global_control_mean=ctrlmu,control_sd=sdv,n_target=sum(uu$n_target),n_control=length(cc),n_gems=nrow(uu),supported_groups_ge5=sum(supported),f_positive=fp,measured_members=size,remaining_members=remain,few_transcript_proxy=remain<3L,target_member=own,target_measured=measured_raw,target_in_universe=in_universe,control_expression=expression,control_detection=detection,estimable=valid,anchor_type=c$anchor_type,technical_group_type=if(study_name%in%c('K562','RPE1'))'GEM'else if(study_name=='Adamson')'capture_lane'else'pooled_context',genuine_GEM_metadata=study_name%in%c('K562','RPE1'))
   # Real per-group contrasts and target-specific global/group control moments.
   nreal<-21L;cs<-contrast[,seq_len(nreal),drop=FALSE];cm<-if(exact)mue else muf[,seq_len(nreal),drop=FALSE]
   ug<-data.table(gene=gene_name,context=context_name,gem=rep(uu$gem,each=nreal),axis=rep(real_axes,times=nrow(uu)),n_target=rep(uu$n_target,each=nreal),n_control=rep(if(c$anchor_type=='matched_GEM')vapply(uu$gem,function(gm)sum(ct$gem==gm),integer(1))else rep(length(cc),nrow(uu)),each=nreal),contrast_z=as.vector(t(cs)),matched_control_mean=as.vector(t(cm)))
   ug[,`:=`(study=study_name,mode=mode,control_sd=rep(sdv[seq_len(nreal)],times=nrow(uu)),genuine_GEM_metadata=study_name%in%c('K562','RPE1'),technical_group_type=if(study_name%in%c('K562','RPE1'))'GEM'else if(study_name=='Adamson')'capture_lane'else'pooled_context')]
   all_gems[[length(all_gems)+1L]]<-ug
   all_controls[[length(all_controls)+1L]]<-data.table(study=study_name,context=context_name,gene=gene_name,axis=real_axes,mode=mode,control_mean=ctrlmu[seq_len(nreal)],control_sd=sdv[seq_len(nreal)],n_control=length(cc))
   lgc<-if(exact)exact_lg else fixed_lg;lgt<-if(exact)member_lg_exact else member_lg_fixed
   member_cm<-do.call(rbind,lapply(uu$gem,function(gm)if(c$anchor_type=='matched_GEM')colMeans(lgc[ct$gem==gm,,drop=FALSE])else colMeans(lgc)))
   member_mt<-colMeans(lgt);member_mc<-colSums(member_cm*weights);member_sd<-sd_columns(lgc)
   all_members[[length(all_members)+1L]]<-data.table(study=study_name,context=context_name,gene=gene_name,member_gene=genes,mode=mode,mean_target_logCP10k=member_mt,mean_matched_control_logCP10k=member_mc,control_sd=member_sd,delta_logCP10k=member_mt-member_mc,delta_z=ifelse(member_sd>0,(member_mt-member_mc)/member_sd,NA_real_),own_target_member=genes==gene_name,n_target=sum(uu$n_target),n_control=length(cc))
  }
  record<-record+1L;if(record%%50L==0L){cat(study_name,'genecontexts',record,'elapsed',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),'s\n');flush.console()}
 }
}
effects<-rbindlist(all_effects);gems<-rbindlist(all_gems);members<-rbindlist(all_members);controls<-rbindlist(all_controls);expr<-rbindlist(all_expr)
fwrite(effects,file.path(dest,paste0(study_name,'_gene_effects.tsv.gz')),sep='\t');fwrite(gems,file.path(dest,paste0(study_name,'_GEM_effects.tsv.gz')),sep='\t');fwrite(members,file.path(dest,paste0(study_name,'_member_effects.tsv.gz')),sep='\t');fwrite(controls,file.path(dest,paste0(study_name,'_real_control_moments.tsv.gz')),sep='\t');fwrite(expr,file.path(dest,paste0(study_name,'_target_control_expression.tsv')),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(dest,paste0(study_name,'_summary_sessionInfo.txt')))
cat('Completed',study_name,'effectrows',nrow(effects),'GEMrows',nrow(gems),'memberrows',nrow(members),'elapsed',as.numeric(difftime(Sys.time(),started,units='secs')),'s\n')
