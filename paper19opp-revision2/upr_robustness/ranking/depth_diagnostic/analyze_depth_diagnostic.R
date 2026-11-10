#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(if(length(args))args[1]else'paper19opp-revision2');out<-file.path(project,'upr_robustness','ranking','depth_diagnostic');options(stringsAsFactors=FALSE,warn=1);set.seed(2026110610);B<-2000L
sha<-function(f)strsplit(system2('shasum',c('-a','256',shQuote(f)),stdout=TRUE),'[[:space:]]+')[[1]][1]
stopifnot(sha(file.path(out,'locked_plan.txt'))=='bbe5f4a0c63ecce6b08f0f04aa6e1c9a901cedab9699eb8148a81b5dfbf79028')
meta_path<-file.path(project,'reconstruction','cross_dataset','phase2_upr_rpe1','source_rank_inputs.tsv');bg_path<-file.path(project,'reconstruction','cross_dataset','phase2_upr_rpe1','frozen_background_matches.tsv');state_path<-file.path(project,'reconstruction','upr_rpe1_validation','phase2_statistical_state.rds');depth_path<-file.path(project,'upr_robustness','specificity','sequencing_depth_covariates.tsv')
files<-c(meta_path,bg_path,state_path,depth_path);stopifnot(all(file.exists(files)))
write.table(data.frame(file=sub(paste0('^',project,'/'),'',files),sha256=vapply(files,sha,'')),file.path(out,'input_hashes.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(paste('First depth-table read UTC',format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'first_read.txt'));cat('Locked two-covariate diagnostic started',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
meta<-read.delim(meta_path,check.names=FALSE);ids<-sort(meta$gene);stopifnot(length(ids)==438L,!anyDuplicated(ids))
s<-readRDS(state_path);Y<-s$primary$point_validation[ids];stopifnot(all(is.finite(Y)))
depth<-read.delim(depth_path,check.names=FALSE)
required<-c('study','context','gene','n_target_cells','mean_log1p_recorded_UMI','mean_log1p_common_library','n_valid_UMI','n_missing_UMI','n_invalid_UMI','n_valid_common_library','n_missing_common_library','n_invalid_common_library');stopifnot(all(required%in%names(depth)))
k<-subset(depth,study=='K562'&context=='K562'&gene%in%ids);r<-subset(depth,study=='RPE1'&context=='RPE1'&gene%in%ids);stopifnot(!anyDuplicated(k$gene),!anyDuplicated(r$gene))
k<-k[match(ids,k$gene),,drop=FALSE];r<-r[match(ids,r$gene),,drop=FALSE]
source_n<-setNames(meta$source_cells,meta$gene)[ids];rpe_n<-setNames(meta$rpe_cells,meta$gene)[ids]
count_ok_source<-!is.na(k$n_target_cells)&k$n_target_cells==source_n;count_ok_RPE1<-!is.na(r$n_target_cells)&r$n_target_cells==rpe_n
stopifnot(all(count_ok_source),all(count_ok_RPE1))
coverage<-data.frame(gene=ids,original_Hallmark90_Y=as.numeric(Y),source_mean_log1p_recorded_UMI=k$mean_log1p_recorded_UMI,RPE1_mean_log1p_recorded_UMI=r$mean_log1p_recorded_UMI,source_mean_log1p_common_library=k$mean_log1p_common_library,RPE1_mean_log1p_common_library=r$mean_log1p_common_library,source_target_cells=k$n_target_cells,RPE1_target_cells=r$n_target_cells,source_valid_UMI=k$n_valid_UMI,RPE1_valid_UMI=r$n_valid_UMI,source_missing_UMI=k$n_missing_UMI,RPE1_missing_UMI=r$n_missing_UMI,source_invalid_UMI=k$n_invalid_UMI,RPE1_invalid_UMI=r$n_invalid_UMI)
coverage$fit_complete_case<-is.finite(coverage$original_Hallmark90_Y)&is.finite(coverage$source_mean_log1p_recorded_UMI)&is.finite(coverage$RPE1_mean_log1p_recorded_UMI)
write.table(coverage,file.path(out,'gene_depth_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
keep<-which(coverage$fit_complete_case);g<-coverage$gene[keep];y<-coverage$original_Hallmark90_Y[keep];names(y)<-g;X<-cbind(intercept=1,sourceMeanLog1pRecordedUMI=coverage$source_mean_log1p_recorded_UMI[keep],RPE1MeanLog1pRecordedUMI=coverage$RPE1_mean_log1p_recorded_UMI[keep]);rownames(X)<-g
fit<-lm.fit(X,y);stopifnot(fit$rank==3L,all(is.finite(fit$coefficients)));residual<-setNames(as.numeric(y-X%*%fit$coefficients),g);stopifnot(abs(mean(residual))<1e-10)
bg<-read.delim(bg_path,check.names=FALSE);bg<-subset(bg,list=='UPR');stopifnot(nrow(bg)==100L,length(unique(bg$selected_gene))==20L,!anyDuplicated(bg$background_gene))
# Preserve source role order; Top10 uses the original nested source roles.
sets<-lapply(c(10L,20L),function(K){b<-bg;if(K==10L)b<-b[b$Top10,,drop=FALSE];list(K=K,selected=unique(b$selected_gene),background=b$background_gene)});names(sets)<-c('Top10','Top20')
weight_mean<-function(z,w)if(sum(w)>0)sum(z*w)/sum(w)else NA_real_
gap<-function(z,weights,A,B){if(!all(c(A,B)%in%names(z)))return(NA_real_);weight_mean(z[A],weights[A])-weight_mean(z[B],weights[B])}
W<-t(rmultinom(B,length(g),rep(1/length(g),length(g))));colnames(W)<-g;coef_draws<-matrix(NA_real_,B,3L,dimnames=list(NULL,colnames(X)));gap_draws<-matrix(NA_real_,B,2L,dimnames=list(NULL,names(sets)));draw_rank<-rep(NA_integer_,B)
for(b in seq_len(B)){
 w<-W[b,];names(w)<-g;fb<-lm.wfit(X,y,w=w);draw_rank[b]<-fb$rank;if(fb$rank!=3L||any(!is.finite(fb$coefficients)))next;coef_draws[b,]<-fb$coefficients;rb<-setNames(as.numeric(y-X%*%fb$coefficients),g)
 for(j in seq_along(sets))gap_draws[b,j]<-gap(rb,w,sets[[j]]$selected,sets[[j]]$background)
}
ci<-function(z){z<-z[is.finite(z)];if(!length(z))return(c(NA_real_,NA_real_));as.numeric(quantile(z,c(.025,.975)))}
coefficients<-data.frame(term=colnames(X),estimate=as.numeric(fit$coefficients),gene_bootstrap_low=apply(coef_draws,2,function(z)ci(z)[1]),gene_bootstrap_high=apply(coef_draws,2,function(z)ci(z)[2]),complete_case_genes=length(g),original_panel_genes=438,valid_draws=colSums(is.finite(coef_draws)),bootstrap_draws=B,model='original_Hallmark90_Y ~ intercept + sourceMeanLog1pRecordedUMI + RPE1MeanLog1pRecordedUMI')
write.table(coefficients,file.path(out,'depth_model_coefficients.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
gaps<-lapply(seq_along(sets),function(j){a<-sets[[j]];valid<-all(c(a$selected,a$background)%in%g);missing<-setdiff(c(a$selected,a$background),g);raw<-mean(Y[a$selected])-mean(Y[a$background]);expected<-if(a$K==10).573095101503677 else .32876313033789;stopifnot(abs(raw-expected)<1e-10);data.frame(K=a$K,selected_genes=length(a$selected),background_genes=length(a$background),complete_case_genes=length(g),literal_list_estimable=valid,missing_genes=paste(missing,collapse=';'),original_unadjusted_gap=raw,depth_adjusted_continuous_gap=if(valid)mean(residual[a$selected])-mean(residual[a$background])else NA_real_,gene_bootstrap_low=ci(gap_draws[,j])[1],gene_bootstrap_high=ci(gap_draws[,j])[2],valid_draws=sum(is.finite(gap_draws[,j])),bootstrap_draws=B)})
write.table(do.call(rbind,gaps),file.path(out,'depth_adjusted_fixed_matched_gaps.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
# Balance uses all original literal IDs, with available-gene counts explicit for every covariate.
covs<-c('source_mean_log1p_recorded_UMI','RPE1_mean_log1p_recorded_UMI','source_mean_log1p_common_library','RPE1_mean_log1p_common_library');balance<-list();bi<-0L
for(a in sets)for(cv in covs){v<-setNames(coverage[[cv]],coverage$gene);sa<-v[a$selected];sb<-v[a$background];sa<-sa[is.finite(sa)];sb<-sb[is.finite(sb)];pool_sd<-sd(v[is.finite(v)]);bi<-bi+1L;balance[[bi]]<-data.frame(K=a$K,covariate=cv,selected_list_genes=length(a$selected),background_list_genes=length(a$background),selected_available=length(sa),background_available=length(sb),selected_mean=if(length(sa))mean(sa)else NA_real_,background_mean=if(length(sb))mean(sb)else NA_real_,mean_difference=if(length(sa)&&length(sb))mean(sa)-mean(sb)else NA_real_,all438_SD=pool_sd,SMD=if(length(sa)&&length(sb)&&is.finite(pool_sd)&&pool_sd>0)(mean(sa)-mean(sb))/pool_sd else NA_real_,selected_median=if(length(sa))median(sa)else NA_real_,background_median=if(length(sb))median(sb)else NA_real_,selected_SD=if(length(sa)>1)sd(sa)else NA_real_,background_SD=if(length(sb)>1)sd(sb)else NA_real_)}
write.table(do.call(rbind,balance),file.path(out,'fixed_list_depth_balance.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
coverage$depth_model_residual<-NA_real_;coverage$depth_model_residual[keep]<-residual
write.table(coverage,file.path(out,'gene_depth_diagnostic_residuals.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
draw_con<-gzfile(file.path(out,'gene_bootstrap_refit_draws.tsv.gz'),'wt');write.table(data.frame(draw=seq_len(B),model_rank=draw_rank,valid_model=draw_rank==3L,coef_draws,gap_draws,check.names=FALSE),draw_con,sep='\t',quote=FALSE,row.names=FALSE,na='NA');close(draw_con)
saveRDS(W,file.path(out,'gene_bootstrap_joint_weights.rds'),compress='gzip')
write.table(data.frame(check=c('source_cell_counts_match_frozen_metadata','RPE1_cell_counts_match_frozen_metadata','unadjusted_Top10_gap_matches_original','unadjusted_Top20_gap_matches_original','OLS_intercept_residual_mean_zero','single_prespecified_model_no_P_no_reranking'),passed=TRUE),file.path(out,'analysis_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));writeLines(c('One retrospective depth diagnostic; original raw Precision/NDCG and source lists are unchanged.','Outcome is the original Hallmark90 fixed_L6780 signed RPE1 contrast, not an exact-denominator or94-member substitute.','Only the two recorded-UMI covariates enter the nuisance model; common-library depths are descriptive balance variables.','All intervals resample whole genes and refit the same nuisance model; no P values or new inferential family.','Post-treatment recorded UMI combines measurement depth and biological RNA-content variation; adjusted residuals do not estimate causal or functional benefit.'),file.path(out,'interpretation_scope.txt'))
cat('Depth diagnostic complete',format(Sys.time(),tz='UTC',usetz=TRUE),'completecases',length(g),'of438; valid modeldraws',sum(draw_rank==3L),'of',B,'\n');print(do.call(rbind,gaps));print(coefficients)
