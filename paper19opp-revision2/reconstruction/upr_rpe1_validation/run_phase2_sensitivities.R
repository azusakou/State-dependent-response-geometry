#!/usr/bin/env Rscript
# Fixed-list descriptive and declared exact-composition checks. No reranking.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
a<-commandArgs(TRUE);out<-if(length(a))a[1]else'paper19opp-revision2/reconstruction/upr_rpe1_validation'
here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))))
source(file.path(here,'phase2_statistics.R'))
z<-readRDS(file.path(out,'phase2_statistical_state.rds'));p<-z$primary;src<-z$source_directory
y<-p$point_validation;x<-p$point_source;sets<-p$sets
summary<-rbindlist(lapply(c(20L,10L),function(k){
 ss<-sets[seq_len(k)];sel<-vapply(ss,`[`,character(1),1L);bg<-lapply(ss,function(v)v[-1]);sv<-y[sel];bv<-vapply(bg,function(g)mean(y[g]),numeric(1));av<-abs(sv);ab<-vapply(bg,function(g)mean(abs(y[g])),numeric(1));ids<-if(k==20L)p$set20_ids else p$set10_ids
 qq<-function(v)phase2_quantile(vapply(ids,function(i)mean(v[i]),numeric(1)))
 qsel<-qq(sv);qbg<-qq(bv);qa<-qq(av);qab<-qq(ab);qag<-qq(av-ab)
 data.table(TopK=k,selected_mean=mean(sv),selected_mean_low=qsel[1],selected_mean_high=qsel[2],background_mean=mean(bv),background_mean_low=qbg[1],background_mean_high=qbg[2],selected_positive=sum(sv>0),selected_negative=sum(sv<0),background_positive=sum(y[unlist(bg)]>0),background_negative=sum(y[unlist(bg)]<0),selected_mean_absolute=mean(av),selected_absolute_low=qa[1],selected_absolute_high=qa[2],background_mean_absolute=mean(ab),background_absolute_low=qab[1],background_absolute_high=qab[2],absolute_matched_gap=mean(av-ab),absolute_gap_low=qag[1],absolute_gap_high=qag[2],interval_unit='resampled intact frozen matched sets; descriptive')
}))
fwrite(summary,file.path(out,'fixed_list_direction_absolute_summaries.tsv'),sep='\t')
gd<-z$source_metadata[,.(gene,source_UPR_deltaZ,source_UPR_f_positive,UPR_priority,source_cells,source_GEM_ge5,rpe_cells,rpe_GEM_ge5)]
gd<-merge(gd,z$ranking,by='gene');gd[,RPE1_UPR_deltaZ:=y[gene]];gd[,source_RPE1_sign_agrees:=sign(source_UPR_deltaZ)==sign(RPE1_UPR_deltaZ)];setorder(gd,rank)
fwrite(gd,file.path(out,'all438_gene_UPR_response_and_fixed_rank.tsv'),sep='\t')
sg<-rbindlist(lapply(list(all438=gd,source_positive=gd[source_UPR_deltaZ>0],Top20=gd[rank<=20],Top10=gd[rank<=10]),function(d)data.table(genes=nrow(d),source_positive=sum(d$source_UPR_deltaZ>0),RPE1_positive=sum(d$RPE1_UPR_deltaZ>0),same_sign=sum(d$source_RPE1_sign_agrees))),idcol='group')
fwrite(sg,file.path(out,'fixed_direction_agreement.tsv'),sep='\t')
loo<-data.table(omitted_gene=names(x),rho=vapply(names(x),function(g)phase2_spearman(x[names(x)!=g],y[names(x)!=g]),numeric(1)))
fwrite(loo,file.path(out,'signed_correlation_leave_one_gene_out.tsv'),sep='\t')
ll<-rbindlist(lapply(c(20L,10L),function(k){ss<-sets[seq_len(k)];gaps<-phase2_set_contrasts(y,ss);data.table(TopK=k,omitted_matched_set=names(ss),matched_gap=vapply(seq_len(k),function(i)mean(gaps[-i]),numeric(1)),selected_mean=vapply(seq_len(k),function(i)mean(y[names(ss)[-i]]),numeric(1)))}))
fwrite(ll,file.path(out,'fixed_list_leave_one_matched_set_out.tsv'),sep='\t')
# Exact L-target sensitivity uses cached raw-count recomputation in both contexts.
sf<-fread(file.path(src,'source_UPR_exact_denominator_gene_effects.tsv'));vf<-fread(file.path(src,'rpe1_UPR_exact_denominator_gene_effects.tsv'))
sx<-copy(z$source);sx[,source_delta_z:=sf$deltaZ[match(gene,sf$gene)]];vv<-vf[,.(gene,validation_delta_z=deltaZ)]
adapt<-function(f){d<-fread(f);d[,.(gene,gem=as.character(gem_group),target_cells=n_target,contrast_z=deltaZ)]}
sgm<-adapt(file.path(src,'source_UPR_exact_denominator_GEM_effects.tsv.gz'));vgm<-adapt(file.path(src,'rpe1_UPR_exact_denominator_GEM_effects.tsv.gz'))
ex<-phase2_primary(sx,vv,z$ranking,z$matches,sgm,vgm,seed=20261110L)
fwrite(ex$endpoints,file.path(out,'exact_denominator_sensitivity_endpoints.tsv'),sep='\t')
ccr<-fread(file.path(src,'complete_covariate_source_rankings.tsv'));ccgenes<-unique(ccr$gene)
base<-phase2_baselines(ex,sx,ccr[,.(method,gene,rank=rank_cc)],ccgenes,'raw_UPR')
fwrite(base,file.path(out,'exact_denominator_baseline_utility.tsv'),sep='\t')
ccnr<-fread(file.path(src,'complete_covariate_null100_rankings.tsv.gz'));ccm<-fread(file.path(src,'complete_covariate_background_matches.tsv'));ey<-ex$point_validation
cal<-rbindlist(lapply(sort(unique(ccnr$module)),function(mm){rr<-ccnr[module==mm,.(gene,rank=rank_cc)];ss<-phase2_make_sets(ccm[list==paste0('CC_',mm),.(selected_gene,background_gene)],rr);data.table(module=mm,matched_gap=mean(phase2_set_contrasts(ey,ss)))}))
fwrite(cal,file.path(out,'exact_denominator_CC100_fixed_null_utility.tsv'),sep='\t')
obs<-ex$endpoints[grepl('Top20',endpoint),estimate];tail<-(1+sum(cal$matched_gap>=obs-1e-12))/101
fwrite(data.table(Top20_real_UPR_exact_denominator_gap=obs,fixed_null_modules=100L,upper_tail_fraction=tail,note='fixed source ranks and matched sets; real UPR exact denominator sensitivity; no null-module correlation claim'),file.path(out,'exact_denominator_CC100_utility_calibration.tsv'),sep='\t')
saveRDS(ex,file.path(out,'exact_denominator_statistical_state.rds'),compress='xz')
print(summary);print(sg);print(ex$endpoints);print(range(loo$rho));print(cal[order(-matched_gap)][1:3]);print(tail)
