#!/usr/bin/env Rscript
# Descriptive post-outcome sensitivity of the fixed readout; no new rank or branch.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
o<-'paper19opp-revision2/reconstruction/cross_dataset/phase2_upr_rpe1';out<-file.path(o,'readout_member_sensitivity');dir.create(out,showWarnings=FALSE)
upr<-fread(file.path(o,'phase2_UPR_members.tsv'))$gene
matches<-fread(file.path(o,'frozen_background_matches.tsv'))[list=='UPR'];selected<-unique(matches$selected_gene);background<-unique(matches$background_gene);used<-c(selected,background)
stopifnot(length(upr)==90L,length(selected)==20L,length(background)==100L,!anyDuplicated(used))
all_gene<-list();all_drop<-list();all_member<-list();all_control<-list()
for(mode in c('source','rpe1')){
 if(mode=='source'){ca<-readRDS(file.path(o,'original91_source_raw_cache.rds'));cells<-as.data.table(ca$cells);lib<-cells$phase2_library;raw<-ca$raw_upr[,upr,drop=FALSE];label<-'K562'}
 else{ca<-readRDS(file.path(o,'RPE1_phase2_cache.rds'));cells<-as.data.table(ca$cell_metadata);lib<-cells$common_library;raw<-ca$cell_UPR_raw_counts[,upr,drop=FALSE];label<-'RPE1'}
 logcounts<-log1p(sweep(as.matrix(raw),1,lib,'/')*1e4)
 ci<-which(cells$is_control);C<-logcounts[ci,,drop=FALSE];cgem<-cells$gem_group[ci];gems<-sort(unique(cgem));cgindex<-match(cgem,gems);cn<-tabulate(cgindex,length(gems));cmeans<-rowsum(C,cgindex,reorder=TRUE)/cn
 csd<-apply(C,2,sd);sumC<-rowSums(C)
 all_control[[mode]]<-data.table(context=label,readout_gene=upr,control_cells=length(ci),mean_control_logCP10k=colMeans(C),sd_control_logCP10k=csd,control_detected_fraction=colMeans(C>0))
 for(g in used){
  ii<-which(!cells$is_control&cells$gene==g);stopifnot(length(ii)>0);T<-logcounts[ii,,drop=FALSE];tg<-match(cells$gem_group[ii],gems);stopifnot(!anyNA(tg));tw<-tabulate(tg,length(gems))/length(ii)
  weightedControl<-as.numeric(crossprod(tw,cmeans));delta<-colMeans(T)-weightedControl
  if(g%in%selected){all_member[[paste(mode,g)]]<-data.table(context=label,assayed_gene=g,readout_gene=upr,treated_cells=length(ii),control_cells=length(ci),target_GEM_labels=uniqueN(cells$gem_group[ii]),mean_treated_logCP10k=colMeans(T),mean_GEM_matched_control_logCP10k=weightedControl,control_global_SD_logCP10k=csd,delta_logCP10k=delta,deltaZ_readout=ifelse(csd>0,delta/csd,NA_real_),readout_is_own_target=upr==g,included_in_primary_UPR=upr!=g)}
  member<-g%in%upr;mg<-90L-as.integer(member);sumT<-rowSums(T);ownC<-if(member)C[,g]else 0;ownT<-if(member)T[,g]else 0
  qC<-(sumC-ownC)/mg;qT<-(sumT-ownT)/mg;cm0<-as.numeric(rowsum(qC,cgindex,reorder=TRUE))/cn;dz0<-(mean(qT)-sum(tw*cm0))/sd(qC)
  DC<-sweep(-C,1,sumC-ownC,'+')/(mg-1L);DT<-sweep(-T,1,sumT-ownT,'+')/(mg-1L)
  if(member){DC[,g]<-qC;DT[,g]<-qT}
  cm<-rowsum(DC,cgindex,reorder=TRUE)/cn;dd<-colMeans(DT)-as.numeric(crossprod(tw,cm));ss<-apply(DC,2,sd);dz<-dd/ss;stopifnot(all(is.finite(dz)))
  all_gene[[paste(mode,g)]]<-data.table(context=label,gene=g,is_selected=g%in%selected,primary_deltaZ=dz0,treated_cells=length(ii))
  all_drop[[paste(mode,g)]]<-data.table(context=label,gene=g,is_selected=g%in%selected,removed_readout_gene=upr,deltaZ=dz,already_own_excluded=upr==g)
 }
 rm(ca,cells,raw,logcounts,C);gc()
}
g<-rbindlist(all_gene);d<-rbindlist(all_drop);members<-rbindlist(all_member);controls<-rbindlist(all_control)
gaps<-d[,.(Top20_mean=mean(deltaZ[is_selected]),background_mean=mean(deltaZ[!is_selected]),gap=mean(deltaZ[is_selected])-mean(deltaZ[!is_selected])),by=.(context,removed_readout_gene)]
base<-g[,.(primary_gap=mean(primary_deltaZ[is_selected])-mean(primary_deltaZ[!is_selected])),by=context]
summary<-merge(gaps[,.(readout_members=.N,minimum_leave_one_gap=min(gap),maximum_leave_one_gap=max(gap),minimum_gap_member=removed_readout_gene[which.min(gap)],maximum_gap_member=removed_readout_gene[which.max(gap)],positive_leave_one_gaps=sum(gap>0)),by=context],base,by='context')
fwrite(members,file.path(out,'fixed_Top20_UPR90_member_contrasts.tsv.gz'),sep='\t');fwrite(controls,file.path(out,'UPR90_control_member_expression.tsv'),sep='\t');fwrite(g,file.path(out,'fixed_Top20_background_gene_effects.tsv'),sep='\t');fwrite(d,file.path(out,'all_leave_one_member_gene_effects.tsv.gz'),sep='\t');fwrite(gaps,file.path(out,'fixed_Top20_leave_one_member_gaps.tsv'),sep='\t');fwrite(summary,file.path(out,'leave_one_member_gap_summary.tsv'),sep='\t')
writeLines(c('Descriptive post-outcome readout sensitivity, performed only after the phase2 primary analysis. The locked Top20 genes and100 frozen matched background genes are unchanged. No re-ranking, alternate branch scoring or readout-member selection is used.',
'Per-member contrasts are target-cell-count-weighted matched-GEM logCP10k changes under fixed L6780, with global NTC gene SD supplied separately. Own-target transcript rows are retained and explicitly marked as excluded from the original program; they are not imputed as zero.',
'For each of all90UPR members, the program excludes that readout member and the assayed target when applicable (a shared exclusion is applied only once). Every gene-specific control mean/SD is recomputed. The primary library denominator remains fixed. Matched Top20-minus-background gaps are recomputed without changing any list or match. All90results are retained.',
'This is a UPR-associated expression readout; the analysis does not identify canonical ER-UPR activation, separate ER-UPR from ISR mechanistically, or establish pathway causality.'),file.path(out,'readout_sensitivity_scope.txt'))
print(summary)
