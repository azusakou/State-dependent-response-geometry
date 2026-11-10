#!/usr/bin/env Rscript
# Evaluate pre-frozen discovery rankings; validation responses never alter score or targets.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);out<-if(length(a))a[1]else'paper19opp-revision2/reconstruction/candidate_validation';B<-2000L;P<-10000L;set.seed(20261106)
lock<-fread(file.path(out,'harmonized_lock_hashes.tsv'));stopifnot(all(vapply(file.path(out,lock$file),digest,character(1),file=TRUE,algo='sha256')==lock$sha256))
rank<-fread(file.path(out,'harmonized_discovery_ranking.tsv'))
integrated_path<-file.path(out,'secondary_integrated_ranking.tsv');if(file.exists(integrated_path))rank<-merge(rank,fread(integrated_path)[,.(gene,integrated_score,pilot_only_score)],by='gene',all.x=TRUE)
depmap_path<-file.path(out,'discovery_generic_essentiality_baseline.tsv');if(file.exists(depmap_path))rank<-merge(rank,fread(depmap_path)[,.(gene,depmap_generic_essentiality)],by='gene',all.x=TRUE)
axes<-c('apoptosis','proliferation','stress');cols<-paste0('score_',axes)
methods<-c(priority_score='Composite priority',primary_magnitude='Dominant-program magnitude',total_response_magnitude='Total response magnitude',median_replicate_cosine='GEM agreement',discovery_control_expression='Control target expression',discovery_cells='Discovery cell count')
if('integrated_score'%in%names(rank))methods<-c(methods,integrated_score='Source-held-out integrated priority',pilot_only_score='Pilot-only priority')
if('depmap_generic_essentiality'%in%names(rank))methods<-c(methods,depmap_generic_essentiality='Generic essentiality')
quant<-function(x)if(sum(is.finite(x))>=100)as.numeric(quantile(x[is.finite(x)],c(.025,.975)))else c(NA_real_,NA_real_)
partialcor<-function(d,score,outcome,extra_magnitude=FALSE){
 co<-c('discovery_control_expression','discovery_cells','validation_cells','validation_control_expression');cc<-copy(d);for(v in co){cc[,(v):=log1p(pmax(get(v),0))];if(any(!is.finite(cc[[v]])))cc[!is.finite(get(v)),(v):=if(any(is.finite(get(v))))median(get(v)[is.finite(get(v))])else 0]}
 X<-model.matrix(~factor(target_program)+factor(target_sign)+discovery_control_expression+discovery_cells+validation_cells+validation_control_expression,data=cc)
 if(extra_magnitude)X<-cbind(X,magnitude_rank=rank(cc$primary_magnitude,ties.method='average'))
 if(nrow(d)<=qr(X)$rank+3)return(NA_real_)
 xx<-lm.fit(X,rank(cc[[score]],ties.method='average'))$residuals;yy<-lm.fit(X,rank(cc[[outcome]],ties.method='average'))$residuals
 if(sd(xx)<1e-12||sd(yy)<1e-12)NA_real_ else cor(xx,yy)
}
metric_rows<-list();gene_rows<-list();choice_rows<-list();equal_rows<-list();contrast_rows<-list();null_rows<-list()
process<-function(cells,expr,study,control_field,min_cells,base_genes=NULL,context_col='context'){
 setDT(cells);setDT(expr)
 for(ctx in unique(cells[[context_col]])){
  z<-cells[get(context_col)==ctx];ctrl<-z[get(control_field)==TRUE&is.finite(score_apoptosis)&is.finite(score_proliferation)&is.finite(score_stress)]
  if(nrow(ctrl)<min_cells)next
  means<-vapply(cols,function(v)mean(ctrl[[v]]),numeric(1));sds<-vapply(cols,function(v)sd(ctrl[[v]]),numeric(1));sds<-pmax(sds,1e-8)
  dd<-z[!get(control_field),c(list(validation_cells=.N),lapply(.SD,mean)),by=gene,.SDcols=cols]
  dd<-merge(rank,dd,by='gene');dd<-dd[metadata_eligible==TRUE&validation_cells>=min_cells]
  if(!is.null(base_genes))dd<-dd[gene%in%base_genes]
  ex<-expr;if('context'%in%names(expr))ex<-expr[context==ctx];if('control_mean_logcp10k'%in%names(ex))setnames(ex,'control_mean_logcp10k','validation_control_expression')
  if(!'validation_control_expression'%in%names(ex))stop('Missing validation baseline expression')
  dd<-merge(dd,ex[,.(gene,validation_control_expression)],by='gene',all.x=TRUE)
  vv<-sweep(sweep(as.matrix(dd[,..cols]),2,means,'-'),2,sds,'/');tt<-match(dd$target_program,axes)
  dd[,directional_effect:=target_sign*vv[cbind(seq_len(.N),tt)]]
  dd[,validation_off_target_norm:=sqrt(rowSums(vv^2)-vv[cbind(seq_len(.N),tt)]^2)]
  dd[,specific_directional_effect:=directional_effect/(1+validation_off_target_norm)]
  dd[,`:=`(study=study,context=ctx,control_cells=nrow(ctrl),validation_objective=paste(target_program,target_sign,sep='_'))]
  gene_rows[[length(gene_rows)+1L]]<<-dd
  n<-nrow(dd);if(n<12)next
  # Gene bootstrap conditions on the pre-frozen rankings and assayed population.
  bootids<-replicate(B,sample.int(n,n,replace=TRUE),simplify=FALSE)
  for(outcome in c('directional_effect','specific_directional_effect')){
   y<-dd[[outcome]]
   for(m in names(methods)){
    if(sum(is.finite(dd[[m]]))<n)next
    rho<-cor(dd[[m]],y,method='spearman');rb<-vapply(bootids,function(ii)suppressWarnings(cor(dd[[m]][ii],y[ii],method='spearman')),numeric(1));qq<-quant(rb)
    pc<-partialcor(dd,m,outcome);pcmag<-partialcor(dd,m,outcome,TRUE)
    metric_rows[[length(metric_rows)+1L]]<<-data.table(study=study,context=ctx,outcome=outcome,method=methods[[m]],genes=n,spearman=rho,conditional_gene_bootstrap_low=qq[1],conditional_gene_bootstrap_high=qq[2],partial_rank_correlation=pc,partial_rank_correlation_plus_magnitude=pcmag)
   }
   for(k in c(10L,20L)){
    if(n<=k)next
    groups<-split(seq_len(n),dd$validation_objective);selections<-list()
    for(m in names(methods)){
     if(sum(is.finite(dd[[m]]))<n)next
     ix<-order(-dd[[m]],dd$gene,na.last=TRUE)[seq_len(k)];selected<-seq_len(n)%in%ix;selections[[m]]<-selected
     gap<-mean(y[selected])-mean(y[!selected]);gb<-vapply(bootids,function(ii){a<-selected[ii];if(!any(a)||all(a))return(NA_real_);mean(y[ii][a])-mean(y[ii][!a])},numeric(1));qg<-quant(gb)
     choice_rows[[length(choice_rows)+1L]]<<-data.table(study=study,context=ctx,outcome=outcome,method=methods[[m]],top_k=k,genes=n,top_mean=mean(y[selected]),rest_mean=mean(y[!selected]),gap=gap,conditional_gene_bootstrap_low=qg[1],conditional_gene_bootstrap_high=qg[2],direction_concordance=mean(dd$directional_effect[selected]>0),selected_genes=paste(dd$gene[ix],collapse=';'))
    }
    sel<-selections[['priority_score']];bm<-selections[['primary_magnitude']]
    delta<-mean(y[sel])-mean(y[bm]);bb<-vapply(bootids,function(ii){aa<-sel[ii];mm<-bm[ii];if(!any(aa)||!any(mm))return(NA_real_);mean(y[ii][aa])-mean(y[ii][mm])},numeric(1));qq<-quant(bb)
    contrast_rows[[length(contrast_rows)+1L]]<<-data.table(study=study,context=ctx,outcome=outcome,top_k=k,composite_minus_magnitude_top_mean=delta,conditional_gene_bootstrap_low=qq[1],conditional_gene_bootstrap_high=qq[2],selection_overlap=sum(sel&bm))
    if('integrated_score'%in%names(selections)){ss<-selections[['integrated_score']];for(refname in c('priority_score','primary_magnitude','pilot_only_score')){rr<-selections[[refname]];dv<-mean(y[ss])-mean(y[rr]);vb<-vapply(bootids,function(ii){u<-ss[ii];v<-rr[ii];if(!any(u)||!any(v))return(NA_real_);mean(y[ii][u])-mean(y[ii][v])},numeric(1));qc<-quant(vb);contrast_rows[[length(contrast_rows)+1L]]<<-data.table(study=study,context=ctx,outcome=outcome,top_k=k,comparison=paste('integrated_minus',refname),integrated_top_mean_difference=dv,conditional_gene_bootstrap_low=qc[1],conditional_gene_bootstrap_high=qc[2],selection_overlap=sum(ss&rr))}}
    # Random subsets preserve selected objective counts and measured ascertainment balance.
    cov<-cbind(log1p(dd$discovery_cells),log1p(dd$validation_cells),log1p(pmax(dd$discovery_control_expression,0)),log1p(pmax(dd$validation_control_expression,0)))
    for(j in seq_len(ncol(cov))){cov[!is.finite(cov[,j]),j]<-if(any(is.finite(cov[,j])))median(cov[is.finite(cov[,j]),j])else 0;ss<-sd(cov[,j]);cov[,j]<-if(ss>0)(cov[,j]-mean(cov[,j]))/ss else 0}
    target_balance<-colMeans(cov[sel,,drop=FALSE]);ns<-table(dd$validation_objective[sel]);null<-numeric();maxsmd<-numeric()
    for(b in seq_len(P)){
     jj<-unlist(lapply(names(ns),function(g)groups[[g]][sample.int(length(groups[[g]]),as.integer(ns[[g]]),replace=FALSE)]),use.names=FALSE)
     balance<-max(abs(colMeans(cov[jj,,drop=FALSE])-target_balance));if(balance<=.5){null<-c(null,mean(y[jj]));maxsmd<-c(maxsmd,balance)}
    }
    null_rows[[length(null_rows)+1L]]<<-data.table(study=study,context=ctx,outcome=outcome,top_k=k,proposals=P,accepted_draws=length(null),balance_limit_SD=.5,selected_mean=mean(y[sel]),matched_random_mean=if(length(null))mean(null)else NA_real_,matched_random_95_low=if(length(null))quantile(null,.025)else NA_real_,matched_random_95_high=if(length(null))quantile(null,.975)else NA_real_,exploratory_upper_tail_p=if(length(null)>=100)(1+sum(null>=mean(y[sel])))/(1+length(null))else NA_real_)
   }
  }
  # A fixed equal-cell count is applied to every accepted target and the controls.
  cellix<-split(seq_len(nrow(z)),z$gene);ctix<-which(z[[control_field]])
  for(b in seq_len(200L)){
   cb<-sample(ctix,min_cells,replace=FALSE);cm<-vapply(cols,function(v)mean(z[[v]][cb]),numeric(1));cs<-vapply(cols,function(v)sd(z[[v]][cb]),numeric(1));cs<-pmax(cs,1e-8)
   yy<-vapply(seq_len(n),function(i){ii<-sample(cellix[[dd$gene[i]]],min_cells,replace=FALSE);v<-(vapply(cols,function(v)mean(z[[v]][ii]),numeric(1))-cm)/cs;dd$target_sign[i]*v[tt[i]]},numeric(1))
   ord<-order(-dd$priority_score,dd$gene);k<-min(10L,n-1L);equal_rows[[length(equal_rows)+1L]]<<-data.table(study=study,context=ctx,draw=b,cells_per_arm=min_cells,genes=n,spearman=cor(dd$priority_score,yy,method='spearman'),top10_mean=mean(yy[ord[seq_len(k)]]),rest_mean=mean(yy[ord[-seq_len(k)]]),gap=mean(yy[ord[seq_len(k)]])-mean(yy[ord[-seq_len(k)]]))
  }
 }
}
# Required primary external tables are produced only after the harmonized rank lock.
adam<-fread(file.path(out,'adamson_cell_scores_harmonized.tsv.gz'));ae<-fread(file.path(out,'adamson_control_target_expression_harmonized.tsv'))
process(adam,ae,'Adamson2016', 'is_control',100L)
covered<-if('integrated_score'%in%names(rank))rank[is.finite(integrated_score)&is.finite(pilot_only_score),gene]else character();if(length(covered))process(adam,ae,'Adamson2016 matched model coverage','is_control',100L,covered)
frangpath<-file.path(out,'frangieh_cell_scores_harmonized.tsv.gz')
if(file.exists(frangpath)){fr<-fread(frangpath);fe<-fread(file.path(out,'frangieh_control_target_expression_harmonized.tsv'));primary_fr<-fread(file.path(out,'harmonized_frangieh_candidates.tsv'))$gene;process(fr,fe,'Frangieh2021','is_control',30L,primary_fr);if(length(covered))process(fr,fe,'Frangieh2021 matched model coverage','is_control',30L,intersect(primary_fr,covered))}
for(z in c('metric_rows','gene_rows','choice_rows','equal_rows','contrast_rows','null_rows'))fwrite(rbindlist(get(z),fill=TRUE),file.path(out,paste0(sub('_rows$','',z),'.tsv.gz')),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'analysis_sessionInfo.txt'))
