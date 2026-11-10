#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);project<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(project,'.Rlib'),.libPaths()));source(file.path(out,'phase2_helpers.R'))
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db)})
e<-fread(file.path(out,'source_phase2_gene_effects.tsv.gz'));eligible<-fread(file.path(out,'phase2_eligible_genes.tsv'))
upr<-e[axis=='UPR'];upr[,UPR_priority:=pmax(deltaZ,0)*f_positive_GEM]
setnames(upr,c('deltaZ','f_positive_GEM'),c('source_UPR_deltaZ','source_UPR_f_positive'))
source_table<-merge(eligible,upr[,.(gene,source_UPR_deltaZ,source_UPR_f_positive,UPR_priority,source_UPR_control_sd=control_sd)],by='gene',all.x=TRUE)
core<-dcast(e[axis%in%paste0('core_',c('apoptosis','proliferation','stress'))],gene~axis,value.var='deltaZ')
core[,generic_three_program_magnitude:=sqrt(core_apoptosis^2+core_proliferation^2+core_stress^2)]
source_table<-merge(source_table,core,by='gene',all.x=TRUE)
ntc<-fread(file.path(out,'source_NTC_gene_statistics.tsv'))
source_table<-merge(source_table,ntc[,.(gene,source_control_expression=source_NTC_mean_logCP10k,source_control_detection=source_NTC_detection)],by='gene',all.x=TRUE)
# Prefer the author's measured-feature Ensembl identity over ambiguous old symbols.
p<-readRDS(file.path(out,'source_preparation.rds'));f<-hdf5r::H5File$new(p$source_path,'r');ens<-as.character(read_annot(f,'var','gene_id'));f$close_all()
author_ids<-unique(data.table(gene=p$raw_gene,ENSEMBL=ens));valid_ens<-intersect(unique(author_ids$ENSEMBL),keys(org.Hs.eg.db,keytype='ENSEMBL'))
ids<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=valid_ens,keytype='ENSEMBL',columns='ENTREZID'))
gene_ids<-merge(author_ids,ids,by='ENSEMBL',allow.cartesian=TRUE)[!is.na(ENTREZID),.(entrez=paste(sort(unique(ENTREZID)),collapse=';'),n_entrez=uniqueN(ENTREZID)),by=gene]
gene_ids[n_entrez!=1L,entrez:=NA_character_]
missing<-setdiff(source_table$gene,gene_ids[!is.na(entrez),gene])
for(keytype in c('SYMBOL','ALIAS')){
 wanted<-intersect(missing,keys(org.Hs.eg.db,keytype=keytype));if(!length(wanted))next
 a<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=wanted,keytype=keytype,columns='ENTREZID'));setnames(a,keytype,'gene')
 z<-a[!is.na(ENTREZID),.(entrez=paste(sort(unique(ENTREZID)),collapse=';'),n_entrez=uniqueN(ENTREZID)),by=gene][n_entrez==1L]
 if(nrow(z)){gene_ids<-rbindlist(list(gene_ids[!gene%in%z$gene],z),fill=TRUE);missing<-setdiff(missing,z$gene)}
}
source_table<-merge(source_table,gene_ids[,.(gene,entrez,n_entrez)],by='gene',all.x=TRUE)
fwrite(source_table[,.(gene,entrez,n_entrez)],file.path(out,'source_gene_identifier_mapping.tsv'),sep='\t')
depmap_file<-'/Volumes/wang/paper4rfv/data/raw/depmap/24Q4/CRISPRGeneEffect.csv'
column_names<-names(fread(depmap_file,nrows=0));col_ids<-sub('^.*\\(([0-9]+)\\)$','\\1',column_names)
select_columns<-which(col_ids%in%source_table$entrez)
dm<-fread(depmap_file,select=select_columns);dep<-data.table(entrez=col_ids[select_columns],depmap_essentiality=-vapply(dm,function(v)mean(v[is.finite(v)]),numeric(1)),depmap_n_lines=vapply(dm,function(v)sum(is.finite(v)),integer(1)))
dep<-dep[,if(.N==1L).SD else NULL,by=entrez]
source_table<-merge(source_table,dep,by='entrez',all.x=TRUE,sort=FALSE)
reference<-readRDS(file.path(project,'paper19opp-revision1/GitHub/reference/published_numerical_reference.rds'))$evidence
annotation<-as.data.table(reference)[,.(gene,genomewide_effect_fdr,genomewide_effect_p_min,genomewide_cells)]
source_table<-merge(source_table,annotation,by='gene',all.x=TRUE)
source_table[,genomewide_annotation_score:=-log10(pmax(genomewide_effect_fdr,.Machine$double.xmin))]
source_table[,source_count_score:=log1p(source_cells)]
source_table[,expression_stratum:=pmin(3L,as.integer((frank(source_control_expression,ties.method='average',na.last='keep')-1)/sum(is.finite(source_control_expression))*3)+1L)]
source_table[,count_stratum:=pmin(3L,as.integer((frank(source_cells,ties.method='average')-1)/.N*3)+1L)]
source_table[,permutation_stratum:=paste(ifelse(is.na(expression_stratum),'unknown',expression_stratum),count_stratum,sep='|')]
stopifnot(nrow(source_table)==438L,all(is.finite(source_table$UPR_priority)))
fwrite(source_table,file.path(out,'source_rank_inputs.tsv'),sep='\t')
methods<-c(UPR_priority='UPR_priority',raw_UPR='source_UPR_deltaZ',generic_magnitude='generic_three_program_magnitude',DepMap='depmap_essentiality',expression='source_control_expression',count='source_count_score',genomewide_annotation='genomewide_annotation_score')
allranks<-list();toplists<-list()
for(m in names(methods)){
 value<-source_table[[methods[[m]]]];z<-data.table(gene=source_table$gene,method=m,score=value,available=is.finite(value),source_UPR_deltaZ=source_table$source_UPR_deltaZ)
 z<-z[order(!available,-score,if(m=='UPR_priority')-source_UPR_deltaZ else gene,gene)]
 z[,rank:=seq_len(.N)];allranks[[m]]<-z
 toplists[[m]]<-z[available==TRUE&rank<=20][,.(gene,method,score,rank,Top10=rank<=10,Top20=TRUE)]
}
ranks<-rbindlist(allranks);tops<-rbindlist(toplists)
fwrite(ranks,file.path(out,'all_source_comparator_rankings.tsv'),sep='\t');fwrite(tops,file.path(out,'source_Top10_Top20_lists.tsv'),sep='\t')
fwrite(ranks[,.(eligible_genes=.N,available_genes=sum(available),missing_genes=paste(gene[!available],collapse=';')),by=method],file.path(out,'source_comparator_coverage.tsv'),sep='\t')
complete_genes<-ranks[,.(complete=all(available)),by=gene][complete==TRUE,gene]
fair<-ranks[gene%in%complete_genes]
fair<-fair[order(method,-score,ifelse(method=='UPR_priority',-source_UPR_deltaZ,0),gene)]
fair[,rank_complete_case:=seq_len(.N),by=method]
fair[,`:=`(cohort='common_complete_comparator_roster',Top10_complete_case=rank_complete_case<=10,Top20_complete_case=rank_complete_case<=20)]
fwrite(fair,file.path(out,'common_complete_case_comparator_rankings.tsv'),sep='\t')
covariates<-c('source_count_score','source_control_expression','generic_three_program_magnitude','depmap_essentiality')
availability<-source_table[,c('gene',covariates),with=FALSE];availability[,complete_matching_covariates:=complete.cases(.SD),.SDcols=covariates]
fwrite(availability,file.path(out,'matching_covariate_availability.tsv'),sep='\t')
match_one<-function(selected,label){
 if(length(selected)!=20L)return(list(matches=NULL,status=data.table(list=label,estimable=FALSE,reason='fewer than20 source-ranked genes',missing_genes=NA_character_)))
 missing<-availability[gene%in%selected&complete_matching_covariates==FALSE,gene]
 if(length(missing))return(list(matches=NULL,status=data.table(list=label,estimable=FALSE,reason='selected member lacks a required matching covariate; no truncation or refill',missing_genes=paste(missing,collapse=';'))))
 eligible_control<-availability[complete_matching_covariates==TRUE&!gene%in%selected,gene]
 if(length(eligible_control)<100L)return(list(matches=NULL,status=data.table(list=label,estimable=FALSE,reason='fewer than100 distinct complete-case background genes',missing_genes=NA_character_)))
 dat<-source_table[match(c(selected,eligible_control),gene)];xx<-as.matrix(dat[,..covariates]);mu<-colMeans(xx);sdv<-apply(xx,2,sd);sdv[sdv<1e-8]<-1;zz<-sweep(sweep(xx,2,mu,'-'),2,sdv,'/');rownames(zz)<-dat$gene
 used<-character();mm<-list()
 for(i in seq_along(selected)){
  gene<-selected[i];pool<-setdiff(eligible_control,used);distance<-rowSums(sweep(zz[pool,,drop=FALSE],2,zz[gene,],'-')^2)
  choices<-data.table(gene=pool,distance=sqrt(distance))[order(distance,gene)][1:5]
  used<-c(used,choices$gene);mm[[i]]<-data.table(list=label,set=i,selected_gene=gene,background_gene=choices$gene,distance=choices$distance,Top10=i<=10)
 }
 list(matches=rbindlist(mm),status=data.table(list=label,estimable=TRUE,reason='all20 members retained; five distinct backgrounds per selected gene',missing_genes=NA_character_))
}
matched<-list();statuses<-list();nullranks<-list()
primary<-tops[method=='UPR_priority'][order(rank),gene];m<-match_one(primary,'UPR');matched$UPR<-m$matches;statuses$UPR<-m$status
for(a in grep('^null_',unique(e$axis),value=TRUE)){
 z<-e[axis==a];z[,score:=pmax(deltaZ,0)*f_positive_GEM];setorder(z,-score,-deltaZ,gene);z[,rank:=.I]
 nullranks[[a]]<-z[,.(gene,module=a,source_deltaZ=deltaZ,f_positive_GEM,score,rank,Top10=rank<=10,Top20=rank<=20)]
 m<-match_one(z[1:20,gene],a);matched[[a]]<-m$matches;statuses[[a]]<-m$status
}
fwrite(rbindlist(nullranks),file.path(out,'null100_source_rankings.tsv.gz'),sep='\t')
status<-rbindlist(statuses);fwrite(status,file.path(out,'matching_estimability_before_RPE1.tsv'),sep='\t')
matches<-rbindlist(matched,fill=TRUE);fwrite(matches,file.path(out,'frozen_background_matches.tsv'),sep='\t')
balance<-list()
for(nm in unique(matches$list)){
 z<-matches[list==nm];for(k in c(10L,20L)){zz<-z[set<=k];sel<-source_table[match(unique(zz$selected_gene),gene)];bg<-source_table[match(zz$background_gene,gene)]
  for(v in covariates){sdpool<-sd(source_table[[v]],na.rm=TRUE);balance[[length(balance)+1L]]<-data.table(list=nm,TopK=k,covariate=v,selected_mean=mean(sel[[v]]),background_mean=mean(bg[[v]]),SMD=(mean(sel[[v]])-mean(bg[[v]]))/sdpool)}
 }
}
fwrite(rbindlist(balance),file.path(out,'matching_balance_before_RPE1.tsv'),sep='\t')
ccgenes<-availability[complete_matching_covariates==TRUE,gene]
fwrite(source_table[gene%in%ccgenes],file.path(out,'complete_covariate_gene_roster.tsv'),sep='\t')
ccr<-ranks[gene%in%ccgenes][order(method,-score,ifelse(method=='UPR_priority',-source_UPR_deltaZ,0),gene)]
ccr[,rank_cc:=seq_len(.N),by=method];ccr[,`:=`(Top10_cc=rank_cc<=10,Top20_cc=rank_cc<=20)]
fwrite(ccr,file.path(out,'complete_covariate_source_rankings.tsv'),sep='\t')
ccnull<-rbindlist(nullranks)[gene%in%ccgenes][order(module,-score,-source_deltaZ,gene)]
ccnull[,rank_cc:=seq_len(.N),by=module];ccnull[,`:=`(Top10_cc=rank_cc<=10,Top20_cc=rank_cc<=20)]
fwrite(ccnull,file.path(out,'complete_covariate_null100_rankings.tsv.gz'),sep='\t')
ccmatches<-list();ccstatus<-list()
ccprimary<-ccr[method=='UPR_priority'&Top20_cc==TRUE][order(rank_cc),gene]
q<-match_one(ccprimary,'CC_UPR');ccmatches$UPR<-q$matches;ccstatus$UPR<-q$status
for(a in unique(ccnull$module)){
 selected<-ccnull[module==a&Top20_cc==TRUE][order(rank_cc),gene]
 q<-match_one(selected,paste0('CC_',a));ccmatches[[a]]<-q$matches;ccstatus[[a]]<-q$status
}
fwrite(rbindlist(ccmatches),file.path(out,'complete_covariate_background_matches.tsv'),sep='\t')
fwrite(rbindlist(ccstatus),file.path(out,'complete_covariate_matching_estimability.tsv'),sep='\t')
writeLines(c('Pre-outcome source-only availability sensitivity authorized by root before RPE1 expression.',
paste('All438 metadata-eligible genes remain the primary roster; source matching covariates are jointly complete for',length(ccgenes),'genes.'),
'The complete-covariate candidate universe is applied uniformly to realUPR, all100 null modules and every comparator; rankings are recomputed only within that same fixed sensitivity universe.',
'Primary full-roster ranks and matching failures remain unchanged and are reported; no missing expression or essentiality value is imputed.',
'This sensitivity does not replace the Holm3 primary family and cannot be selected using RPE1 outcomes.'),file.path(out,'complete_covariate_sensitivity_amendment.txt'))
fwrite(source_table[,.(genes=.N,available=sum(is.finite(depmap_essentiality)),missing=paste(gene[!is.finite(depmap_essentiality)],collapse=';'))],file.path(out,'DepMap_availability.tsv'),sep='\t')
writeLines(c('Annotation-derived comparator uses the existing gene-level Bonferroni/BH aggregation of deposited transcript energy-test P values; it is not mislabeled as an untouched original author FDR.',
'QARS is resolved through the author feature Ensembl identifier and raw DepMap Entrez5859 column, avoiding an ambiguous QARS alias.',
'No absent DepMap value or unmeasured source expression is imputed as a known negative.',
'Primary shared roster remains438 regardless of comparator availability; fair complete-comparator analyses have separately recorded coverage.'),file.path(out,'comparator_provenance_notes.txt'))
print(status[list=='UPR']);print(source_table[gene%in%primary,.(gene,UPR_priority,source_UPR_deltaZ,depmap_essentiality,source_control_expression)])
if(!status[list=='UPR',estimable]){cat('PREFREEZE BLOCK: literal primaryTop20 matching unestimable. RPE1 must remain unopened until root reviews.\n');quit(status=3L)}
cat('Source ranking/matching preflight complete. Final hash lock is a separate step; RPE1 X unread.\n')
