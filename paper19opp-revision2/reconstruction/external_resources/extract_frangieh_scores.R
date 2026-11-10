#!/usr/bin/env Rscript
# Independent study validation scoring. Selection/signs remain discovery-frozen.
suppressPackageStartupMessages({library(data.table);library(hdf5r);library(yaml)})
out<-'paper19opp-revision2/reconstruction/external_resources'
sourcefile<-'/Volumes/wang/paper13/raw/perturbation/SCP1064/FrangiehIzar2021_RNA.h5ad'
programfile<-'paper19opp-revision1/GitHub/upstream/program_gene_sets.yaml'
lockfile<-'paper19opp-revision2/reconstruction/candidate_validation/discovery_lock_hashes.tsv'
stopifnot(file.exists(lockfile))
f<-H5File$new(sourcefile,'r');on.exit(f$close_all())
readcol<-function(name){o<-f[['obs']][[name]];if(inherits(o,'H5Group')){catg<-o[['categories']][];codes<-o[['codes']][];return(catg[codes+1L])};o[]}
d<-data.table(cell_id=readcol('cell_name'),gene=readcol('perturbation'),context=readcol('perturbation_2'),sgRNA=readcol('sgRNA'),guide_id=readcol('guide_id'),nperts=readcol('nperts'),library_size=readcol('ncounts'),detected_genes=readcol('ngenes'),percent_mito=readcol('percent_mito'))
d[,is_control:=gene=='control']
overlap<-fread(file.path(out,'frangieh_discovery_overlap_metadata.tsv'))[in_replogle==TRUE,gene]
keep<-which(d$is_control|d$gene%in%overlap)
d<-d[keep];stopifnot(!anyDuplicated(d$cell_id),all(d$library_size>0),all(d$nperts%in%c(0L,1L)))
ng<-f[['var']][['gene_symbol']][];stopifnot(!anyDuplicated(ng))
ptr<-f[['X']][['indptr']][];indices<-f[['X']][['indices']];values<-f[['X']][['data']]
# CSC by gene: only frozen program and common target columns are accessed.
get_gene<-function(gene){j<-match(gene,ng);if(is.na(j))return(rep(NA_real_,nrow(d)));lo<-ptr[j]+1;hi<-ptr[j+1];v<-numeric(length(readcol('cell_name')));if(hi>=lo){ii<-indices[lo:hi]+1;vv<-values[lo:hi];stopifnot(all(is.finite(vv)),all(vv>=0),max(abs(vv-round(vv)))<1e-5);v[ii]<-vv};v[keep]}
common<-fread('paper19opp-revision2/reconstruction/cross_dataset/final_gene_universe.tsv')$gene
stopifnot(length(common)==7848,all(common%in%ng))
axis_table<-fread('paper19opp-revision2/reconstruction/cross_dataset/final_axis_membership.tsv')
programs<-lapply(split(axis_table[in_common==TRUE],by='program'),function(x)list(positive_genes=x[direction=='positive_genes',gene],negative_genes=x[direction=='negative_genes',gene]))
stopifnot(length(programs)==12L)
allgenes<-unique(c(unlist(lapply(programs,function(x)c(x$positive_genes,x$negative_genes))),overlap))
coverage<-list();cache<-list();nall<-length(readcol('cell_name'))
# Re-declare without repeated metadata reads for each gene.
get_gene<-function(gene){j<-match(gene,ng);if(is.na(j))return(rep(NA_real_,nrow(d)));lo<-ptr[j]+1;hi<-ptr[j+1];v<-numeric(nall);if(hi>=lo){ii<-indices[lo:hi]+1;vv<-values[lo:hi];stopifnot(all(is.finite(vv)),all(vv>=0),max(abs(vv-round(vv)))<1e-5);v[ii]<-vv};v[keep]}
for(g in allgenes){raw<-get_gene(g);cache[[g]]<-raw}
# Exactly the metadata-declared shared measured universe; no missing genes imputed.
rowmap<-integer(nall);rowmap[keep]<-seq_along(keep);common_library<-numeric(nrow(d))
for(g in common){j<-match(g,ng);lo<-ptr[j]+1;hi<-ptr[j+1];if(hi>=lo){ii<-rowmap[indices[lo:hi]+1L];vv<-values[lo:hi];ok<-ii>0L;common_library[ii[ok]]<-common_library[ii[ok]]+vv[ok]}}
stopifnot(all(common_library>0),all(common_library<=d$library_size+1e-5))
d[,library_size_common:=common_library]
for(a in names(programs)){
 p<-intersect(programs[[a]]$positive_genes,ng);n<-intersect(programs[[a]]$negative_genes,ng)
 for(denom_name in c('library_size_common','library_size')){denom<-d[[denom_name]]
 score<-Reduce('+',lapply(p,function(g)log1p(1e4*cache[[g]]/denom)))/length(p)
 if(length(n))score<-score-Reduce('+',lapply(n,function(g)log1p(1e4*cache[[g]]/denom)))/length(n)
 field<-if(denom_name=='library_size_common')paste0('score_',a)else paste0('score_full_library_',a)
 d[,(field):=score]
 }
 coverage[[a]]<-data.table(program=a,positive_expected=length(programs[[a]]$positive_genes),positive_present=length(p),negative_expected=length(programs[[a]]$negative_genes),negative_present=length(n),positive_symbols=paste(p,collapse=','),negative_symbols=paste(n,collapse=','))
}
d[,`:=`(target_count=NA_real_,target_logcp10k=NA_real_,target_logcp10k_full_library=NA_real_,simple_single_guide=!is.na(sgRNA)&sgRNA!='nan'&sgRNA==guide_id)]
for(g in overlap){idx<-which(d$gene==g);if(!length(idx))next;raw<-cache[[g]][idx];d[idx,`:=`(target_count=raw,target_logcp10k=log1p(1e4*raw/library_size_common),target_logcp10k_full_library=log1p(1e4*raw/library_size))]}
controls<-rbindlist(lapply(overlap,function(g){raw<-cache[[g]];if(all(is.na(raw)))return(NULL);rbindlist(lapply(unique(d$context),function(co){idx<-which(d$is_control&d$context==co);z<-log1p(1e4*raw[idx]/d$library_size_common[idx]);zf<-log1p(1e4*raw[idx]/d$library_size[idx]);data.table(gene=g,context=co,control_cells=length(idx),control_mean_count=mean(raw[idx]),control_mean_logcp10k=mean(z),control_sd_logcp10k=sd(z),control_mean_logcp10k_full_library=mean(zf),control_sd_logcp10k_full_library=sd(zf),detected_fraction=mean(raw[idx]>0))}))}))
fwrite(d,file.path(out,'frangieh_program_scored_cells.tsv.gz'),sep='\t')
fwrite(controls,file.path(out,'frangieh_control_target_expression.tsv'),sep='\t')
fwrite(rbindlist(coverage),file.path(out,'frangieh_program_coverage.tsv'),sep='\t')
fwrite(d[,.(cells=.N,known_simple_single_guide_cells=sum(simple_single_guide)),by=.(gene,context)],file.path(out,'frangieh_scoring_accounting.tsv'),sep='\t')
writeLines(c('Source: independent patient-derived melanoma CRISPR-Cas9 screen, Frangieh et al. Nat Genet 2021, SCP1064; 218,331 published QC cells.',
'Selected: all 51 perturbation genes overlapping source Replogle K562, plus every non-targeting/intergenic cell annotated control in harmonized object. Selection used labels/counts only.',
'Inputs: unmodified integer count X, 23,712 unique gene symbols, CSC by gene. The metadata-locked 7,848 genes measured across Replogle, Norman, MIX-Seq, Adamson and Frangieh provide the common library denominator; frozen modules are restricted to these shared measured members.',
'Primary per-cell normalization: log1p(1e4*count/library_size_common), positive-gene mean minus negative-gene mean. Original full-library normalization with the same matched module members is retained in score_full_library_* columns. Twelve axes include three matched original programs and nine metadata-frozen Hallmark programs. No validation-derived program, sign or gene selection.',
'Controls are context matched. target_count/target_logcp10k is missing for control cells; per-target context-specific control reference expression is in a separate table.',
'Guide labels nan or multiguide construct IDs remain explicit; simple_single_guide records known sgRNA==guide_id only. These are not biological replicate labels.',
'No pretreatment expression is paired to perturbed cells. Control, IFNgamma and co-culture are experimental treatment contexts, not per-cell measured baseline states.',
'External outcome interpretation is held until the candidate agent refreezes discovery ranking using the same 7,848-gene normalization; the original ranking is retained separately as a sensitivity result. Primary inferential units for external ranking are held-out genes; contexts share genes and source experiment. Surviving-cell transcriptomic scores do not directly measure killing, division or fate probabilities.',
 paste('Executed cells',nrow(d),'target genes',length(overlap))),file.path(out,'frangieh_scoring_provenance.txt'))
f$close_all();cat('Wrote',nrow(d),'cells;51common genes and controls; no response effects summarized.\n')
