#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(hdf5r);library(digest)})
a<-commandArgs(TRUE);root<-if(length(a))normalizePath(a[1])else normalizePath('.')
out<-file.path(root,'paper19opp-revision2/reconstruction/candidate_validation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-as.data.table(readRDS(file.path(root,'paper19opp-revision1/GitHub/inputs/replogle_scores.rds')))
x<-x[passes_qc_final==TRUE&scDblFinder.class=='singlet'];x[,gene:=as.character(gene)]
axes<-c('apoptosis','proliferation','stress');cols<-paste0('score_',axes)
ct<-x[gene=='non-targeting'];mu<-vapply(cols,function(v)mean(ct[[v]]),numeric(1));sdv<-vapply(cols,function(v)sd(ct[[v]]),numeric(1))
z<-x[gene!='non-targeting',as.list((vapply(.SD,mean,numeric(1))-mu)/sdv),by=gene,.SDcols=cols]
setnames(z,cols,paste0('z_',axes));zmat<-as.matrix(z[,paste0('z_',axes),with=FALSE]);p<-max.col(abs(zmat),ties.method='first')
z[,target_program:=axes[p]];z[,target_sign:=sign(zmat[cbind(seq_len(nrow(z)),p)])]
z[,primary_magnitude:=abs(zmat[cbind(seq_len(nrow(z)),p)])]
z[,off_target_magnitude:=sqrt(rowSums(zmat^2)-primary_magnitude^2)]
z[,total_response_magnitude:=sqrt(rowSums(zmat^2))]
q<-fread(file.path(root,'results/task_registry/replogle_task_qc.tsv'))
z<-merge(z,q[,.(gene,discovery_cells=target_cell_count,discovery_gem_groups=replicate_groups,median_replicate_cosine,original_effect_eligible=passes_task_qc)],by='gene')
z[,priority_score:=primary_magnitude*pmax(median_replicate_cosine,0)/(1+off_target_magnitude)]
z[,metadata_eligible:=discovery_cells>=100&discovery_gem_groups>=3&is.finite(priority_score)]
h<-H5File$new('/Volumes/wang/paper19/raw/K562_gwps_raw_bulk_01.h5ad','r')
ann<-function(f){nd<-h[['obs']][[f]];if(inherits(nd,'H5Group')){co<-nd[['codes']][];cc<-nd[['categories']][];v<-rep(NA_character_,length(co));v[co>=0]<-cc[co[co>=0]+1];v}else nd[]}
expr<-data.table(gene=sub('^[^_]+_([^_]+)_.*$','\\1',ann('gene_transcript')),discovery_control_expression=ann('control_expr'));expr<-expr[,.(discovery_control_expression=mean(discovery_control_expression)),by=gene];h$close_all()
z<-merge(z,expr,by='gene',all.x=TRUE);setorderv(z,c('priority_score','gene'),c(-1,1));z[,global_discovery_rank:=seq_len(.N)]
fwrite(z,file.path(out,'locked_discovery_ranking.tsv'),sep='\t')
meta<-fread(file.path(root,'paper19opp-revision2/reconstruction/external_resources/frangieh_metadata_counts.tsv'))
meta<-meta[perturbation_2=='Control',.(gene=perturbation,validation_cells=cells,guide_labels,guide_ids)]
v<-merge(z,meta,by='gene');v<-v[metadata_eligible==TRUE&validation_cells>=30]
setorderv(v,c('priority_score','gene'),c(-1,1));v[,candidate_rank:=seq_len(.N)];v[,`:=`(top10=candidate_rank<=10,top20=candidate_rank<=20)]
fwrite(v,file.path(out,'locked_frangieh_candidates.tsv'),sep='\t')
cat('Discovery metadata eligible',sum(z$metadata_eligible),'Frangieh primary overlap',nrow(v),'\n');print(v[,.(gene,candidate_rank,priority_score,target_program,target_sign,validation_cells,discovery_cells)])
writeLines(c('Discovery: Replogle2022 K562 CRISPRi only; no MultiPerturb, Norman, Frangieh or validation effect quantities enter ranking.','Primary score = absolute dominant discovery program z-shift * max(0, median discovery GEM response cosine) / (1 + Euclidean movement on the other two programs).','Program and sign are selected entirely in Replogle; three original frozen program gene sets are unchanged.','Discovery metadata eligibility: at least100 retained cells, at least3 GEM groups, finite score. No discovery effect significance gate.','Validation: FrangiehIzar2021 SCP1064 melanoma CRISPRknockout, independent study and cellular context.','Primary population: gene present in discovery and validation metadata, discovery metadata-eligible, at least30 cells in unstimulated Control context before reading expression outcomes. All such genes are evaluated.','Primary outcome: validation signed z-shift on each gene\x27s discovery-selected program, scaled by matched non-targeting-control mean and SD in Control context. No validation-dependent program selection.','Secondary outcome: signed target z-shift divided by (1 + off-target program norm). IFNg and Co-culture contexts reported separately with the same locked rank and target definitions.','Primary ranking assessments: all-gene Spearman correlation; top10 and top20 mean directional outcome, direction concordance and gap versus unselected genes.','Simple baselines: discovery dominant-program magnitude, total3program response magnitude, GEMagreement alone and discovery baseline target-expression rank. Ties broken lexicographically by gene.','All target sets fixed before validation expression read. Equal-cell30 resampling and gene bootstrap intervals, matched non-targeting controls, matched-background/random-rank tests on discovery baseline expression and cell abundance assess uncertainty/ascertainment; no claim of biological replication from cell bootstrap.','Sensitivity: validation metadata count gates20,50,100 with full denominators; no gate selected using response effects.','This is exploratory validation proposed during revision, not preregistration before data generation. It can validate transcriptomic prioritization, not functional benefit, stable reachability or control-node status.'),file.path(out,'locked_protocol.txt'))
files<-file.path(out,c('locked_protocol.txt','locked_discovery_ranking.tsv','locked_frangieh_candidates.tsv'));fwrite(data.table(file=basename(files),sha256=vapply(files,digest,character(1),file=TRUE,algo='sha256')),file.path(out,'discovery_lock_hashes.tsv'),sep='\t')
