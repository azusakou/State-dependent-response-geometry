#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(BiocNeighbors)})
args<-commandArgs(trailingOnly=TRUE)
pkg<-if(length(args))args[1]else'paper19opp-revision1/GitHub'
out<-if(length(args)>1)args[2]else'paper19opp-revision2/reproducibility/lineage_simulation'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache<-readRDS(file.path(pkg,'inputs','lineage_inputs.rds'))
base_cov<-c('cell_cycle_state','stress_state','chromatin_regulator_proxy','metabolic_state','mitochondrial_fraction')
stress_axes<-c('hallmark__HYPOXIA','hallmark__REACTIVE_OXYGEN_SPECIES_PATHWAY','hallmark__UNFOLDED_PROTEIN_RESPONSE')
for(e in parse(file.path(pkg,'code','10_lineage.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c('fit_one','neighbours','geometry'))eval(e)
for(nm in names(cache$data)){
 d<-cache$data[[nm]];feat<-setdiff(intersect(cache$features,names(d)),stress_axes);branch<-grep('^p_branch_',names(d),value=TRUE)
 nn<-neighbours(d,feat)
 rhs<-c(base_cov,'log_baseline_n','log_descendant_n','log_detected_genes')
 if(uniqueN(d$future_condition)>1)rhs<-c(rhs,'factor(future_condition)')
 if(uniqueN(d$dominant_library)>1)rhs<-c(rhs,'factor(dominant_library)')
 x<-model.matrix(as.formula(paste('~',paste(rhs,collapse='+'))),d)
 stresscol<-match('stress_state',colnames(x));x0<-x[,-stresscol,drop=FALSE]
 sres<-lm.fit(x0,d$stress_state)$residuals
 fwrite(data.table(clone_id=d$clone_id,condition=d$future_condition,descendants=d$descendant_cells,stress=d$stress_state,stress_residual=sres,observed_rank=d$local_effective_rank_stress_state_loo),file.path(out,paste0(nm,'_design.tsv')),sep='\t')
 fwrite(data.table(row=rep(seq_len(nrow(d)),lengths(nn$index)),neighbor=unlist(nn$index)),file.path(out,paste0(nm,'_neighbors.tsv')),sep='\t')
 fwrite(data.table(row=seq_len(nrow(d)),as.data.table(d[,..branch])),file.path(out,paste0(nm,'_branches.tsv')),sep='\t')
 # Rebuilt geometry for selected deeper donor/query population is a separate design.
 ds<-d[descendant_cells>=10];nns<-neighbours(ds,feat)
 rhs2<-c(base_cov,'log_baseline_n','log_descendant_n','log_detected_genes')
 if(uniqueN(ds$future_condition)>1)rhs2<-c(rhs2,'factor(future_condition)')
 if(uniqueN(ds$dominant_library)>1)rhs2<-c(rhs2,'factor(dominant_library)')
 xx<-model.matrix(as.formula(paste('~',paste(rhs2,collapse='+'))),ds);jj<-match('stress_state',colnames(xx));sres2<-lm.fit(xx[,-jj,drop=FALSE],ds$stress_state)$residuals
 fwrite(data.table(original_row=which(d$descendant_cells>=10),stress_residual=sres2),file.path(out,paste0(nm,'_deep_design.tsv')),sep='\t')
 fwrite(data.table(row=rep(seq_len(nrow(ds)),lengths(nns$index)),neighbor=unlist(nns$index)),file.path(out,paste0(nm,'_deep_neighbors.tsv')),sep='\t')
 # Verify the optimized coefficient agrees with the original regression.
 rank<-geometry(as.matrix(d[,..branch]),nn)$rank
 stopifnot(max(abs(rank-d$local_effective_rank_stress_state_loo))<1e-8)
 slope<-sum(sres*rank)/sum(sres^2)
 stopifnot(abs(slope-fit_one(d,rank,TRUE)$beta)<1e-9)
 cat(nm,'observed adjusted slope',slope,'\n')
}
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
