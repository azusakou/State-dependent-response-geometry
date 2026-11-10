#!/usr/bin/env Rscript
args<-commandArgs(TRUE);project<-normalizePath(if(length(args))args[1]else'paper19opp-revision2');p<-file.path(project,'upr_robustness','ranking');historical<-file.path(p,'alternative_definitions');replay<-file.path(p,'source_reproduction');dir.create(replay,showWarnings=FALSE)
# Execute the historically locked source-only calculation into a clean, separate output
# directory. The historical manifest and first-read provenance are never overwritten.
expressions<-parse(file.path(p,'freeze_alternative_source_rankings.R'));env<-new.env(parent=globalenv())
for(e in expressions){if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('out')))e<-substitute(out<-PATH,list(PATH=replay));eval(e,envir=env)}
files<-c('source_input_hashes.tsv','source_GEM_consistency_QA.tsv','source_alternative_rankings.tsv.gz','source_alternative_Top10_Top20.tsv','source_parent_specific_null_rankings.tsv.gz','source_rank_stability_vs_current.tsv','source_alternative_background_matches.tsv','source_alternative_matching_coverage.tsv')
checks<-lapply(files,function(n){a<-data.table::fread(file.path(historical,n));b<-data.table::fread(file.path(replay,n));data.frame(file=n,passed=isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))})
checks<-do.call(rbind,checks);write.table(checks,file.path(replay,'source_reproduction_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE);stopifnot(all(checks$passed));cat('Historical source ranks/matches reproduced in separate directory; all8 table checks passed.\n')
