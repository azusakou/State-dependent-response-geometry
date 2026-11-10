#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);project<-normalizePath(args[1]);out<-normalizePath(args[2])
.libPaths(c(file.path(project,'.Rlib'),.libPaths()));suppressPackageStartupMessages({library(data.table);library(digest)})
stopifnot(!file.exists(file.path(out,'RPE1_expression_read_started.txt')))
proposal<-readLines(file.path(dirname(out),'phase2_upr_rpe1_protocol_proposal.txt'))
proposal[grepl('^STATUS:',proposal)]<-'STATUS: APPROVED FOLLOW-UP ANALYSIS, LOCKED BEFORE RPE1 EXPRESSION ACCESS. This is outcome-blind retrospective validation, not a prospective experiment or a public preregistration.'
writeLines(c(proposal,'',
'Pre-outcome authorized comparator amendment: annotation-derived genome-wide response significance is added to the exploratory comparator family; its gene-level FDR aggregates deposited energy-test P values and is not described as an untouched author FDR.',
'Pre-outcome authorized availability sensitivity: all438 primary genes remain. The370-gene complete four-covariate roster is applied uniformly to realUPR, all100 null modules and all comparators; original84 unmatched null flags are retained.',
'Pre-outcome conditional-GEM refinement: shared Rademacher GEM residual perturbations preserve sparse genes and shared control references instead of deleting genes that lose all sampled GEMs. Ordinary whole-gene and matched-set intervals are also retained. These are conditional technical checks, not independent biological replication.'),file.path(out,'locked_analysis_protocol.txt'))
status<-fread(file.path(out,'matching_estimability_before_RPE1.tsv'));stopifnot(status[list=='UPR',estimable])
cc<-fread(file.path(out,'complete_covariate_matching_estimability.tsv'));stopifnot(nrow(cc)==101L,all(cc$estimable))
genes<-fread(file.path(out,'phase2_eligible_genes.tsv'));stopifnot(nrow(genes)==438L)
nulls<-fread(file.path(out,'null100_module_membership.tsv'));stopifnot(nrow(nulls)==9000L,uniqueN(nulls$module)==100L)
lists<-fread(file.path(out,'source_Top10_Top20_lists.tsv'))[method=='UPR_priority'];stopifnot(nrow(lists)==20L,sum(lists$Top10)==10L)
matches<-fread(file.path(out,'frozen_background_matches.tsv'))[list=='UPR'];stopifnot(nrow(matches)==100L,uniqueN(matches$selected_gene)==20L,uniqueN(matches$background_gene)==100L)
stat_files<-c('../../upr_rpe1_validation/phase2_statistics.R','../../upr_rpe1_validation/run_phase2_statistics.R','../../upr_rpe1_validation/phase2_statistics_preflight.R','../../upr_rpe1_validation/STATISTICS_INTERFACE.md')
stopifnot(all(file.exists(file.path(out,stat_files))))
local_files<-list.files(out,full.names=FALSE)
local_files<-local_files[grepl('\\.(R|rds|tsv|tsv.gz|txt)$',local_files)&!grepl('^(source_lock_manifest|source_lock_record|RPE1_)',local_files)]
# The metadata-only RPE1 roster and the unopened reader are themselves locked.
local_files<-unique(c(local_files,'RPE1_author_roster_metadata_only.tsv.gz','RPE1_extract_once.R'))
paths<-unique(c(local_files,stat_files));stopifnot(all(file.exists(file.path(out,paths))))
hashes<-vapply(file.path(out,paths),digest,character(1),file=TRUE,algo='sha256')
manifest<-data.table(file=paths,sha256=hashes);fwrite(manifest,file.path(out,'source_lock_manifest.tsv'),sep='\t')
manifest_hash<-digest(file.path(out,'source_lock_manifest.tsv'),file=TRUE,algo='sha256')
writeLines(c(paste('UTC',format(Sys.time(),tz='UTC',format='%Y-%m-%dT%H:%M:%SZ')),
paste('manifest_sha256',manifest_hash),paste('locked_files',nrow(manifest)),
'RPE1 expression matrix and expression summaries have not been accessed before this lock.',
'Source roster438, UPR90, denominator6780, null100, all source ranks, literalTop10/20, all matches and370-gene uniform sensitivity were locked together.',
'No missing selected-UPR matching covariates; all20 original members retained. No RPE1-based scoring, matching, cutoff or model choice is permitted.'),file.path(out,'source_lock_record.txt'))
cat('SOURCE_LOCK_COMPLETE',manifest_hash,'files',nrow(manifest),'\n')
