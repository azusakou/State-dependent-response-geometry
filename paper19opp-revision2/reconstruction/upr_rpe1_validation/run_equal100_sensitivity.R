#!/usr/bin/env Rscript
# Declared equal-target-cell support sensitivity, with frozen ranks and control references.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
a<-commandArgs(TRUE);out<-if(length(a))a[1]else'paper19opp-revision2/reconstruction/upr_rpe1_validation'
here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))))
source(file.path(here,'phase2_statistics.R'))
z<-readRDS(file.path(out,'phase2_statistical_state.rds'));src<-z$source_directory;genes<-z$source$gene;B<-200L;depth<-100L;set.seed(20261111L)
spec<-readRDS(file.path(src,'phase2_axis_spec.rds'));members<-spec$UPR$positive
if(is.null(members))members<-spec$UPR$pos
if(is.null(members))stop('UPR member schema not recognized')
stopifnot(length(members)==90L)
r<-readRDS(file.path(src,'original91_source_raw_cache.rds'));m<-as.data.table(r$cells);u<-r$raw_upr[,members,drop=FALSE];stopifnot(all(m$phase2_library>0))
u@x<-log1p(u@x*10000/m$phase2_library[u@i+1L]);us<-Matrix::rowSums(u);own<-m$gene%in%members
s<-as.numeric((us-ifelse(own,log1p(10000*m$target_raw_count/m$phase2_library),0))/(90L-as.integer(own)))
sg<-fread(file.path(src,'source_phase2_GEM_effects.tsv.gz'),select=c('gene','axis','gem_group','treated_mean','deltaZ','control_sd'))[axis=='UPR'&gene%in%genes]
meancheck<-data.table(gene=m$gene,gem_group=m$gem_group,score=s)[gene%in%genes,.(new_mean=mean(score)),by=.(gene,gem_group)]
chk<-merge(meancheck,sg,by=c('gene','gem_group'));stopifnot(nrow(chk)==nrow(sg),max(abs(chk$new_mean-chk$treated_mean))<1e-10)
make_values<-function(meta,score,gem,point,sd){
 d<-data.table(gene=meta$gene,gem=as.character(meta$gem_group),score=score)[gene%in%genes]
 d[,unit_mean:=mean(score),by=.(gene,gem)];gg<-as.data.table(gem)[,.(gene,gem=as.character(gem),contrast_z)];d<-merge(d,gg,by=c('gene','gem'),all.x=TRUE,sort=FALSE)
 stopifnot(!anyNA(d$contrast_z));d[,value:=(score-unit_mean)/sd[gene]+contrast_z]
 vals<-split(d$value,d$gene);stopifnot(all(lengths(vals)>=depth));check<-vapply(vals,mean,numeric(1));stopifnot(max(abs(check-point[names(check)]))<1e-9);vals
}
sdsrc<-setNames(z$source_metadata$source_UPR_control_sd,z$source_metadata$gene)
sv<-make_values(m,s,z$source_gem,z$primary$point_source,sdsrc);rm(r,u,m,s);gc(FALSE)
v<-readRDS(z$validation_cache);vm<-as.data.table(v$cell_metadata);sdval<-setNames(z$validation_gene_effects[axis=='UPR',control_sd],z$validation_gene_effects[axis=='UPR',gene]);vv<-make_values(vm,v$cell_scores[,'UPR'],z$validation_gem,z$primary$point_validation,sdval);rm(v,vm);gc(FALSE)
sample_means<-function(vals){ans<-t(vapply(genes,function(g){v<-vals[[g]];vapply(seq_len(B),function(b)mean(v[sample.int(length(v),depth,replace=FALSE)]),numeric(1))},numeric(B)));rownames(ans)<-genes;ans}
xs<-sample_means(sv);ys<-sample_means(vv)
dd<-rbindlist(lapply(seq_len(B),function(b){y<-setNames(ys[,b],genes);gap<-phase2_set_contrasts(y,z$primary$sets);data.table(draw=b,rho=phase2_spearman(xs[,b],ys[,b]),Top20_matched_gap=mean(gap),Top10_matched_gap=mean(gap[1:10]),Top20_mean=mean(y[names(z$primary$sets)]),Top10_mean=mean(y[names(z$primary$sets)[1:10]]))}))
fwrite(dd,file.path(out,'equal100_target_cell_draws.tsv'),sep='\t')
ss<-rbindlist(lapply(setdiff(names(dd),'draw'),function(n){x<-dd[[n]];data.table(endpoint=n,median=median(x),resampling_low=quantile(x,.025),resampling_high=quantile(x,.975),draws=B,target_cells_each_context=depth,positive_draw_fraction=mean(x>0))}))
ss[,note:='without-replacement target-cell subsampling; original GEM-matched NTC references and global target-specific SD fixed; frozen lists']
fwrite(ss,file.path(out,'equal100_target_cell_summary.tsv'),sep='\t')
ccr<-fread(file.path(src,'complete_covariate_source_rankings.tsv'));cmp<-rbindlist(lapply(c(20L,10L),function(k){raw<-ccr[method=='raw_UPR'&rank_cc<=k,gene];rbindlist(lapply(unique(ccr$method),function(mm){gg<-ccr[method==mm&rank_cc<=k,gene];v<-colMeans(ys[gg,,drop=FALSE])-colMeans(ys[raw,,drop=FALSE]);data.table(method=mm,TopK=k,mean_RPE1=mean(colMeans(ys[gg,,drop=FALSE])),minus_raw_UPR_median=median(v),resampling_low=quantile(v,.025),resampling_high=quantile(v,.975))}))}))
fwrite(cmp,file.path(out,'equal100_fixed_baseline_sensitivity.tsv'),sep='\t')
saveRDS(list(source=xs,validation=ys,genes=genes,depth=depth,draws=B),file.path(out,'equal100_target_cell_effect_vectors.rds'),compress='xz')
print(ss);print(cmp)
