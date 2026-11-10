suppressPackageStartupMessages(library(data.table));set.seed(2026110703);w<-'independent_validation';d<-fread(file.path(w,'HepG2_cell_effects.tsv.gz'));v<-fread(file.path(w,'HepG2_gene_effects.tsv'));s<-fread('reconstruction/cross_dataset/phase2_upr_rpe1/source_rank_inputs.tsv');rk<-fread('reconstruction/cross_dataset/phase2_upr_rpe1/all_source_comparator_rankings.tsv')[method=='UPR_priority'][order(rank),gene];ans<-list();draws<-list()
for(mode_ in unique(v$mode)){
 z<-d[mode==mode_];gg<-sort(unique(z$gem));gs<-sort(unique(z$gene));a<-z[,.(mean=mean(effect),n=.N),by=.(gene,gem)];base<-v[mode==mode_];a<-merge(a,base[,.(gene,overall=effect,total=cells)],by='gene');R<-matrix(0,length(gs),length(gg),dimnames=list(gs,gg));R[cbind(match(a$gene,gs),match(a$gem,gg))]<-a$n/a$total*(a$mean-a$overall);W<-matrix(sample(c(-1,1),length(gg)*2000,TRUE),length(gg),2000);Y<-R%*%W+base$effect[match(gs,base$gene)];rownames(Y)<-gs
 for(dep in c(30L,10L,100L)){
  eligible<-base[cells>=dep,gene];eligible<-intersect(eligible,s$gene);xx<-s$source_UPR_deltaZ[match(eligible,s$gene)];r<-apply(Y[eligible,,drop=FALSE],2,function(y)cor(xx,y,method='spearman'));dd<-data.table(mode=mode_,minimum_cells=dep,draw=seq_len(2000),rho=r)
  for(K in c(10L,20L)){ids<-head(rk,K);other<-setdiff(eligible,ids);gap<-if(all(ids%in%eligible))colMeans(Y[ids,,drop=FALSE])-colMeans(Y[other,,drop=FALSE]) else rep(NA_real_,2000);dd[[paste0('gap',K)]]<-gap}
  draws[[length(draws)+1]]<-dd
  for(nm in c('rho','gap10','gap20')){q<-dd[[nm]];q<-q[is.finite(q)];ans[[length(ans)+1]]<-data.table(mode=mode_,minimum_cells=dep,metric=nm,finite_draws=length(q),low=if(length(q))quantile(q,.025) else NA_real_,high=if(length(q))quantile(q,.975) else NA_real_)}
 }
}
fwrite(rbindlist(ans),file.path(w,'shared_GEM_intervals.tsv'),sep='\t');fwrite(rbindlist(draws),file.path(w,'shared_GEM_draws.tsv.gz'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(w,'GEM_sessionInfo.txt'))
