suppressPackageStartupMessages({library(data.table);library(digest)})
w<-'independent_validation';set.seed(2026110702);B<-4999L;P<-10000L
v<-fread(file.path(w,'HepG2_gene_effects.tsv'));s<-fread('reconstruction/cross_dataset/phase2_upr_rpe1/source_rank_inputs.tsv');rk<-fread('reconstruction/cross_dataset/phase2_upr_rpe1/all_source_comparator_rankings.tsv');nr<-fread('reconstruction/cross_dataset/phase2_upr_rpe1/null100_source_rankings.tsv.gz');src<-setNames(s$source_UPR_deltaZ,s$gene)
methods<-unique(rk$method);CC<-Reduce(intersect,lapply(methods,function(m)rk[method==m&available==TRUE,gene]));ranks<-lapply(methods,function(m)rk[method==m&gene%in%CC][order(rank),gene]);names(ranks)<-methods
pri<-rk[method=='UPR_priority'][order(rank),gene];rawname<-methods[grepl('raw.*UPR|UPR.*raw',methods)][1];stopifnot(length(pri)==438,!is.na(rawname));literal<-list(priority=pri,raw=rk[method==rawname][order(rank),gene]);q95<-function(x){x<-x[is.finite(x)];if(!length(x))return(c(NA_real_,NA_real_));as.numeric(quantile(x,c(.025,.975)))}
metric<-function(ids,y,weights=NULL){
 K<-length(ids);if(!all(ids%in%names(y)))return(c(Precision=NA_real_,NDCG=NA_real_,gap=NA_real_));o<-setdiff(names(y),ids);if(!length(o))return(rep(NA_real_,3));idcg<-sum(sort(pmax(y,0),decreasing=TRUE)[seq_len(K)]/log2(seq_len(K)+1));if(is.null(weights))weights<-setNames(rep(1,length(y)),names(y));x<-y[ids];ww<-weights[ids];wo<-weights[o];if(sum(ww)==0||sum(wo)==0)return(c(Precision=NA_real_,NDCG=NA_real_,gap=NA_real_));c(Precision=sum(ww*(x>0))/sum(ww),NDCG=if(idcg>0)sum(ww*pmax(x,0)/log2(seq_len(K)+1))/idcg else NA_real_,gap=sum(ww*x)/sum(ww)-sum(wo*y[o])/sum(wo))
}
allm<-alln<-allr<-allb<-allw<-list();test<-list();bi<-0L
for(mode_ in c('exact_Lminus_target','fixed_L'))for(depth in c(30L,10L,100L)){
 vv<-v[mode==mode_&cells>=depth&is.finite(effect)];ids<-intersect(s$gene,vv$gene);y<-setNames(vv$effect,vv$gene)[ids];x<-src[ids];n<-length(y);rho<-cor(x,y,method='spearman');bootrho<-replicate(B,{ii<-sample.int(n,n,TRUE);cor(x[ii],y[ii],method='spearman')});ci<-q95(bootrho);p<-if(depth==30&&mode_=='exact_Lminus_target')(1+sum(abs(replicate(P,cor(x,sample(y),method='spearman')))>=abs(rho)))/(P+1) else NA_real_
 allr[[length(allr)+1]]<-data.table(mode=mode_,minimum_cells=depth,genes=n,rho=rho,low=ci[1],high=ci[2],raw_p=p,sign_concordance=mean(sign(x)==sign(y)))
 if(depth==30&&mode_=='exact_Lminus_target')test[[1]]<-data.table(endpoint='signed_UPR_rho',estimate=rho,raw_p=p,estimable=TRUE)
 yc<-y[intersect(CC,names(y))]
 for(m in methods)for(K in c(10L,20L)){
  chosen<-head(ranks[[m]],K);met<-metric(chosen,yc);allm[[length(allm)+1]]<-data.table(mode=mode_,minimum_cells=depth,cohort='fixed_common_source_covariate_panel',method=m,K=K,cohort_genes=length(yc),list_genes_available=sum(chosen%in%names(yc)),estimable=all(chosen%in%names(yc)),Precision=met[1],NDCG=met[2],gap=met[3])
 }
 # Primary literal lists retain all438-source membership, never refill.
 for(K in c(10L,20L)){
  chosen<-head(pri,K);met<-metric(chosen,y);if(depth==30&&mode_=='exact_Lminus_target')test[[length(test)+1]]<-data.table(endpoint=paste0('literal_Top',K,'_signed_gap'),estimate=met[3],raw_p=if(is.na(met[3]))1 else (1+sum(abs(replicate(P,metric(sample(names(y),K),y)[3]))>=abs(met[3])))/(P+1),estimable=is.finite(met[3]))
  rawids<-head(literal$raw,K);delta<-metric(chosen,y)[2]-metric(rawids,y)[2]
  if(depth==30&&mode_=='exact_Lminus_target')test[[length(test)+1]]<-data.table(endpoint=paste0('weighted_minus_raw_NDCG',K),estimate=delta,raw_p=if(is.na(delta))1 else (1+sum(abs(replicate(P,{yp<-setNames(sample(y),names(y));metric(chosen,yp)[2]-metric(rawids,yp)[2]}))>=abs(delta)))/(P+1),estimable=is.finite(delta))
  bb<-matrix(NA_real_,B,4);colnames(bb)<-c('Precision','NDCG','gap','weighted_minus_raw_NDCG');for(b in seq_len(B)){ww<-setNames(tabulate(sample.int(n,n,TRUE),nbins=n),names(y));a<-metric(chosen,y,ww);bb[b,]<-c(a,a[2]-metric(rawids,y,ww)[2])};allb[[length(allb)+1]]<-cbind(data.table(mode=mode_,minimum_cells=depth,K=K,metric=colnames(bb),estimate=c(met,delta)),as.data.table(t(apply(bb,2,q95))));setnames(allb[[length(allb)]],c('V1','V2'),c('low','high'))
 }
 # Random module priorities predict the SAME frozen UPR endpoint. Their validation
 # expression sets are not redefined, and unavailable literal lists remain NA.
 for(mm in unique(nr$module))for(K in c(10L,20L)){chosen<-nr[module==mm&gene%in%CC][order(rank),head(gene,K)];met<-metric(chosen,yc);alln[[length(alln)+1]]<-data.table(mode=mode_,minimum_cells=depth,module=mm,K=K,available=sum(chosen%in%names(yc)),Precision=met[1],NDCG=met[2],gap=met[3])}
 for(K in c(10L,20L)){u<-replicate(P,metric(sample(names(yc),K),yc));allw[[length(allw)+1]]<-data.table(mode=mode_,minimum_cells=depth,K=K,draw=seq_len(P),Precision=u[1,],NDCG=u[2,],gap=u[3,])}
 cat(mode_,depth,'complete; genes',n,'rho',rho,'\n');flush.console()
}
t<-rbindlist(test);stopifnot(nrow(t)==5);t[,Holm5:=p.adjust(raw_p,'holm',n=5)];fwrite(t,file.path(w,'primary_Holm5.tsv'),sep='\t');fwrite(rbindlist(allr),file.path(w,'all_ordering_results.tsv'),sep='\t');fwrite(rbindlist(allm),file.path(w,'all_fixed_baseline_metrics.tsv'),sep='\t');fwrite(rbindlist(allb),file.path(w,'literal_list_gene_intervals.tsv'),sep='\t');fwrite(rbindlist(alln),file.path(w,'all100_module_priority_metrics.tsv'),sep='\t');fwrite(rbindlist(allw),file.path(w,'uniform10000_reference_draws.tsv.gz'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(w,'statistics_sessionInfo.txt'));print(t)
