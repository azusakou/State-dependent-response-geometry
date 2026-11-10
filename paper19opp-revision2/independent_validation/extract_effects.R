suppressPackageStartupMessages({library(data.table);library(hdf5r);library(Matrix);library(digest)})
w<-'independent_validation';stopifnot(file.exists(file.path(w,'DOWNLOAD_COMPLETE')))
source('upr_robustness/modules/scoring_helpers.R');set.seed(2026110701)
md<-readRDS(file.path(w,'cell_metadata.rds'));u0<-fread('reconstruction/cross_dataset/phase2_upr_rpe1/phase2_gene_universe.tsv')[[1]];sy<-md$symbol;u<-u0[u0%in%sy];sp<-readRDS('reconstruction/cross_dataset/phase2_upr_rpe1/phase2_axis_spec.rds');A<-sp$UPR$positive;stopifnot(length(A)==90,all(A%in%sy))
st<-readRDS('reconstruction/upr_rpe1_validation/phase2_statistical_state.rds');targets<-st$source$gene;sel<-which(md$gene=='non-targeting'|md$gene%in%targets);ci<-which(md$gene[sel]=='non-targeting');stopifnot(length(ci)>=60)
# Fixed membership is preserved. Assayed normalization coverage is determined only
# from feature metadata; missing features are not imputed or replaced.
ui<-which(sy%in%u);M<-sparseMatrix(i=match(sy[ui],u),j=seq_along(ui),x=1,dims=c(length(u),length(ui)));ai<-match(A,u);ti<-match(targets,u)
X<-matrix(NA_real_,length(u),length(sel),dimnames=list(u,NULL));f<-H5File$new(file.path(w,'HepG2_complete.h5ad'),'r')
for(b in split(seq_along(sel),ceiling(seq_along(sel)/512))){v<-f[['X']][,sel[b],drop=FALSE];stopifnot(all(is.finite(v)),all(v>=0),max(abs(v-round(v)))<1e-6);X[,b]<-as.matrix(M%*%v[ui,,drop=FALSE]);cat('cells extracted',max(b),'/',length(sel),'\n');flush.console()};f$close_all()
L<-colSums(X);stopifnot(all(L>0));gene<-md$gene[sel];gem<-md$gem[sel];control<-X[,ci,drop=FALSE];Lc<-L[ci];gc<-gem[ci];rows<-list();cells<-list()
for(j in seq_along(targets)){
 g<-targets[j];it<-which(gene==g);if(!length(it))next;t<-ti[j];den<-L[it];denc<-Lc
 if(!is.na(t)){den<-den-X[t,it];denc<-denc-control[t,]};stopifnot(all(den>0),all(denc>0));keep<-which(A!=g)
 for(mode in c('exact_Lminus_target','fixed_L')){
  lt<-if(mode=='fixed_L')L[it] else den;lc<-if(mode=='fixed_L')Lc else denc
  q<-colMeans(log1p(sweep(X[ai[keep],it,drop=FALSE],2,lt,'/')*1e4));qc<-colMeans(log1p(sweep(control[ai[keep],,drop=FALSE],2,lc,'/')*1e4));sig<-sd(qc);mu<-tapply(qc,gc,mean);stopifnot(all(gem[it]%in%names(mu)),is.finite(sig),sig>0)
  dzcell<-(q-as.numeric(mu[gem[it]]))/sig;rows[[length(rows)+1]]<-data.table(gene=g,mode=mode,cells=length(it),GEMs=uniqueN(gem[it]),effect=mean(dzcell),control_sd=sig)
  cells[[length(cells)+1]]<-data.table(gene=g,mode=mode,gem=gem[it],effect=dzcell)
 }
}
fwrite(rbindlist(rows),file.path(w,'HepG2_gene_effects.tsv'),sep='\t');fwrite(rbindlist(cells),file.path(w,'HepG2_cell_effects.tsv.gz'),sep='\t');fwrite(data.table(gene=u),file.path(w,'assayed_frozen_universe.tsv'),sep='\t');saveRDS(list(L=L,gene=gene,gem=gem,source_rows=sel,control_indices=ci),file.path(w,'cell_accounting.rds'));writeLines(capture.output(sessionInfo()),file.path(w,'extraction_sessionInfo.txt'));cat('extraction complete\n')
