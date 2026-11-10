suppressPackageStartupMessages(library(data.table))

sp_finite_range<-function(v){v<-v[is.finite(v)];if(length(v))range(v)else c(NA_real_,NA_real_)}
sp_ci<-function(v,min_valid=1000L){v<-v[is.finite(v)];if(length(v)<min_valid)c(NA_real_,NA_real_)else as.numeric(quantile(v,c(.025,.975)))}
sp_rho<-function(x,y){if(length(x)<3L||any(!is.finite(x))||any(!is.finite(y))||sd(x)==0||sd(y)==0)NA_real_ else suppressWarnings(cor(x,y,method='spearman'))}
sp_p<-function(observed,null){if(!is.finite(observed)||!all(is.finite(null)))return(NA_real_);mu<-mean(null);(1+sum(abs(null-mu)>=abs(observed-mu)-1e-12))/(length(null)+1)}
sp_adjust<-function(p,n=length(p),method='BH'){a<-p.adjust(ifelse(is.finite(p),p,1),method,n=n);a[!is.finite(p)]<-NA_real_;a}
sp_agreement<-function(x,y){
 ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];keep<-x!=0&y!=0;xx<-x[keep]>0;yy<-y[keep]>0;n<-length(xx)
 if(!n)return(c(concordance=NA_real_,kappa=NA_real_,MCC=NA_real_,nonzero_pairs=0))
 a<-as.numeric(sum(xx&yy));b<-as.numeric(sum(xx&!yy));c<-as.numeric(sum(!xx&yy));d<-as.numeric(sum(!xx&!yy));acc<-(a+d)/n;pe<-((a+b)*(a+c)+(c+d)*(b+d))/n^2;den<-sqrt((a+b)*(a+c)*(d+b)*(d+c))
 c(concordance=acc,kappa=if(pe<1)(acc-pe)/(1-pe)else NA_real_,MCC=if(den>0)(a*d-b*c)/den else NA_real_,nonzero_pairs=n)
}
sp_perm_index<-function(n,P){if(n<2L)return(matrix(seq_len(n),n,P));replicate(P,sample.int(n,n,replace=FALSE))}
sp_perm_metrics<-function(x,y,index){
 n<-length(x);P<-ncol(index);out<-matrix(NA_real_,P,3,dimnames=list(NULL,c('rho','direction','C')))
 if(n<1L)return(out)
 if(n>=3L&&sd(x)>0&&sd(y)>0){rx<-rank(x)-mean(rank(x));ry<-rank(y)-mean(rank(y));out[,'rho']<-as.numeric(crossprod(rx,matrix(ry[index],n,P)))/sqrt(sum(rx^2)*sum(ry^2))}
 out[,'C']<-as.numeric(crossprod(sign(x),matrix(y[index],n,P)))/n
 nz<-x!=0&y!=0
 if(sum(nz)>0){sx<-sign(x[nz]);sy<-sign(y[nz]);ip<-if(all(nz))index else sp_perm_index(sum(nz),P);out[,'direction']<-colMeans(matrix(sy[ip],sum(nz),P)==sx)}
 out
}
sp_pair<-function(d,B=4999L,P=19999L,permutation_index=NULL,bootstrap_index=NULL){
 d<-as.data.table(d);d<-d[is.finite(source_deltaZ)&is.finite(validation_deltaZ)];setorder(d,gene);x<-d$source_deltaZ;y<-d$validation_deltaZ;n<-nrow(d)
 point_ag<-sp_agreement(x,y);point<-c(rho=sp_rho(x,y),direction=unname(point_ag['concordance']),C=if(n)mean(sign(x)*y)else NA_real_,kappa=unname(point_ag['kappa']),MCC=unname(point_ag['MCC']))
 if(is.null(permutation_index))permutation_index<-sp_perm_index(n,P)
 if(is.null(bootstrap_index))bootstrap_index<-if(n)replicate(B,sample.int(n,n,replace=TRUE))else matrix(integer(),0,B)
 boot<-t(vapply(seq_len(B),function(b){ii<-bootstrap_index[,b];xx<-x[ii];yy<-y[ii];a<-sp_agreement(xx,yy);c(rho=sp_rho(xx,yy),direction=unname(a['concordance']),C=if(n)mean(sign(xx)*yy)else NA_real_,kappa=unname(a['kappa']),MCC=unname(a['MCC']))},numeric(5)))
 null<-sp_perm_metrics(x,y,permutation_index)
 rows<-rbindlist(lapply(names(point),function(k){q<-sp_ci(boot[,k]);data.table(metric=k,estimate=point[k],gene_bootstrap_low=q[1],gene_bootstrap_high=q[2],finite_bootstrap_draws=sum(is.finite(boot[,k])),raw_permutation_p=if(k%in%colnames(null))sp_p(point[k],null[,k])else NA_real_)}))
 rows[,`:=`(genes=n,source_positive=sum(x>0),source_negative=sum(x<0),source_zero=sum(x==0),validation_positive=sum(y>0),validation_negative=sum(y<0),validation_zero=sum(y==0),nonzero_pairs=point_ag['nonzero_pairs'])]
 loo<-if(n>=4L)rbindlist(lapply(seq_len(n),function(i)data.table(omitted_gene=d$gene[i],rho=sp_rho(x[-i],y[-i]),C=mean(sign(x[-i])*y[-i]))))else data.table(omitted_gene=character(),rho=numeric(),C=numeric())
 list(results=rows,bootstrap=as.data.table(cbind(draw=seq_len(B),boot)),permutation=as.data.table(cbind(draw=seq_len(P),null)),leave_one_out=loo,data=d)
}

sp_partial_fit<-function(x,y,covariates){
 ok<-is.finite(x)&is.finite(y)&apply(covariates,1,function(v)all(is.finite(v)));x<-x[ok];y<-y[ok];cc<-covariates[ok,,drop=FALSE];n<-length(x)
 if(n<3L)return(list(rho=NA_real_,n=n,df=NA_integer_))
 X<-cbind(1,apply(cc,2,rank));qq<-qr(X);df<-n-qq$rank
 rx<-qr.resid(qq,rank(x));ry<-qr.resid(qq,rank(y));den<-sqrt(sum(rx^2)*sum(ry^2))
 rr<-if(df>=3L&&den>1e-12)sum(rx*ry)/den else NA_real_
 list(rho=rr,n=n,df=df,x=x,y=y,covariates=cc,qr=qq,rx=rx,ry=ry)
}
sp_partial<-function(x,y,covariates,B=4999L,P=19999L){
 f<-sp_partial_fit(x,y,covariates);if(!is.finite(f$rho))return(list(summary=data.table(estimate=NA_real_,genes=f$n,residual_df=f$df,gene_bootstrap_low=NA_real_,gene_bootstrap_high=NA_real_,finite_bootstrap_draws=0L,raw_permutation_p=NA_real_),bootstrap=data.table(),permutation=data.table()))
 n<-f$n;idx<-sp_perm_index(n,P);E<-matrix(f$ry[idx],n,P);Q<-qr.Q(f$qr)[,seq_len(f$qr$rank),drop=FALSE];norm_y<-pmax(colSums(E^2)-colSums(crossprod(Q,E)^2),0);null<-as.numeric(crossprod(f$rx,E))/sqrt(sum(f$rx^2)*norm_y);null[norm_y<1e-20]<-NA_real_
 # Degenerate Freedman-Lane draws are recorded and omitted from the finite null.
 finite_null<-null[is.finite(null)];pv<-if(length(finite_null)>=P*.95)sp_p(f$rho,finite_null)else NA_real_
 bs<-vapply(seq_len(B),function(b){ii<-sample.int(n,n,replace=TRUE);sp_partial_fit(f$x[ii],f$y[ii],f$covariates[ii,,drop=FALSE])$rho},numeric(1));ci<-sp_ci(bs)
 list(summary=data.table(estimate=f$rho,genes=n,residual_df=f$df,gene_bootstrap_low=ci[1],gene_bootstrap_high=ci[2],finite_bootstrap_draws=sum(is.finite(bs)),raw_permutation_p=pv,finite_permutation_draws=length(finite_null)),bootstrap=data.table(draw=seq_len(B),estimate=bs),permutation=data.table(draw=seq_len(P),estimate=null))
}

sp_meta<-function(pairs,B=4999L){
 # pairs: gene, comparison, source_deltaZ, validation_deltaZ. Frangieh has3 contexts.
 comps<-c('RPE1','Adamson','Frangieh_Control','Frangieh_IFNg','Frangieh_Coculture');p<-as.data.table(pairs)[is.finite(source_deltaZ)&is.finite(validation_deltaZ)];gg<-sort(unique(p$gene));A<-O<-matrix(0,length(gg),5,dimnames=list(gg,comps))
 if(!length(gg))return(list(summary=data.table(estimate=NA_real_,gene_bootstrap_low=NA_real_,gene_bootstrap_high=NA_real_,finite_bootstrap_draws=0L),draws=data.table()))
 for(k in seq_along(comps)){d<-p[comparison==comps[k]];ii<-match(d$gene,gg);A[ii,k]<-sign(d$source_deltaZ)*d$validation_deltaZ;O[ii,k]<-1}
 means<-colSums(A)/colSums(O);fam<-c(RPE1=means[1],Adamson=means[2],Frangieh=mean(means[3:5]));point<-if(all(is.finite(fam)))mean(fam)else NA_real_
 w<-replicate(B,tabulate(sample.int(length(gg),length(gg),replace=TRUE),nbins=length(gg)));num<-crossprod(A,w);den<-crossprod(O,w);vv<-num/den;bs<-colMeans(rbind(vv[1,],vv[2,],colMeans(vv[3:5,,drop=FALSE])));ci<-sp_ci(bs)
 list(summary=data.table(estimate=point,gene_bootstrap_low=ci[1],gene_bootstrap_high=ci[2],finite_bootstrap_draws=sum(is.finite(bs)),RPE1_C=fam[1],Adamson_C=fam[2],Frangieh_mean_C=fam[3],leave_one_family_min=min(vapply(seq_len(3),function(k)mean(fam[-k]),numeric(1))),leave_one_family_max=max(vapply(seq_len(3),function(k)mean(fam[-k]),numeric(1)))),draws=data.table(draw=seq_len(B),C=bs))
}
