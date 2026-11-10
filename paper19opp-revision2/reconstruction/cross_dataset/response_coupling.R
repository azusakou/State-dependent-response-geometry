#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
.libPaths(c(file.path(root,'.Rlib'),.libPaths()));suppressPackageStartupMessages(library(data.table))
pairs<-list(c('stress','proliferation'),c('stress','apoptosis'),c('apoptosis','proliferation'))
robust<-function(a,xaxis,yaxis,half_x,half_y,magnitude_adjust=FALSE,anchor=FALSE){
 xname<-paste0('response_',xaxis,'_',half_x);yname<-paste0('response_',yaxis,'_',half_y)
 a<-copy(a);a[,x:=get(xname)];a[,y:=get(yname)]
 a[,`:=`(log_n=log(n_treated_A+n_treated_B),library_A=mean_log_library_treated_A,library_B=mean_log_library_treated_B)]
 a<-a[context%in%a[,.N,by=context][N>=2,context]]
 if(magnitude_adjust)a<-a[is.finite(independent_context_magnitude)]
 a[,x:=as.numeric(scale(x))]
 form<-if(magnitude_adjust)y~x+log_n+library_A+library_B+independent_context_magnitude+factor(context)else y~x+log_n+library_A+library_B+factor(context)
 f<-lm(form,data=a);keep<-names(coef(f))[is.finite(coef(f))]
 xx<-model.matrix(f)[,keep,drop=FALSE];u<-residuals(f);ii<-match('x',keep)
 cluster<-a$perturbation
 if(anchor)cluster<-vapply(strsplit(cluster,'_',fixed=TRUE),function(g)sort(g[!grepl('^NegCtrl',g)])[1],character(1))
 score<-rowsum(xx*u,cluster,reorder=FALSE);g<-nrow(score);n<-nrow(xx);k<-ncol(xx)
 bread<-solve(crossprod(xx));vv<-bread%*%crossprod(score)%*%bread*g/(g-1)*(n-1)/(n-k)
 se<-sqrt(vv[ii,ii]);b<-coef(f)['x'];ci<-b+c(-1,1)*qt(.975,g-1)*se
 data.table(n_units=n,n_contexts=uniqueN(a$context),intervention_clusters=g,beta_per_sd_x=as.numeric(b),cluster_se=se,ci_low=ci[1],ci_high=ci[2],p=2*pt(-abs(b/se),g-1),raw_pearson=cor(a[[xname]],a[[yname]]),effect='score units of independent-half Y per within-source SD of X',adjustment=if(magnitude_adjust)'context FE, counts, libraries, leave-context-out intervention response magnitude'else'context FE, counts, libraries')
}
results<-list()
for(rep in 1:10){
 d<-fread(file.path(out,paste0('analysis_units_partition',rep,'.tsv.gz')))
 responses<-grep('^response_.*_[AB]$',names(d),value=TRUE)
 # Magnitude uses OTHER contexts of the same recorded intervention; it excludes
 # all observations, response anchors and control pools in the query context.
 mag<-d[,.(observed_magnitude=sqrt(mean(as.matrix(.SD)^2))),by=.(dataset,context,perturbation),.SDcols=responses]
 mag[,`:=`(intervention_sum=sum(observed_magnitude),intervention_n=.N),by=.(dataset,perturbation)]
 mag[,independent_context_magnitude:=(intervention_sum-observed_magnitude)/pmax(intervention_n-1L,1L)]
 mag[intervention_n<2L,independent_context_magnitude:=NA_real_]
 d<-merge(d,mag[,.(dataset,context,perturbation,independent_context_magnitude)],by=c('dataset','context','perturbation'))
 for(s in unique(d$dataset))for(pair in pairs)for(swap in c(FALSE,TRUE))for(adj in c(FALSE,TRUE)){
  if(rep>1L&&(swap||adj))next
  z<-robust(d[dataset==s],pair[1],pair[2],if(swap)'B'else'A',if(swap)'A'else'B',adj)
  z[,`:=`(partition_id=rep,dataset=s,x_axis=pair[1],y_axis=pair[2],half_direction=if(swap)'B_to_A'else'A_to_B',magnitude_adjustment=adj,cluster_definition='complete_condition')]
  results[[length(results)+1L]]<-z
  if(rep==1L&&s=='Norman'&&!swap&&!adj){
   z<-robust(d[dataset==s],pair[1],pair[2],'A','B',FALSE,TRUE)
   z[,`:=`(partition_id=rep,dataset=s,x_axis=pair[1],y_axis=pair[2],half_direction='A_to_B',magnitude_adjustment=FALSE,cluster_definition='deterministic_anchor_gene')]
   results[[length(results)+1L]]<-z
  }
 }
}
r<-rbindlist(results)
# Nine primary comparisons. Swapped-half, adjusted-magnitude and alternative
# partition results are sensitivities, with their own explicitly separate BH sets.
r[,q_BH:=p.adjust(p,'BH'),by=.(partition_id,half_direction,magnitude_adjustment,cluster_definition)]
fwrite(r,file.path(out,'cross_program_coupling.tsv'),sep='\t')
print(r[partition_id==1&half_direction=='A_to_B'&magnitude_adjustment==FALSE,.(dataset,x_axis,y_axis,n_units,intervention_clusters,beta_per_sd_x,ci_low,ci_high,q_BH)])
