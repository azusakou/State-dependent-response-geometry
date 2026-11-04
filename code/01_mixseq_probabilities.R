suppressPackageStartupMessages(library(data.table))
out<-'results/mixseq_probabilities';dir.create(out,recursive=TRUE,showWarnings=FALSE)
s<-fread('inputs/mixseq_program_cell_scores.tsv.gz');aux<-readRDS('inputs/auxiliary.rds')
original<-readRDS('reference/published_numerical_reference.rds')$mixseq_probabilities
id_map<-aux$mixseq_id_map
rank<-c(dmso=0,untreated=1,vehicle=2,control=3,sglacz=0,sgor2j2=1);s[,control_rank:=unname(rank[tolower(perturbation)])]
keys<-c('dataset_id','experiment','time_hours','task','cell_line');endpoints<-c('apoptosis','proliferation','stress');rows<-list();s[,context_id:=.GRP,by=keys]
for(ctx in unique(s$context_id)){
 a<-s[context_id==ctx];controls<-a[is.finite(control_rank)];if(!nrow(controls))next;controls<-controls[control_rank==min(control_rank)];if(nrow(controls)<20)next
 for(pt in a[is.na(control_rank),unique(perturbation)]){
 b<-a[perturbation==pt];if(nrow(b)<20)next
 z<-b[1,c(keys,'perturbation'),with=FALSE];z[,`:=`(baseline=controls$perturbation[1],baseline_cells=nrow(controls),treated_cells=nrow(b))]
 for(e in endpoints){v<-controls[[paste0('score_',e)]];v<-v[is.finite(v)];th<-median(v)+1.4826*median(abs(v-median(v)));t<-b[[paste0('score_',e)]];t<-t[is.finite(t)];successes<-sum(t>th);p<-(successes+1)/(length(t)+2);z[,(paste0(e,'_threshold')):=th];z[,(paste0('p_',e)):=p];z[,(paste0('p_',e,'_se')):=sqrt(p*(1-p)/(length(t)+3))];z[,(paste0('n_',e)):=length(t)];z[,(paste0('successes_',e)):=successes]}
 rows[[length(rows)+1L]]<-z
 }
}
r<-rbindlist(rows);r<-merge(r,id_map,by=c(keys,'perturbation','baseline'),sort=FALSE)
fwrite(r,file.path(out,'pairs.tsv'),sep='\t')
cmp<-merge(r,original,by='pair_id',suffixes=c('_new','_original'));fields<-c('baseline_cells','treated_cells',paste0('p_',endpoints),paste0('p_',endpoints,'_se'))
checks<-rbindlist(lapply(fields,function(nm){d<-cmp[[paste0(nm,'_new')]]-cmp[[paste0(nm,'_original')]];data.table(field=nm,n=nrow(cmp),max_absolute_difference=max(abs(d),na.rm=TRUE),n_different=sum(abs(d)>1e-10,na.rm=TRUE))}))
fwrite(r,file.path(out,'mixseq_probability_sufficient_counts.tsv'),sep='\t');fwrite(checks,file.path(out,'mixseq_probability_replication.tsv'),sep='\t');print(checks);cat('Rows reconstructed',nrow(r),'original',nrow(original),'matched',nrow(cmp),'\n');stopifnot(nrow(r)==nrow(original),nrow(cmp)==nrow(original),all(checks$max_absolute_difference<1e-10))
