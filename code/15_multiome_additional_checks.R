#!/usr/bin/env Rscript
# Secondary depth, on-target transcript, negative-control and TSS sensitivities.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
a<-commandArgs(TRUE);pkg<-if(length(a)&&dir.exists(file.path(a[1],'inputs')))a[1]else'.';out<-if(length(a)>1&&!grepl('^--',a[2]))a[2]else file.path(pkg,'results/multiome_validation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(file.path(pkg,'inputs/multiome_count_inputs.rds'));g<-as.data.table(x$guide);tg<-as.data.table(x$targets);base_input<-if(is.null(x$original_rna_response))readRDS(file.path(pkg,'inputs/multiome.rds'))else NULL;R<-if(is.null(x$original_rna_response))base_input$rna_response else x$original_rna_response;A<-if(is.null(x$original_atac_response))base_input$atac_response else x$original_atac_response
set.seed(20261104);B<-500L
rms<-function(v)sqrt(mean(v*v));co<-function(v,w){d<-sqrt(sum(v*v)*sum(w*w));if(!is.finite(d)||d<1e-12)return(NA_real_);sum(v*w)/d}
logo<-function(M){if(ncol(M)<2)return(NA_real_);median(vapply(seq_len(ncol(M)),function(j)co(M[,j],rowMeans(M[,-j,drop=FALSE])),numeric(1)),na.rm=TRUE)}
fastlogo<-function(G,j){s<-G[j,j,drop=FALSE];v<-(rowSums(s)-diag(s))/sqrt(pmax(diag(s),1e-12)*pmax(sum(s)-2*rowSums(s)+diag(s),1e-12));median(v,na.rm=TRUE)}
index_target<-function(gene)which(g$context=='ATRT_human'&g$guide_gene==gene)
summarize<-function(RR,AA,window='all_annotated'){
 z<-rbindlist(lapply(which(tg$context=='ATRT_human'),function(i){t<-tg[i];jj<-index_target(t$guide_gene);data.table(gene=t$guide_gene,n_guides=length(jj),n_cells=t$n_cells,count_quality_gate=t$count_quality_gate,rna_cosine=logo(RR[,jj,drop=FALSE]),atac_cosine=logo(AA[,jj,drop=FALSE]),rna_movement=rms(rowMeans(RR[,jj,drop=FALSE])),atac_movement=rms(rowMeans(AA[,jj,drop=FALSE])))}))
 z[,passes:=count_quality_gate&is.finite(rna_cosine)&is.finite(atac_cosine)&rna_cosine>.1&atac_cosine>.1];z[,window:=window]
 z[,role:='insufficient_support'];if(any(z$passes)){rc<-quantile(z$rna_movement[z$passes],.75);ac<-quantile(z$atac_movement[z$passes],.75);z[passes==TRUE,role:=fifelse(rna_movement>=rc&atac_movement>=ac,'coupled_mover',fifelse(rna_movement>=rc,'rna_dominant',fifelse(atac_movement>=ac,'atac_dominant','weak_or_contextual')))]};z
}
primary<-summarize(R,A)
# 1. Alternative restrictions of existing peak-to-gene annotations.
tss<-rbindlist(c(list(primary),lapply(names(x$atac_tss_responses),function(w)summarize(R,x$atac_tss_responses[[w]],paste0('TSS_',w)))))
tss[,atac_genes:=c(all_annotated=nrow(A),TSS_2000=nrow(x$atac_tss_responses[['2000']]),TSS_10000=nrow(x$atac_tss_responses[['10000']]))[window]]
fwrite(tss,file.path(out,'tss_task_sensitivity.tsv'),sep='\t');fwrite(tss[,.(eligible=sum(passes),coupled=sum(role=='coupled_mover'),primary_eligible_retained=sum(passes&gene%in%primary[passes==TRUE,gene]),atac_genes=unique(atac_genes)),by=window],file.path(out,'tss_summary.tsv'),sep='\t')
# 2. Actual non-targeting-guide construct null, excluding each null group from its reference.
nt<-which(g$context=='ATRT_human'&g$is_control&g$n_cells>=15)
ntnull<-list()
for(k in c(2L,3L)){
 groups<-combn(nt,k,simplify=FALSE)
 for(jj in groups){ref<-setdiff(nt,jj);rr<-as.numeric(R[,ref,drop=FALSE]%*%g$n_cells[ref]/sum(g$n_cells[ref]));ar<-as.numeric(A[,ref,drop=FALSE]%*%g$n_cells[ref]/sum(g$n_cells[ref]));rc<-logo(sweep(R[,jj,drop=FALSE],1,rr,'-'));ac<-logo(sweep(A[,jj,drop=FALSE],1,ar,'-'));ntnull[[length(ntnull)+1L]]<-data.table(n_guides=k,guide_ids=paste(g$guide_id[jj],collapse=';'),rna_cosine=rc,atac_cosine=ac,joint=min(rc,ac))}
}
ntnull<-rbindlist(ntnull);fwrite(ntnull,file.path(out,'non_targeting_guide_null.tsv'),sep='\t')
nteval<-primary[count_quality_gate==TRUE];nteval[,nt_construct_p:=vapply(seq_len(.N),function(i){z<-ntnull[n_guides==nteval$n_guides[i],joint];(1+sum(z>=min(nteval$rna_cosine[i],nteval$atac_cosine[i])))/(1+length(z))},numeric(1))];nteval[,nt_construct_fdr:=p.adjust(nt_construct_p,'BH')]
fwrite(nteval,file.path(out,'non_targeting_target_comparison.tsv'),sep='\t');fwrite(ntnull[,.(groups=.N,joint_null_95=quantile(joint,.95),joint_null_median=median(joint),false_agreement_fraction=mean(rna_cosine>.1&atac_cosine>.1)),by=n_guides],file.path(out,'non_targeting_null_summary.tsv'),sep='\t')
# 3. A target-RNA depletion proxy; not a direct guide-efficiency measurement.
om<-as.data.table(x$on_target_cell_metadata);nnt<-length(x$control_rna_library)
w<-replicate(B,tabulate(sample.int(nnt,nnt,replace=TRUE),nbins=nnt));cb<-as.matrix(x$control_target_counts%*%w);lb<-as.numeric(crossprod(x$control_rna_library,w));cpoint<-Matrix::rowSums(x$control_target_counts);lpoint<-sum(x$control_rna_library)
ot<-rbindlist(lapply(which(g$context=='ATRT_human'&!g$is_control&g$n_cells>=15),function(j){
 d<-om[guide_index==j];i<-match(g$guide_gene[j],x$target_genes);if(is.na(i)||!nrow(d))return(NULL)
 lp<-log2(((sum(d$target_rna_count)+.5)/(sum(d$rna_library)+1))/((cpoint[i]+.5)/(lpoint+1)))
 vb<-vapply(seq_len(B),function(b){ii<-sample.int(nrow(d),nrow(d),replace=TRUE);log2(((sum(d$target_rna_count[ii])+.5)/(sum(d$rna_library[ii])+1))/((cb[i,b]+.5)/(lb[b]+1)))},numeric(1))
 data.table(guide_id=g$guide_id[j],gene=g$guide_gene[j],n_cells=nrow(d),control_target_counts=cpoint[i],control_target_cpm=1e6*cpoint[i]/lpoint,estimable_control=cpoint[i]>=20,target_rna_log2_ratio=lp,cell_bootstrap_low=quantile(vb,.025),cell_bootstrap_high=quantile(vb,.975),rna_movement=g$rna_movement[j],atac_movement=g$atac_movement[j],mean_rna_library=g$mean_rna_library[j],mean_atac_library=g$mean_atac_library[j],primary_target=g$guide_gene[j]%in%primary[passes==TRUE,gene])
}))
ot[,target_in_rna_annotation:=x$target_gene_feature_count[match(gene,x$target_genes)]>0]
ot[,proxy_evaluation_status:=fifelse(!target_in_rna_annotation,'not_measured_symbol_absent',fifelse(!estimable_control,'insufficient_control_target_counts','available_transcript_proxy'))]
fwrite(ot,file.path(out,'target_rna_depletion_by_guide.tsv'),sep='\t')
ots<-ot[,.(guides=.N,genes=uniqueN(gene),estimable=sum(estimable_control),lower_rna_point=sum(estimable_control&target_rna_log2_ratio<0),upper_interval_below_zero=sum(estimable_control&cell_bootstrap_high<0),median_log2_ratio=median(target_rna_log2_ratio[estimable_control])),by=primary_target];fwrite(ots,file.path(out,'target_rna_depletion_summary.tsv'),sep='\t')
otcor<-rbindlist(lapply(c('all_supported_guides','primary_guides'),function(scope){z<-ot[estimable_control==TRUE&(scope=='all_supported_guides'|primary_target)];rbindlist(lapply(c('rna_movement','atac_movement','n_cells','mean_rna_library','mean_atac_library'),function(v)data.table(scope=scope,n=nrow(z),variable=v,spearman=cor(z$target_rna_log2_ratio,z[[v]],method='spearman'))))}));fwrite(otcor,file.path(out,'target_rna_depletion_correlations.tsv'),sep='\t')
# Amplitude and depth matched label-null (orientation diagnostic, not a biological test).
pool<-which(g$context=='ATRT_human'&!g$is_control&g$n_cells>=15&g$nuisance_separation_check)
cov<-scale(cbind(log1p(g$n_cells[pool]),log1p(g$mean_rna_library[pool]),log1p(g$mean_atac_library[pool]),log1p(g$rna_movement[pool]),log1p(g$atac_movement[pool])))
GR<-crossprod(R);GA<-crossprod(A);set.seed(20261105)
matched<-rbindlist(lapply(primary[count_quality_gate==TRUE,gene],function(gene_nm){jj<-index_target(gene_nm);neighbors<-lapply(jj,function(j){q<-match(j,pool);ds<-rowSums(sweep(cov,2,cov[q,],'-')^2);ii<-order(ds);ii<-ii[g$guide_gene[pool[ii]]!=gene_nm];pool[head(ii,20)]});sm<-replicate(999,{sel<-integer();for(z in neighbors){z<-z[!g$guide_gene[z]%in%g$guide_gene[sel]];sel<-c(sel,sample(z,1))};c(fastlogo(GR,sel),fastlogo(GA,sel))});obs<-primary[gene==gene_nm];data.table(gene=gene_nm,rna_cosine=obs$rna_cosine,atac_cosine=obs$atac_cosine,p=(1+sum(pmin(sm[1,],sm[2,])>=min(obs$rna_cosine,obs$atac_cosine)))/1000,null95=quantile(pmin(sm[1,],sm[2,]),.95),null_pass=mean(sm[1,]>.1&sm[2,]>.1))}))
matched[,fdr:=p.adjust(p,'BH')];fwrite(matched,file.path(out,'depth_amplitude_matched_label_null.tsv'),sep='\t')
cat('Completed TSS, NT construct null, target-RNA depletion and matched-label checks.\n');print(fread(file.path(out,'tss_summary.tsv')));print(fread(file.path(out,'non_targeting_null_summary.tsv')));print(ots)
# 4. Equal-cell counts: every guide contributes m cells sampled without replacement.
set.seed(20261106)
cell_by_guide<-split(seq_along(x$cell_guide),x$cell_guide)
pbulk<-function(M,L,ii)log1p(as.numeric(Matrix::rowSums(M[,ii,drop=FALSE]))/max(sum(L[ii]),1)*1e6)
# Standard control anchor uses the original cell-weighted guide log-CPM rule.
ctrlR<-matrix(0,nrow(R),length(nt));ctrlA<-matrix(0,nrow(A),length(nt));for(j in seq_along(nt)){ii<-cell_by_guide[[as.character(nt[j])]];ctrlR[,j]<-pbulk(x$rna_counts,x$rna_library,ii);ctrlA[,j]<-pbulk(x$atac_counts,x$atac_library,ii)}
refR<-as.numeric(ctrlR%*%g$n_cells[nt]/sum(g$n_cells[nt]));refA<-as.numeric(ctrlA%*%g$n_cells[nt]/sum(g$n_cells[nt]))
rs<-list();null_bal<-list()
for(m in c(15L,30L,50L)){
 genes<-tg[count_quality_gate==TRUE&minimum_guide_cells>=m,guide_gene]
 # 499 balanced negative-control groups for each guide count, holding all selected null guides out of the reference.
 for(k in c(2L,3L))for(b in seq_len(499)){
  jj<-sample(nt,k,replace=FALSE);other<-match(setdiff(nt,jj),nt);cr<-as.numeric(ctrlR[,other,drop=FALSE]%*%g$n_cells[nt[other]]/sum(g$n_cells[nt[other]]));ca<-as.numeric(ctrlA[,other,drop=FALSE]%*%g$n_cells[nt[other]]/sum(g$n_cells[nt[other]]));vR<-vA<-NULL
  for(j in jj){ii<-sample(cell_by_guide[[as.character(j)]],m,replace=FALSE);vR<-cbind(vR,pbulk(x$rna_counts,x$rna_library,ii)-cr);vA<-cbind(vA,pbulk(x$atac_counts,x$atac_library,ii)-ca)}
  nr<-logo(vR);na<-logo(vA);null_bal[[length(null_bal)+1L]]<-data.table(m=m,n_guides=k,replicate=b,rna_cosine=nr,atac_cosine=na,joint=min(nr,na))
 }
 for(b in seq_len(100))for(gene in genes){
  jj<-index_target(gene);vR<-vA<-NULL;Lr<-La<-numeric()
  for(j in jj){ii<-sample(cell_by_guide[[as.character(j)]],m,replace=FALSE);vR<-cbind(vR,pbulk(x$rna_counts,x$rna_library,ii)-refR);vA<-cbind(vA,pbulk(x$atac_counts,x$atac_library,ii)-refA);Lr<-c(Lr,sum(x$rna_library[ii]));La<-c(La,sum(x$atac_library[ii]))}
  rr<-logo(vR);aa<-logo(vA);rs[[length(rs)+1L]]<-data.table(m=m,replicate=b,gene=gene,n_guides=length(jj),rna_cosine=rr,atac_cosine=aa,passes=rr>.1&aa>.1,rna_movement=rms(rowMeans(vR)),atac_movement=rms(rowMeans(vA)),min_rna_molecules=min(Lr),min_atac_molecules=min(La))
 }
 cat('Equal-cell m',m,'targets',length(genes),'completed\n')
}
rs<-rbindlist(rs);nb<-rbindlist(null_bal);fwrite(rs,file.path(out,'equal_cell_replicates.tsv'),sep='\t');fwrite(nb,file.path(out,'equal_cell_negative_control_null.tsv'),sep='\t')
# Roles are evaluated within the same primary ten targets when supported by a depth setting.
rs[,primary_target:=gene%in%primary[passes==TRUE,gene]];rs[,role:='not_primary_or_failed']
rs[primary_target==TRUE&passes==TRUE,role:={rh<-rna_movement>=quantile(rna_movement,.75);ah<-atac_movement>=quantile(atac_movement,.75);fifelse(rh&ah,'coupled_mover',fifelse(rh,'rna_dominant',fifelse(ah,'atac_dominant','weak_or_contextual')))},by=.(m,replicate)]
fr<-rs[,.(resamples=.N,agreement_pass_fraction=mean(passes),median_rna_cosine=median(rna_cosine),median_atac_cosine=median(atac_cosine),coupled_fraction=mean(role=='coupled_mover'),primary_target=unique(primary_target)),by=.(m,gene)]
fwrite(fr,file.path(out,'equal_cell_target_stability.tsv'),sep='\t')
summary<-rs[,.(eligible=sum(passes),primary_retained=sum(passes&primary_target)),by=.(m,replicate)][,.(target_draws=.N,median_eligible=median(eligible),eligible_min=min(eligible),eligible_max=max(eligible),median_primary_retained=median(primary_retained)),by=m]
fwrite(summary,file.path(out,'equal_cell_summary.tsv'),sep='\t')
nullsum<-nb[,.(draws=.N,null95=quantile(joint,.95),null_median=median(joint),pass_fraction=mean(rna_cosine>.1&atac_cosine>.1)),by=.(m,n_guides)];fwrite(nullsum,file.path(out,'equal_cell_null_summary.tsv'),sep='\t')
# One seeded balanced draw per target is compared with the matched-cell-count null.
first<-rs[replicate==1];first[,nt_p:=vapply(seq_len(.N),function(i){zz<-nb[m==first$m[i]&n_guides==first$n_guides[i],joint];(1+sum(zz>=min(first$rna_cosine[i],first$atac_cosine[i])))/(1+length(zz))},numeric(1))];first[,nt_fdr:=p.adjust(nt_p,'BH'),by=m];fwrite(first,file.path(out,'equal_cell_first_draw_null_test.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));print(summary);print(nullsum);cat('Multiome validation completed.\n')


# Additional count- and library-matched contrasts with independent controls per guide.
suppressPackageStartupMessages({library(data.table);library(Matrix);library(BiocNeighbors)})
a<-commandArgs(TRUE);pkg<-if(length(a)&&dir.exists(file.path(a[1],'inputs')))a[1]else'.';out<-if(length(a)>1&&!grepl('^--',a[2]))a[2]else file.path(pkg,'results/multiome_validation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(file.path(pkg,'inputs/multiome_count_inputs.rds'));g<-as.data.table(x$guide);tg<-as.data.table(x$targets);primary<-fread(file.path(out,'tss_task_sensitivity.tsv'))[window=='all_annotated'&passes==TRUE,gene]
set.seed(20261107)
ntg<-which(g$context=='ATRT_human'&g$is_control&g$n_cells>=15);ntcells<-which(x$cell_guide%in%ntg)
Z<-cbind(log1p(x$rna_library),log1p(x$atac_library));mu<-colMeans(Z[ntcells,,drop=FALSE]);sig<-apply(Z[ntcells,,drop=FALSE],2,sd);Z<-sweep(sweep(Z,2,mu,'-'),2,sig,'/')
knn<-BiocNeighbors::queryKNN(Z[ntcells,,drop=FALSE],Z,k=50,BNPARAM=KmknnParam())$index
byguide<-split(seq_along(x$cell_guide),x$cell_guide)
matchNT<-function(ii,blocked=integer()){
 picked<-integer();for(j in ii){pool<-ntcells[knn[j,]];pool<-pool[!pool%in%c(blocked,picked)];if(!length(pool))pool<-setdiff(ntcells,c(blocked,picked));picked<-c(picked,pool[sample.int(length(pool),1)])};picked
}
pb<-function(M,L,ii)log1p(as.numeric(Matrix::rowSums(M[,ii,drop=FALSE]))/max(sum(L[ii]),1)*1e6)
cos<-function(v,w){dd<-sqrt(sum(v*v)*sum(w*w));if(dd<1e-12)return(NA_real_);sum(v*w)/dd}
lg<-function(M)median(vapply(seq_len(ncol(M)),function(j)cos(M[,j],rowMeans(M[,-j,drop=FALSE])),numeric(1)),na.rm=TRUE)
contrast<-function(sets,nt_source=FALSE){vr<-va<-NULL;blocked<-if(nt_source)unlist(sets)else integer();dif<-numeric()
 for(ii in sets){ct<-matchNT(ii,blocked);blocked<-c(blocked,ct);vr<-cbind(vr,pb(x$rna_counts,x$rna_library,ii)-pb(x$rna_counts,x$rna_library,ct));va<-cbind(va,pb(x$atac_counts,x$atac_library,ii)-pb(x$atac_counts,x$atac_library,ct));dif<-c(dif,mean(abs(Z[ii,]-Z[ct,])))}
 c(rna=lg(vr),atac=lg(va),rna_movement=sqrt(mean(rowMeans(vr)^2)),atac_movement=sqrt(mean(rowMeans(va)^2)),mean_scaled_library_mismatch=mean(dif))
}
rr<-nn<-list()
for(m in c(15L,30L,50L)){
 for(k in c(2L,3L))for(b in seq_len(199)){
  cells<-sample(ntcells,m*k,replace=FALSE);sets<-split(cells,rep(seq_len(k),each=m));v<-contrast(sets,TRUE);nn[[length(nn)+1L]]<-data.table(m=m,n_guides=k,replicate=b,rna_cosine=v[1],atac_cosine=v[2],joint=min(v[1:2]),library_mismatch=v[5])
 }
 genes<-tg[count_quality_gate==TRUE&minimum_guide_cells>=m,guide_gene]
 for(b in seq_len(50))for(gene in genes){jj<-which(g$context=='ATRT_human'&g$guide_gene==gene);sets<-lapply(jj,function(j)sample(byguide[[as.character(j)]],m,replace=FALSE));v<-contrast(sets);rr[[length(rr)+1L]]<-data.table(m=m,replicate=b,gene=gene,n_guides=length(jj),rna_cosine=v[1],atac_cosine=v[2],passes=v[1]>.1&v[2]>.1,rna_movement=v[3],atac_movement=v[4],library_mismatch=v[5])}
 cat('Independent-reference m',m,'completed\n')
}
r<-rbindlist(rr);n<-rbindlist(nn);fwrite(r,file.path(out,'independent_reference_replicates.tsv'),sep='\t');fwrite(n,file.path(out,'independent_reference_negative_control_null.tsv'),sep='\t')
r[,primary_target:=gene%in%primary];sm<-r[,.(eligible=sum(passes),primary_retained=sum(passes&primary_target)),by=.(m,replicate)][,.(median_eligible=median(eligible),minimum_eligible=min(eligible),maximum_eligible=max(eligible),median_primary_retained=median(primary_retained)),by=m]
fwrite(sm,file.path(out,'independent_reference_summary.tsv'),sep='\t')
fwrite(r[,.(pass_fraction=mean(passes),median_rna_cosine=median(rna_cosine),median_atac_cosine=median(atac_cosine),mean_library_mismatch=mean(library_mismatch),primary_target=unique(primary_target)),by=.(m,gene)],file.path(out,'independent_reference_target_stability.tsv'),sep='\t')
ns<-n[,.(draws=.N,joint_null95=quantile(joint,.95),pass_fraction=mean(rna_cosine>.1&atac_cosine>.1),mean_library_mismatch=mean(library_mismatch)),by=.(m,n_guides)];fwrite(ns,file.path(out,'independent_reference_null_summary.tsv'),sep='\t')
f<-r[replicate==1];f[,nt_p:=vapply(seq_len(.N),function(i){z<-n[m==f$m[i]&n_guides==f$n_guides[i],joint];(1+sum(z>=min(f$rna_cosine[i],f$atac_cosine[i])))/(1+length(z))},numeric(1))];f[,nt_fdr:=p.adjust(nt_p,'BH'),by=m];fwrite(f,file.path(out,'independent_reference_first_draw_null_test.tsv'),sep='\t');print(sm);print(ns)


suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);pkg<-if(length(a)&&dir.exists(file.path(a[1],'inputs')))a[1]else'.';out<-if(length(a)>1&&!grepl('^--',a[2]))a[2]else file.path(pkg,'results/multiome_validation')
e<-as.data.table(readRDS(file.path(pkg,'reference/published_numerical_reference.rds'))$evidence);nt<-fread(file.path(out,'non_targeting_target_comparison.tsv'));ntpass<-nt[nt_construct_fdr<.05,gene]
z<-copy(e);z[,primary_class:=evidence_class];z[,multiperturb_eligible:=multiperturb_eligible&gene%in%ntpass]
z[,datasets_with_causal_support:=as.integer(replogle_eligible)+as.integer(norman_single_eligible)+as.integer(multiperturb_eligible)]
z[,independent_contexts_with_causal_support:=as.integer(replogle_eligible|norman_single_eligible)+as.integer(multiperturb_eligible)]
z[,causal_regulator:=replogle_eligible|norman_single_eligible|multiperturb_eligible]
z[,validated_control_node_candidate:=causal_regulator&independent_contexts_with_causal_support>=2&combination_consistent&path_specific]
cut<-quantile(z$observed_combinations[z$observed_combinations>0],.9)
z[,evidence_class:=fifelse(validated_control_node_candidate,'validated_control_node_candidate',fifelse(causal_regulator,'causal_regulator',fifelse(observed_combinations>0&observed_combinations>=cut,'association_hub',fifelse(genomewide_screened&is.finite(genomewide_effect_fdr)&genomewide_effect_fdr<.05,'state_marker','insufficient_causal_evidence'))))]
summary<-z[,.N,by=evidence_class];summary[,sensitivity:='MultiPerturb requires NT-construct-null BH q<0.05'];fwrite(summary,file.path(out,'nt_required_evidence_class_counts.tsv'),sep='\t')
changed<-z[evidence_class!=primary_class,.(gene,primary_class,sensitivity_class=evidence_class,genomewide_screened,genomewide_effect_fdr)];fwrite(changed,file.path(out,'nt_required_evidence_changed_genes.tsv'),sep='\t')
fwrite(data.table(metric=c('MPS_eligible','causal_eligible','two_context_eligible','control_candidates'),value=c(sum(z$multiperturb_eligible),sum(z$causal_regulator),sum(z$independent_contexts_with_causal_support>=2),sum(z$validated_control_node_candidate))),file.path(out,'nt_required_evidence_summary.tsv'),sep='\t')
print(summary);print(changed)


#!/usr/bin/env Rscript
# Four display states without changing any primary gene class or eligibility flag.
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(TRUE);pkg<-if(length(a)&&dir.exists(file.path(a[1],'inputs')))a[1]else'.';out<-if(length(a)>1&&!grepl('^--',a[2]))a[2]else file.path(pkg,'results/multiome_validation');dir.create(out,recursive=TRUE,showWarnings=FALSE)
r<-readRDS(file.path(pkg,'reference/published_numerical_reference.rds'));aux<-readRDS(file.path(pkg,'inputs/auxiliary.rds'));x<-readRDS(file.path(pkg,'inputs/multiome_count_inputs.rds'))
# Recompute the source's component gates using the unchanged Replogle code prefix.
s<-readLines(file.path(pkg,'code/04_task_registry.R'));s<-s[seq_len(which(grepl('^# Norman:',s))[1]-1L)];s<-paste(s,collapse='\n');s<-gsub('source("code/common.R")','source(file.path(pkg,"code/common.R"))',s,fixed=TRUE);s<-gsub('out_dir <- "results/task_registry"','out_dir <- out',s,fixed=TRUE);s<-gsub('"inputs/replogle_scores.rds"','file.path(pkg,"inputs/replogle_scores.rds")',s,fixed=TRUE);eval(parse(text=s))
rep<-rep_summary[,.(gene,rep_cells=target_cell_count,rep_gems=replicate_groups,rep_nuisance=nuisance_separation_check,rep_met=passes_task_qc)]
rep[,rep_assayed:=TRUE];rep[,rep_evaluable:=rep_cells>=100&rep_gems>=3&rep_nuisance]
n<-as.data.table(aux$norman_index);q<-as.data.table(aux$norman_qc);n<-merge(n,q,by=c('condition_id','raw_guide_identity','n_cells','condition_type'))
ns<-n[condition_type=='single'&!is.na(single_gene),.(norm_assayed=TRUE,norm_evaluable=any(n_cells>=100&n_gem_groups>=3&nuisance_separation_check)),by=.(gene=single_gene)]
ms<-as.data.table(x$targets)[context=='ATRT_human',.(gene=guide_gene,multi_assayed=TRUE,multi_evaluable=count_quality_gate)]
e<-as.data.table(r$evidence);orig<-e[,.(gene,evidence_class,validated_control_node_candidate)]
e<-merge(e,rep,by='gene',all.x=TRUE);e<-merge(e,ns,by='gene',all.x=TRUE);e<-merge(e,ms,by='gene',all.x=TRUE)
for(v in c('rep_assayed','rep_evaluable','norm_assayed','norm_evaluable','multi_assayed','multi_evaluable'))e[is.na(get(v)),(v):=FALSE]
e[,`:=`(k562_assayed=rep_assayed|norm_assayed,k562_evaluable=rep_evaluable|norm_evaluable,datasets_assayed=as.integer(rep_assayed)+as.integer(norm_assayed)+as.integer(multi_assayed),datasets_evaluable=as.integer(rep_evaluable)+as.integer(norm_evaluable)+as.integer(multi_evaluable))]
state<-function(met,measured,evaluable)fifelse(met,'Met',fifelse(!measured,'Not measured',fifelse(!evaluable,'Insufficient support','Not met')))
e[,`:=`(causal_state=state(causal_regulator,rep_assayed|norm_assayed|multi_assayed,rep_evaluable|norm_evaluable|multi_evaluable),multiple_datasets_state=state(datasets_with_causal_support>=2,datasets_assayed>=2,datasets_evaluable>=2),multiple_contexts_state=state(independent_contexts_with_causal_support>=2,k562_assayed&multi_assayed,k562_evaluable&multi_evaluable),specificity_state=state(path_specific,k562_assayed,replogle_eligible|norman_single_eligible),combination_state=state(combination_consistent,observed_combinations>0,observed_combinations>0))]
e[,candidate_state:=vapply(seq_len(.N),function(i){ss<-c(multiple_contexts_state[i],specificity_state[i],combination_state[i]);if(validated_control_node_candidate[i])'Met'else if(any(ss=='Not measured'))'Not measured'else if(any(ss=='Insufficient support'))'Insufficient support'else'Not met'},character(1))]
e[,evaluation_note:='Not met denotes a measured numerical criterion that was not passed; it is not proof of a biological absence. Insufficient support includes inadequate counts/guides or unresolved nuisance separation.']
keep<-c('gene','evidence_class','causal_state','multiple_datasets_state','multiple_contexts_state','specificity_state','combination_state','candidate_state','rep_assayed','rep_evaluable','norm_assayed','norm_evaluable','multi_assayed','multi_evaluable','datasets_assayed','datasets_evaluable','observed_combinations','evaluation_note')
fwrite(e[,..keep],file.path(out,'gene_criterion_four_state.tsv'),sep='\t')
m<-melt(e[,.(gene,`Multiple datasets`=multiple_datasets_state,`Multiple contexts`=multiple_contexts_state,`Path specific`=specificity_state,`Pair consistent`=combination_state,`Control candidate`=candidate_state)],id.vars='gene',variable.name='criterion',value.name='state');fwrite(m,file.path(out,'gene_criterion_four_state_long.tsv.gz'),sep='\t')
ct<-m[,.N,by=.(criterion,state)];fwrite(ct,file.path(out,'gene_criterion_four_state_counts.tsv'),sep='\t')
z<-merge(orig,e[,.(gene,evidence_class,validated_control_node_candidate)],by='gene');stopifnot(all(z$evidence_class.x==z$evidence_class.y),all(z$validated_control_node_candidate.x==z$validated_control_node_candidate.y))
print(ct);cat('Primary evidence classes and candidate flags unchanged.\n')
