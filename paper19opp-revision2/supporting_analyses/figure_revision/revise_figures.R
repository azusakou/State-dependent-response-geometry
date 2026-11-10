#!/usr/bin/env Rscript
# Reconstruct revised figures from preserved source data and explicit specificity tests.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(Matrix)})
args<-commandArgs(TRUE)
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
if(length(file_arg)!=1)stop('Run this script with Rscript.')
script_dir<-dirname(normalizePath(sub('^--file=','',file_arg)))
root<-if(length(args)>=1)normalizePath(args[1],mustWork=TRUE)else file.path(script_dir,'source_inputs')
target<-if(length(args)>=2)args[2]else file.path(script_dir,'reproduced')
dir.create(target,recursive=TRUE,showWarnings=FALSE)
target<-normalizePath(target,mustWork=TRUE)
work<-target
legacy<-root
plot_support<-file.path(script_dir,'plot_support')
dir.create(file.path(target,'figure'),recursive=TRUE,showWarnings=FALSE)
source(file.path(plot_support,'figure_theme.R'))
theme_set(theme_nature_control())
# Core conclusion: guide agreement alone does not distinguish target-specific joint responses.
# Evidence chain: raw agreement -> control construct null -> adjusted target evidence.
# Archetype: quantitative grid. R only; 183 mm; embedded editable vector text.
# Figure 8 describes coverage rather than inferring biological absence from missing assays.
save_revision<-function(p,name,w=183,h=130){
  stem<-file.path(target,'figure',name)
  if(Sys.info()[['sysname']]=='Darwin')grDevices::quartz(type='pdf',file=paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica')else grDevices::cairo_pdf(file=paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='Helvetica');print(p);dev.off()
  svglite::svglite(file.path(work,paste0(name,'.svg')),width=w/25.4,height=h/25.4);print(p);dev.off()
}
method_colors<-c(training_mean='#A0A0A0',ridge_static='#7A7A7A',elastic_net='#4C78A8',local_response_ridge='#E3A857',kernel_transition='#2A9D8F',conditional_rate_network='#7B6BA8',kernel_regression='#2A9D8F',additive_null='#A0A0A0',nonlinear_xgboost='#C95A4B')
# Figure 3: existing points and assignments, terminology update only.
s<-readLines(file.path(plot_support,'11_plot_main_figures.R'))
i<-grep('^# Figure 3 ',s)[1];j<-grep('^# Figure 4 ',s)[1]
block<-s[(i+1):(j-1)]
block<-gsub('"results/',paste0('"',root,'/results/'),block,fixed=TRUE)

for(expr in parse(text=paste(block,collapse='\n'))){if(is.call(expr)&&as.character(expr[[1]])%in%c('fwrite','save_pub_r'))next;eval(expr)}
p3d<-p3d+labs(x='Positive non-additivity AUROC')
p3b<-p3b+labs(y='Program-endpoint RMSE')
p3c<-p3c+labs(x='Binned program-tail error')
fig3<-(p3a|p3b)/(p3c|p3d)+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(size=8,face='bold'))
save_revision(fig3,'Fig3_model_benchmark',h=183)
# Same correction in extended-baseline supplement; no metric changed.
d<-fread(file.path(root,'results/figure_source_data/SuppFig1_extended_baselines.tsv'))
d[,method_label:=factor(method_label,levels=rev(unique(method_label)))]
colours<-c('Leave gene'='#4C78A8','Leave PFAM family'='#C95A4B')
mp<-function(metric,label,reference)ggplot(d,aes(.data[[metric]],method_label,colour=split_label,shape=split_label))+geom_vline(xintercept=reference,linetype=2,colour='#BDBDBD')+geom_line(aes(group=method_label),colour='#D9D9D9',linewidth=.4)+geom_point(size=1.8)+scale_colour_manual(values=colours)+labs(x=label,y=NULL,colour=NULL,shape=NULL)
p<-mp('synergy_auc','Positive non-additivity AUROC',.5)+mp('interaction_spearman','Spearman correlation',0)+plot_layout(guides='collect')+plot_annotation(tag_levels='a')&theme(legend.position='bottom')
save_revision(p,'SuppFig1_extended_baselines',h=90)
# Recompute all construct-level negative controls from saved response vectors.
x<-readRDS(file.path(legacy,'inputs/multiome_count_inputs.rds'));g<-as.data.table(x$guide);tg<-as.data.table(x$targets)
R<-x$original_rna_response;A<-x$original_atac_response
if(is.null(R)){base<-readRDS(file.path(legacy,'inputs/multiome.rds'));R<-base$rna_response;A<-base$atac_response}
cosine<-function(v,w){den<-sqrt(sum(v*v)*sum(w*w));if(!is.finite(den)||den<1e-12)return(NA_real_);sum(v*w)/den}
logo<-function(M)if(ncol(M)<2)NA_real_ else median(vapply(seq_len(ncol(M)),function(j)cosine(M[,j],rowMeans(M[,-j,drop=FALSE])),numeric(1)),na.rm=TRUE)
primary<-rbindlist(lapply(which(tg$context=='ATRT_human'),function(i){t<-tg[i];jj<-which(g$context=='ATRT_human'&g$guide_gene==t$guide_gene);data.table(gene=t$guide_gene,n_guides=length(jj),count_quality_gate=t$count_quality_gate,rna_cosine=logo(R[,jj,drop=FALSE]),atac_cosine=logo(A[,jj,drop=FALSE]))}))
primary[,agreement_pass:=count_quality_gate&rna_cosine>.1&atac_cosine>.1]
nt<-which(g$context=='ATRT_human'&g$is_control&g$n_cells>=15)
ntnull<-rbindlist(lapply(c(2L,3L),function(k)rbindlist(lapply(combn(nt,k,simplify=FALSE),function(jj){ref<-setdiff(nt,jj);rr<-as.numeric(R[,ref,drop=FALSE]%*%g$n_cells[ref]/sum(g$n_cells[ref]));ar<-as.numeric(A[,ref,drop=FALSE]%*%g$n_cells[ref]/sum(g$n_cells[ref]));rc<-logo(sweep(R[,jj,drop=FALSE],1,rr,'-'));ac<-logo(sweep(A[,jj,drop=FALSE],1,ar,'-'));data.table(n_guides=k,rna_cosine=rc,atac_cosine=ac,joint=min(rc,ac))}))))
ntnull[,agreement_pass:=rna_cosine>.1&atac_cosine>.1]
ntsum<-ntnull[,.(groups=.N,passes=sum(agreement_pass),fraction=mean(agreement_pass),joint95=quantile(joint,.95)),by=n_guides]
evaltarget<-primary[count_quality_gate==TRUE]
evaltarget[,nt_construct_p:=vapply(seq_len(.N),function(i){z<-ntnull[n_guides==evaltarget$n_guides[i],joint];(1+sum(z>=min(evaltarget$rna_cosine[i],evaltarget$atac_cosine[i])))/(1+length(z))},numeric(1))]
evaltarget[,nt_construct_fdr:=p.adjust(nt_construct_p,'BH')]
stopifnot(nrow(evaltarget)==27,sum(evaltarget$agreement_pass)==10,sum(evaltarget$nt_construct_fdr<.05)==0)
fwrite(ntnull,file.path(work,'non_targeting_constructs.tsv'),sep='\t');fwrite(ntsum,file.path(work,'non_targeting_summary.tsv'),sep='\t');fwrite(evaltarget,file.path(work,'target_specificity.tsv'),sep='\t')
allpoints<-rbindlist(list(ntnull[,.(rna_cosine,atac_cosine,group='Non-targeting groups')],evaltarget[,.(rna_cosine,atac_cosine,group='Supported targets')]))
p5a<-ggplot(allpoints,aes(rna_cosine,atac_cosine,colour=group))+geom_hline(yintercept=.1,linetype=2,linewidth=.35,colour='#777777')+geom_vline(xintercept=.1,linetype=2,linewidth=.35,colour='#777777')+geom_point(size=1.6,alpha=.75)+scale_colour_manual(values=c('Non-targeting groups'='#B7B7B7','Supported targets'='#4C78A8'))+coord_equal()+labs(x='RNA guide-agreement cosine',y='ATAC guide-agreement cosine',colour=NULL)+theme(legend.position='bottom',legend.text=element_text(size=6))
ntsum[,guides:=factor(paste(n_guides,'guides'),levels=c('2 guides','3 guides'))]
p5b<-ggplot(ntsum,aes(guides,fraction))+geom_col(width=.55,fill='#A0A0A0')+geom_text(aes(label=sprintf('%d/%d\n(%.1f%%)',passes,groups,100*fraction)),vjust=-.5,size=2.5)+scale_y_continuous(limits=c(0,.85),labels=scales::percent_format(accuracy=1),expand=expansion(mult=c(0,.02)))+labs(x='Non-targeting constructs',y='Fraction passing both 0.10 cutoffs')
evaltarget[,gene:=factor(gene,levels=gene[order(nt_construct_fdr,rna_cosine)])]
p5c<-ggplot(evaltarget,aes(nt_construct_fdr,gene,colour=agreement_pass))+geom_vline(xintercept=.05,linetype=2,colour='#777777',linewidth=.35)+geom_point(size=1.8)+scale_colour_manual(values=c('FALSE'='#B7B7B7','TRUE'='#4C78A8'),labels=c('Other supported targets','Passes agreement cutoffs'))+scale_x_continuous(limits=c(0,1),breaks=c(.05,.5,1))+labs(x='Non-targeting-construct BH q',y=NULL,colour=NULL)+theme(legend.position='bottom',legend.text=element_text(size=5.8),axis.text.y=element_text(size=6.2))
fig5<-p5a+p5b+p5c+plot_layout(design='AC\nBC',widths=c(1.2,1),heights=c(1.15,1))+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(size=8,face='bold'))
save_revision(fig5,'Fig5_multiome_control_roles',h=150)
# Figure 6: same graph statistics with component terminology.
i<-grep('^# Figure 6 ',s)[1];j<-grep('^# Figure 7 ',s)[1];block<-s[(i+1):(j-1)]
block<-gsub('"results/',paste0('"',root,'/results/'),block,fixed=TRUE)

block<-gsub('Effective number of response paths','Effective graph components',block,fixed=TRUE)
block<-gsub('Observed response-path diversity','Effective graph components',block,fixed=TRUE)
block<-gsub('Effective paths','Effective\ncomponents',block,fixed=TRUE)
for(expr in parse(text=paste(block,collapse='\n'))){if(is.call(expr)&&as.character(expr[[1]])%in%c('fwrite','save_pub_r'))next;eval(expr)}
save_revision(fig6,'Fig6_control_degeneracy',h=145)
# Figure 7: existing matched population data, explicitly nonexclusive tail fractions.
old<-readLines(file.path(plot_support,'plot_panels.R'));i<-grep('^# Figure 7 ',old)[1];j<-grep('^# Figure 8 ',old)[1];block<-old[(i+1):(j-1)]
block<-gsub('"results/',paste0('"',root,'/results/'),block,fixed=TRUE);
block<-gsub('Maximum program probability','Maximum program-tail fraction',block,fixed=TRUE)
block<-gsub('maximum program probability','maximum program-tail fraction',block,fixed=TRUE)
block<-gsub('colour="red"','colour="black"',block,fixed=TRUE)
for(expr in parse(text=paste(block,collapse='\n'))){if(is.call(expr)&&as.character(expr[[1]])%in%c('fwrite','save_pub_r'))next;eval(expr)}
save_revision(fig7,'Fig7_observed_support_sequences',h=130)
# Figure 8: specificity-aware primary accounting, no unsupported candidate classification.
ref<-readRDS(file.path(legacy,'reference/published_numerical_reference.rds'))
e<-copy(as.data.table(ref$evidence));e[,historical_class:=evidence_class]
e[,multiperturb_eligible:=multiperturb_eligible&gene%in%evaltarget[nt_construct_fdr<.05,gene]]
e[,datasets_with_causal_support:=as.integer(replogle_eligible)+as.integer(norman_single_eligible)+as.integer(multiperturb_eligible)]
e[,independent_contexts_with_causal_support:=as.integer(replogle_eligible|norman_single_eligible)+as.integer(multiperturb_eligible)]
e[,causal_regulator:=datasets_with_causal_support>=1]
e[,validated_control_node_candidate:=causal_regulator&independent_contexts_with_causal_support>=2&combination_consistent&path_specific]
cut<-quantile(e$observed_combinations[e$observed_combinations>0],.9)
e[,evidence_class:=fifelse(validated_control_node_candidate,'validated_control_node_candidate',fifelse(causal_regulator,'causal_regulator',fifelse(observed_combinations>0&observed_combinations>=cut,'association_hub',fifelse(genomewide_screened&is.finite(genomewide_effect_fdr)&genomewide_effect_fdr<.05,'state_marker','insufficient_causal_evidence'))))]
ec<-e[,.N,by=evidence_class];stopifnot(ec[evidence_class=='causal_regulator',N]==389,ec[evidence_class=='state_marker',N]==4128,ec[evidence_class=='insufficient_causal_evidence',N]==5377)
tasks<-copy(as.data.table(ref$tasks));tasks[dataset_id=='yan_2024_multiperturbseq'&status=='eligible',`:=`(status='excluded',exclusion_reason='non_targeting_construct_BH_q_not_passed')]
ts<-tasks[,.(eligible_tasks=sum(status=='eligible')),by=dataset_id]
ts[,dataset_label:=c(mixseq_2020='MIX-Seq',norman_2019_combinatorial='Norman',replogle_2022_k562_essential='Replogle',yan_2024_multiperturbseq='MultiPerturb')[dataset_id]]
ts<-ts[!is.na(dataset_label)];stopifnot(sum(ts$eligible_tasks)==705)
p8a<-ggplot(ts,aes(eligible_tasks,reorder(dataset_label,eligible_tasks),fill=dataset_label))+geom_col(width=.65)+geom_text(aes(label=eligible_tasks),hjust=-.25,size=2.5)+scale_x_continuous(expand=expansion(mult=c(0,.13)))+scale_fill_manual(values=c(Replogle='#4C78A8',Norman='#7B6BA8','MIX-Seq'='#E3A857',MultiPerturb='#2A9D8F'))+labs(x='Eligible endpoint tasks',y=NULL)+theme(legend.position='none')
p8b<-ggplot(ec,aes(N,reorder(evidence_class,N),fill=evidence_class))+geom_col(width=.65)+geom_text(aes(label=format(N,big.mark=',')),hjust=-.15,size=2.5)+scale_x_continuous(expand=expansion(mult=c(0,.2)))+scale_fill_manual(values=c(causal_regulator='#4C78A8',state_marker='#E3A857',insufficient_causal_evidence='#D9D9D9'))+scale_y_discrete(labels=c(causal_regulator='Perturbation-\nsupported',state_marker='State marker',insufficient_causal_evidence='Insufficient\nperturbation evidence'))+labs(x='Gene symbols',y=NULL)+theme(legend.position='none')
gc<-as.data.table(ref$gene_criterion_status)
coverage<-data.table(criterion=c('Both contexts assayed','Combination assayed','Eligible specificity contrast','All candidate evidence available'),available=c(sum((gc$replogle_assayed|gc$norman_assayed)&gc$multiome_assayed),sum(gc$observed_combinations>0),sum(gc$replogle_eligible|gc$norman_single_eligible),sum(gc$all_candidate_evidence_available)))
coverage[,missing:=nrow(gc)-available]
coverage[,criterion:=factor(criterion,levels=rev(criterion))]
cl<-melt(coverage,id.vars='criterion',variable.name='status',value.name='N')
p8c<-ggplot(cl,aes(N,criterion,fill=status))+geom_col(width=.6,position=position_stack(reverse=TRUE))+geom_text(data=coverage,aes(x=9894,y=criterion,label=paste0(available,' / 9,894 available')),inherit.aes=FALSE,hjust=1.02,size=2.4)+scale_fill_manual(values=c(available='#6C9FBF',missing='#EEEEEE'),labels=c(available='Available',missing='Unavailable'))+scale_x_continuous(limits=c(0,9894),breaks=c(0,2500,5000,7500,9894),labels=scales::comma)+labs(x='Assay coverage for candidate assessment (gene symbols)',y=NULL,fill=NULL)+theme(legend.position='bottom')
fig8<-p8a+p8b+p8c+plot_layout(design='AB\nCC',heights=c(1.1,1))+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(size=8,face='bold'))
save_revision(fig8,'Fig8_atlas_boundaries',h=130)
fwrite(e,file.path(work,'specificity_aware_gene_evidence.tsv.gz'),sep='\t')
fwrite(tasks,file.path(work,'specificity_aware_task_catalogue.tsv.gz'),sep='\t')
fwrite(ec,file.path(work,'specificity_aware_evidence_counts.tsv'),sep='\t')
fwrite(ts,file.path(work,'specificity_aware_task_counts.tsv'),sep='\t')
fwrite(coverage,file.path(work,'candidate_assay_coverage.tsv'),sep='\t')
print(ntsum);print(evaltarget[,.(targets=.N,agreement_pass=sum(agreement_pass),specificity_pass=sum(nt_construct_fdr<.05),minimum_q=min(nt_construct_fdr))]);print(ts);print(ec);print(coverage)
writeLines(capture.output(sessionInfo()),file.path(work,'R_session_info.txt'))
# Figure 1: inventory and inferential denominators, without workload thresholds.
cells<-fread(file.path(root,'results/figure_source_data/fig_phase0_cells.tsv'))
cells[,dataset:=factor(dataset,levels=rev(dataset))]
p1a<-ggplot(cells,aes(n_cells/1000,dataset,fill=role_group))+geom_col(width=.66)+geom_text(aes(label=label),hjust=-.15,size=2.15,family='Helvetica')+scale_fill_manual(values=palette_control[c('genetic','combination','multiome','dynamics')])+scale_x_continuous(expand=expansion(mult=c(0,.14)),labels=scales::label_number(big.mark=','))+labs(x='Inventoried cells (thousands)',y=NULL,fill='Evidence role')+theme(legend.position='bottom',legend.key.width=grid::unit(3.5,'mm'),axis.line.y=element_blank(),axis.ticks.y=element_blank())
inventory<-data.table(label=c('Inventoried single cells','Breadth-screen genes','Eligible endpoint tasks','LARRY study / arms','Joint RNA + ATAC resource','MIX-Seq tasks outside\ncell inventory'),value=c(format(sum(cells$n_cells),big.mark=','),'9,867',as.character(sum(ts$eligible_tasks)),'1 / 2','1','151'))
inventory[,label:=factor(label,levels=rev(label))]
p1b<-ggplot(inventory,aes(y=label))+geom_text(aes(x=.12,label=value),hjust=0,size=2.5,family='Helvetica')+scale_x_continuous(limits=c(0,1),expand=c(0,0))+labs(x=NULL,y=NULL)+theme_void(base_size=6.5,base_family='Helvetica')+theme(axis.text.y=element_text(size=6.1,colour='black',hjust=1),plot.margin=margin(5.5,5.5,5.5,8))
fig1<-p1a+p1b+plot_layout(widths=c(1.45,1))+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(size=8,face='bold'))
save_revision(fig1,'Fig1_data_backbone',h=95)
fwrite(inventory,file.path(work,'source_inventory_denominators.tsv'),sep='\t')
fwrite(cells,file.path(work,'source_inventory_cells.tsv'),sep='\t')
# Record the deterministic specificity diagnostic explicitly.
fwrite(evaltarget[,.(supported_targets=.N,agreement_passes=sum(agreement_pass),specificity_passes=sum(nt_construct_fdr<.05),minimum_BH_q=min(nt_construct_fdr))],file.path(work,'target_specificity_summary.tsv'),sep='\t')
