#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1]else'paper19opp-revision2/reproducibility/lineage_simulation'
fig<-if(length(args)>1)args[2]else'paper19opp-revision2/figure/SuppFigS6_lineage_sampling_null.pdf'
# Figure contract: Under known stress-independent intrinsic probabilities,
# descendant sampling changes covariance rank; a deep-clone restriction changes
# coefficient variability. Quantitative grid, same source-specific units, 183x151 mm.
d<-fread(file.path(out,'simulation_replicates.tsv'))[latent_logit_sd==1]
keys<-c('oracle_probability','recorded_depth','equal_depth_10','permuted_depth','random_neighbors_recorded_depth','deep_oracle_probability','deep_recorded_depth')
labels<-c('Intrinsic probabilities','Recorded descendants','Ten descendants','Permuted depths','Random donor neighbours','Deep clones: intrinsic','Deep clones: recorded')
d<-d[scenario%in%keys];d[,scenario:=factor(scenario,levels=rev(keys),labels=rev(labels))]
d[,type:=ifelse(grepl('Intrinsic|intrinsic',as.character(scenario)),'Intrinsic','Sampled')]
s<-d[,.(rank_median=median(median_rank),rank_lower=quantile(median_rank,.025),rank_upper=quantile(median_rank,.975),slope_median=median(slope),slope_lower=quantile(slope,.025),slope_upper=quantile(slope,.975)),by=.(dataset,scenario,type)]
cols<-c(Intrinsic='#7A7A7A',Sampled='#4C78A8')
base<-theme_classic(base_size=7,base_family='Helvetica')+theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),axis.text=element_text(color='#2B2B2B'),axis.title=element_text(size=7),plot.title=element_text(size=8,face='bold'),plot.subtitle=element_text(size=6.7),plot.margin=margin(5,8,5,4),legend.position='none',panel.grid.major.x=element_line(color='#EEEEEE',linewidth=.25))
make_plot<-function(nm,metric,title,subtitle){
 a<-s[dataset==nm];a[,`:=`(value=get(paste0(metric,'_median')),lower=get(paste0(metric,'_lower')),upper=get(paste0(metric,'_upper')))]
 g<-ggplot(a,aes(y=scenario,x=value,color=type))+geom_segment(aes(x=lower,xend=upper,yend=scenario),linewidth=.6)+geom_point(size=1.6)+scale_color_manual(values=cols)+labs(x=if(metric=='rank')'Median local effective rank'else'Adjusted stress coefficient',y=NULL,title=title,subtitle=subtitle)+base
 if(metric=='slope')g<-g+geom_vline(xintercept=0,color='#BBBBBB',linewidth=.3,linetype=2)
 g
}
p1<-make_plot('larry_cytokine','rank','Cytokine resource','Full: 1,268 units; deep: 92 units')
p2<-make_plot('larry_in_vitro','rank','In vitro resource','Full: 1,243 units; deep: 284 units')
p3<-make_plot('larry_cytokine','slope','Cytokine resource','Known-null intrinsic probabilities')
p4<-make_plot('larry_in_vitro','slope','In vitro resource','Known-null intrinsic probabilities')
p<-(p1|p2)/(p3|p4)+plot_annotation(tag_levels='a',theme=theme(plot.tag=element_text(size=9,face='bold',family='Helvetica')))
# Export vector PDF and an R-rendered preview; no cross-backend rendering.
grDevices::quartz(type='pdf',file=fig,width=183/25.4,height=151/25.4,family='Helvetica');print(p);dev.off()
ragg::agg_png(file.path(out,'figure_preview.png'),width=183,height=151,units='mm',res=180);print(p);dev.off()
svglite::svglite(file.path(out,'figure_editable.svg'),width=183/25.4,height=151/25.4);print(p);dev.off()
fwrite(s,file.path(out,'figure_source_data.tsv'),sep='\t')
