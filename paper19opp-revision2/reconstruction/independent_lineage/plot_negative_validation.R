#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
o<-'paper19opp-revision2/reconstruction/independent_lineage';fig<-'paper19opp-revision2/figure/SuppFig_independent_lineage.pdf'
# Figure contract: earlier clone states add no demonstrated improvement beyond
# treatment/sampling information in an independent held-lineage benchmark.
# Schematic-led composite, R-only,183x102mm, editable embedded PDF text.
base<-theme_classic(base_size=7,base_family='Helvetica')+theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),axis.text=element_text(color='#2B2B2B'),plot.title=element_text(size=8,face='bold'),plot.margin=margin(5,5,4,5),legend.position='none')
boxes<-data.table(xmin=c(.17,.01,.55,.55),xmax=c(.84,.44,.99,.99),ymin=c(.84,.57,.57,.31),ymax=c(.98,.75,.75,.49),x=c(.505,.225,.77,.77),y=c(.91,.66,.66,.40),label=c('Barcoded clone after first treatment','Earlier RNA sample\nmeasured in relatives','Remaining relatives\nknown second treatment','Later RNA sample\nmeasured descendants'))
p1<-ggplot()+geom_rect(data=boxes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='#F4F5F6',color='#7A7A7A',linewidth=.35)+geom_text(data=boxes,aes(x=x,y=y,label=label),size=2.35,family='Helvetica',lineheight=1.05)+
 annotate('segment',x=.42,xend=.225,y=.84,yend=.76,linewidth=.35,arrow=grid::arrow(length=grid::unit(1.2,'mm'),type='closed'))+
 annotate('segment',x=.59,xend=.77,y=.84,yend=.76,linewidth=.35,arrow=grid::arrow(length=grid::unit(1.2,'mm'),type='closed'))+
 annotate('segment',x=.77,xend=.77,y=.57,yend=.50,linewidth=.35,arrow=grid::arrow(length=grid::unit(1.2,'mm'),type='closed'))+
 annotate('segment',x=.225,xend=.225,y=.57,yend=.40,linewidth=.35,linetype=2,color='#4C78A8')+
 annotate('segment',x=.225,xend=.54,y=.40,yend=.40,linewidth=.35,linetype=2,color='#4C78A8',arrow=grid::arrow(length=grid::unit(1.2,'mm'),type='closed'))+
 annotate('text',x=.29,y=.50,label='Prediction',size=2.2,hjust=0,color='#4C78A8',family='Helvetica')+
 annotate('text',x=.50,y=.265,label='Different cells; shared ancestral lineage',size=2.35,family='Helvetica')+
 annotate('text',x=.50,y=.165,label='Primary: \u22653 earlier / \u22652 later cells\n257 paired tasks; 166 lineages',size=2.45,color='#4C78A8',family='Helvetica',lineheight=1.08)+
 annotate('text',x=.50,y=.045,label='Strict: \u226510 earlier / \u22655 later cells\n85 paired tasks; 52 lineages',size=2.45,color='#7A7A7A',family='Helvetica',lineheight=1.08)+
 coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip='off')+labs(title='Measured clone states and later relatives')+theme_void(base_family='Helvetica')+theme(plot.title=element_text(size=8,face='bold'),plot.margin=margin(5,10,4,4))
d<-fread(file.path(o,'paired_model_comparisons.tsv'));keep<-c('state_ridge_minus_condition_mean','state_ridge_minus_counts_ridge','state_ridge_minus_parent_persistence');d<-d[comparison%in%keep]
d[,baseline:=factor(comparison,levels=rev(keep),labels=rev(c('Treatment mean','Treatment/counts','Parent persistence')))];d[,support:=factor(setting,levels=c('primary','strict'),labels=c('Primary','Strict'))]
p2<-ggplot(d,aes(x=difference,y=baseline,color=support))+geom_vline(xintercept=0,linetype=2,color='#BBBBBB',linewidth=.35)+geom_segment(aes(x=lower,xend=upper,yend=baseline),linewidth=.65)+geom_point(size=1.8)+facet_grid(support~.,scales='free_y',space='free_y')+
 scale_color_manual(values=c(Primary='#4C78A8',Strict='#7A7A7A'))+scale_x_continuous(breaks=c(-.4,-.2,0,.1),limits=c(-.43,.11))+
 labs(title='Added states versus simple baselines',x='State-model minus baseline RMSE\n(training-SD units)',y=NULL)+base+theme(strip.background=element_blank(),strip.text.y=element_text(size=7,face='bold'),axis.title.x=element_text(size=7),axis.text.y=element_text(size=6.8))
p<-(p1|p2)+plot_layout(widths=c(.94,1.06))+plot_annotation(tag_levels='a',theme=theme(plot.tag=element_text(size=9,face='bold',family='Helvetica')))
grDevices::quartz(type='pdf',file=fig,width=183/25.4,height=102/25.4,family='Helvetica');print(p);dev.off()
ragg::agg_png(file.path(o,'independent_lineage_figure_preview.png'),width=183,height=102,units='mm',res=200);print(p);dev.off()
svglite::svglite(file.path(o,'independent_lineage_figure_editable.svg'),width=183/25.4,height=102/25.4);print(p);dev.off()
fwrite(d,file.path(o,'independent_lineage_figure_source_data.tsv'),sep='\t')
writeLines('Measured earlier clone states did not improve held-lineage descendant-state prediction beyond treatment and sampling information. a, Experimental and predictive design in the independent barcoded melanoma study. Earlier and later RNA profiles come from different sampled relatives of an ancestral lineage; the dashed arrow denotes prediction. All histories and branches of each lineage remain in the same validation fold. Primary support requires at least three earlier and two later cells; strict support requires ten and five. b, State-model minus simple-baseline joint RMSE. Outcomes are the three apoptosis-, proliferation- and stress-associated expression scores, each scaled by its outer-training SD; negative differences favour the state model. Points are five-fold held-lineage estimates with each lineage receiving total weight one. Lines span 95% percentile intervals from 249 whole-lineage bootstrap draws that repeat feature/outcome scaling, model fitting, inner tuning and comparator selection within fixed grouped folds. The training-selected simple comparator was treatment/count ridge in all five original folds of both support settings. These are conditional computational uncertainty intervals within one source experiment, not between-culture biological replication. Descendant sampled-cell count is a post-sampling covariate. No additional predictive value was demonstrated, and the intervals do not establish biological equivalence.',file.path(o,'independent_lineage_figure_legend.txt'))
