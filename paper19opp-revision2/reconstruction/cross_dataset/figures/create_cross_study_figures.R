#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);project<-normalizePath(args[1]);base<-file.path(project,'paper19opp-revision2')
.libPaths(c(file.path(project,'.Rlib'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(svglite);library(ragg)})
data_dir<-file.path(base,'reconstruction/cross_dataset');archive<-file.path(data_dir,'figures');src<-file.path(data_dir,'figure_source_data');dir.create(src,showWarnings=FALSE,recursive=TRUE)
palette<-c(Replogle='#4C78A8',Norman='#7B6BA8','MIX-Seq'='#B8873A')
grey<-'#737373';light<-'#EDF1F5';ink<-'#282828'
font<-'Helvetica'
writeLines(c('Figure 1: A common measurement definition supports a source-held-out response benchmark with independent context-state and response samples.',
'Archetype: schematic-led composite; source units -> common genes/independent splits -> held-out study.',
'Figure 2: Adding contemporaneous control-state features does not improve the prespecified cross-study pilot-response reconstruction.',
'Archetype: asymmetric quantitative composite; hero paired error difference -> absolute matched baselines -> computational split sensitivity.',
'Both use183mm width,7pt body and 8pt panel labels, editable PDF/SVG; no images or biological causal arrows.'),file.path(archive,'figure_contract.txt'))
save_draw<-function(draw,stem,width_mm,height_mm){
 w<-width_mm/25.4;h<-height_mm/25.4
 svglite::svglite(file.path(archive,paste0(stem,'.svg')),width=w,height=h);draw();dev.off()
 if(Sys.info()[['sysname']]=='Darwin')grDevices::quartz(type='pdf',file=file.path(base,'figure',paste0(stem,'.pdf')),width=w,height=h,family=font)else grDevices::cairo_pdf(file.path(base,'figure',paste0(stem,'.pdf')),width=w,height=h,family=font)
 draw();dev.off()
 ragg::agg_png(file.path(archive,paste0(stem,'_preview.png')),width=w,height=h,units='in',res=300);draw();dev.off()
 ragg::agg_png(file.path(archive,paste0(stem,'_actual_size.png')),width=w,height=h,units='in',res=150);draw();dev.off()
}
# Design schematic: coordinates in final-publication millimetres.
text_mm<-function(label,x,y,size=7,bold=FALSE,just='centre',colour=ink){grid.text(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,gp=gpar(fontfamily=font,fontsize=size,fontface=if(bold)'bold'else'plain',col=colour,lineheight=1.16))}
box_mm<-function(x,y,w,h,fill='white',border='#C4CAD0'){grid.rect(x=unit(x,'mm'),y=unit(y,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),gp=gpar(fill=fill,col=border,lwd=.55))}
arrow_mm<-function(x0,y0,x1,y1,colour='#7B838B'){grid.lines(x=unit(c(x0,x1),'mm'),y=unit(c(y0,y1),'mm'),arrow=arrow(type='closed',length=unit(1.1,'mm')),gp=gpar(col=colour,fill=colour,lwd=.6))}
design<-function(){
 grid.newpage();pushViewport(viewport(width=unit(183,'mm'),height=unit(156,'mm')))
 text_mm('a',3,151,8,TRUE,just='left');text_mm('Source-specific response units',9,151,7,TRUE,just='left')
 cards<-data.table(x=c(31,91.5,152),study=c('Replogle','Norman','MIX-Seq'),assay=c('K562 · CRISPRi','K562 · CRISPRa','Chemical perturbations'),cells=c('276,108 QC-retained cells','104,507 QC-retained cells','120,003 chemical cells'),units=c('2,151 gene × GEM units','2,180 condition × GEM units','527 drug × context units'),support=c('358 genes · 48 technical GEMs','280 conditions · 8 technical GEMs','13 drugs · 109 cell lines'),extra=c('one published experiment','one published experiment','169 matched assay contexts'))
 for(i in seq_len(nrow(cards))){
  cc<-cards[i];box_mm(cc$x,126,55,36,fill='#F7F8FA');grid.lines(x=unit(c(cc$x-27.5,cc$x+27.5),'mm'),y=unit(c(144,144),'mm'),gp=gpar(col=palette[cc$study],lwd=1.7))
  text_mm(cc$study,cc$x,140,7,TRUE);text_mm(cc$assay,cc$x,135);text_mm(cc$cells,cc$x,129.5)
  text_mm(cc$units,cc$x,123.5,7,TRUE);text_mm(cc$support,cc$x,117.5);text_mm(cc$extra,cc$x,112)
 }
 text_mm('Three held-out studies; GEM partitions describe within-experiment sampling.',91.5,103.5)
 text_mm('b',3,96,8,TRUE,just='left');text_mm('Shared measurements and independent cell samples',9,96,7,TRUE,just='left')
 box_mm(40.5,70,73,42,fill='#F7F8FA')
 text_mm('Metadata-defined gene intersection',40.5,87,7,TRUE)
 text_mm('Replogle · Norman · MIX-Seq',40.5,81.8);text_mm('Adamson · Frangieh',40.5,77.2)
 text_mm('7,848 measured genes',40.5,70.5,7,TRUE)
 text_mm('Duplicate counts summed before log1p',40.5,64.8);text_mm('CP10k denominator: the same 7,848 genes',40.5,60)
 text_mm('3 signed modules + 9 Hallmark summaries',40.5,54)
 text_mm('Contemporaneous matched controls',132,87,7,TRUE)
 for(x in c(100,129,158)){arrow_mm(129,83,x,78)}
 for(i in 1:3){box_mm(c(100,129,158)[i],73,25,9,fill=if(i==3)'#E7EFF7'else'#F1F1F1');text_mm(paste('control',LETTERS[i]),c(100,129,158)[i],73)}
 box_mm(100,57.5,25,15);box_mm(129,57.5,25,15);text_mm('treated A\n− control A',100,57.5);text_mm('treated B\n− control B',129,57.5)
 arrow_mm(100,50,100,47);arrow_mm(129,50,129,47)
 text_mm('response A',100,43.5,7,TRUE);text_mm('response B',129,43.5,7,TRUE)
 arrow_mm(100,68.5,100,65.5);arrow_mm(129,68.5,129,65.5)
 arrow_mm(158,68.5,158,62.5);box_mm(158,55,25,14,fill='#E7EFF7');text_mm('control-state\nfeatures C',158,55)
 text_mm('Disjoint cell subsets',131,36.5)
 text_mm('c',3,34,8,TRUE,just='left');text_mm('Study-held-out prediction',9,34,7,TRUE,just='left')
 rows<-data.table(y=c(26,18,10),train=c('Norman + MIX-Seq','Replogle + MIX-Seq','Replogle + Norman'),test=c('Replogle','Norman','MIX-Seq'))
 for(i in 1:3){box_mm(42,rows$y[i],66,6.5,fill='#F5F5F5');text_mm(paste('Train:',rows$train[i]),42,rows$y[i]);arrow_mm(76,rows$y[i],91,rows$y[i]);box_mm(119,rows$y[i],53,6.5,fill='white',border=palette[rows$test[i]]);text_mm(paste('Test:',rows$test[i]),119,rows$y[i])}
 text_mm('Scaling, ridge penalties and model fitting use training studies only.',91.5,3)
 popViewport()
}
save_draw(design,'Fig_cross_study_design',183,156)
cohorts<-fread(file.path(data_dir,'analysis_denominators.tsv'));fwrite(cohorts,file.path(src,'Fig_cross_study_design_cohorts.tsv'),sep='\t')
fwrite(fread(file.path(data_dir,'final_module_gene_coverage.tsv'))[,.(measured=sum(in_common),specified=.N),by=.(program,direction)],file.path(src,'Fig_cross_study_design_coverage.tsv'),sep='\t')
theme_fig<-function(){theme_classic(base_size=7,base_family=font)+theme(axis.line=element_line(linewidth=.3,colour=ink),axis.ticks=element_line(linewidth=.3,colour=ink),axis.ticks.length=unit(1.1,'mm'),axis.text=element_text(size=7,colour=ink),axis.title=element_text(size=7,colour=ink),strip.background=element_blank(),strip.text=element_text(size=7,face='bold'),panel.grid=element_blank(),legend.position='none',plot.margin=margin(4,6,4,4),plot.tag=element_text(size=8,face='bold'),plot.tag.position=c(0,1))}
diff<-fread(file.path(data_dir,'direct_state_value_cluster_comparisons.tsv'))[task=='response'&axis_group=='core3'&bootstrap=='crossed_intervention_context'&intervention_definition=='complete_condition']
diff[,dataset:=factor(dataset,levels=c('MIX-Seq','Norman','Replogle'))]
diff[,source_label:=c('MIX-Seq  (527 units)','Norman  (2,180 units)','Replogle  (2,151 units)')[match(as.character(dataset),c('MIX-Seq','Norman','Replogle'))]]
fwrite(diff,file.path(src,'Fig_cross_study_transfer_a.tsv'),sep='\t')
hero<-ggplot(diff,aes(x=rmse_difference,y=dataset,colour=dataset))+geom_vline(xintercept=0,colour='#929292',linewidth=.35,linetype=2)+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',width=.13,linewidth=.55)+geom_point(size=2,shape=16)+geom_text(aes(x=ci_high+.0012,label=sprintf('%+.4f',rmse_difference)),colour=ink,size=7/.pt,hjust=0)+scale_colour_manual(values=palette)+scale_y_discrete(labels=setNames(diff$source_label,as.character(diff$dataset)),expand=expansion(add=.55))+scale_x_continuous(limits=c(-.0045,.0495),breaks=c(0,.01,.02,.03,.04),labels=function(x)sprintf('%.2f',x),expand=expansion(mult=0))+labs(x='Core-program RMSE difference (pilot + control state − pilot ridge)',y=NULL)+theme_fig()+theme(plot.margin=margin(5,8,5,4))
absolute<-fread(file.path(data_dir,'primary_response_rmse_by_axis_group.tsv'))[axis_group=='core3'&model%in%c('training_mean','pilot_identity','pilot_response','pilot_plus_state')]
method_labels<-c(training_mean='Training mean',pilot_identity='Pilot A unchanged',pilot_response='Pilot ridge',pilot_plus_state='Pilot + control state')
absolute[,method:=factor(method_labels[model],levels=rev(unname(method_labels)))];absolute[,dataset:=factor(dataset,levels=c('Replogle','Norman','MIX-Seq'))]
fwrite(absolute,file.path(src,'Fig_cross_study_transfer_b.tsv'),sep='\t')
comparators<-ggplot(absolute,aes(x=rmse,y=method,colour=dataset))+geom_segment(aes(x=0,xend=rmse,yend=method),colour='#D9DDE1',linewidth=.4)+geom_point(aes(shape=model),size=1.8)+scale_shape_manual(values=c(training_mean=1,pilot_identity=2,pilot_response=16,pilot_plus_state=18))+scale_colour_manual(values=palette)+facet_grid(dataset~.,switch='y')+scale_x_continuous(limits=c(0,.27),breaks=c(0,.1,.2),expand=expansion(mult=c(0,.03)))+labs(x='Held-out core-program RMSE',y=NULL)+theme_fig()+theme(strip.placement='outside',strip.text.y.left=element_text(angle=0,size=7),panel.spacing.y=unit(5,'mm'),plot.margin=margin(6,4,4,3))
am<-fread(file.path(data_dir,'response_metrics_by_axis.tsv'))[axis%in%c('apoptosis','proliferation','stress')&model%in%c('pilot_response','pilot_plus_state')]
cm<-am[,.(rmse=sqrt(mean(rmse^2))),by=.(partition_id,dataset,model)]
splits<-dcast(cm,partition_id+dataset~model,value.var='rmse');splits[,difference:=pilot_plus_state-pilot_response];splits[,dataset:=factor(dataset,levels=c('MIX-Seq','Norman','Replogle'))]
splits[,y:=as.numeric(dataset)+(partition_id-5.5)*.033]
summ<-splits[,.(median=median(difference),y=as.numeric(dataset[1])),by=dataset]
fwrite(splits,file.path(src,'Fig_cross_study_transfer_c.tsv'),sep='\t')
sensitivity<-ggplot(splits,aes(x=difference,y=y,colour=dataset))+geom_vline(xintercept=0,linewidth=.35,linetype=2,colour='#929292')+geom_point(size=1.4,alpha=.75)+geom_segment(data=summ,aes(x=median,xend=median,y=y-.23,yend=y+.23),linewidth=.75,colour=ink,inherit.aes=FALSE)+scale_colour_manual(values=palette)+scale_x_continuous(limits=c(-.0045,.0425),breaks=c(0,.02,.04),labels=function(x)sprintf('%.2f',x),expand=expansion(mult=0))+scale_y_continuous(limits=c(.45,3.55),breaks=1:3,labels=c('MIX-Seq','Norman','Replogle'),expand=expansion(mult=0))+labs(x='Core-program RMSE difference',y=NULL)+theme_fig()+theme(plot.margin=margin(6,6,4,5))
bottom<-comparators+sensitivity+plot_layout(widths=c(1.6,1))
transfer<-(hero/bottom)+plot_layout(heights=c(1,1.7))+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(size=8,face='bold',family=font))
save_draw(function()print(transfer),'Fig_cross_study_transfer',183,150)
writeLines(c('Fig. 1 | Common measurements and source-held-out response evaluation.',
'a, Discovery cohorts and metadata-eligible analysis units. Replogle and Norman each represent one published K562 experiment; their 48 and 8 GEM partitions are technical contexts, not independent biological replicates. MIX-Seq contributes 527 chemical drug/context tasks from 13 drugs and 109 cell lines. Units require at least 10 treated and 60 matched-control cells, without response-magnitude filtering. Cell counts refer to reconstructed QC-retained source cells; the MIX-Seq count includes chemical-screen cells only, with ancillary genetic cells excluded. Not every reconstructed chemical context meets the benchmark cell-support requirements.',
'b, The measured gene intersection was frozen using metadata from the three discovery studies and Adamson and Frangieh before external response evaluation. Duplicated gene-symbol counts were summed before log1p CP10k normalization with the same 7,848-gene denominator. Three signed modules and nine Hallmark mean-expression scores were defined from their measured members. Controls were randomly divided into disjoint A/B reference subsets and C context-feature cells; treated A/B cells were also disjoint. Response A and response B are treated-minus-matched-control mean-score contrasts. C describes contemporaneous matched-control state, not paired individual pretreatment state.',
'c, Each published discovery study is held out in turn. Feature centring/scaling, study-balanced ridge fitting and penalty selection use only the other two studies, with inner source-held-out tuning. A pilot-response reconstruction uses observed A responses, while the strict control-state attempt does not. The independent-cell split assesses sampling precision within the observed assay, separately from source-level transfer.'),file.path(archive,'Fig_cross_study_design_legend.txt'))
writeLines(c('Fig. 2 | Contemporaneous control-state features do not improve the prespecified pilot-response reconstruction.',
'a, Added-state minus pilot-ridge RMSE for the three apoptosis-, proliferation- and stress-associated modules on identical held-out units. Negative values favour adding control state. Points are primary-split differences; horizontal bars are 95% percentile intervals from 2,000 crossed intervention/context bootstrap draws conditional on fixed predictions and the cell split. Interventions and technical GEMs are resampled in Replogle (358 genes, 48 GEMs) and Norman (280 conditions, 8 GEMs); MIX-Seq resamples 13 drugs and 109 complete cell lines. All feature scaling, model fitting and ridge-penalty selection exclude the held-out study. These intervals do not quantify independent biological replication or refit the models. A Norman anchor-gene sensitivity is supplied separately.',
'b, Absolute RMSE of the training-study mean, unchanged pilot response A, pilot ridge and pilot plus control-state models. Each observed unit and each core module have equal evaluation weight within its held-out study. Training mean is estimated in the two training studies. Pilot A unchanged directly predicts B from its independent A response without cross-study fitting. Control features include twelve C-state means and variances, with sampling-support terms also present in the pilot-ridge comparator. The chemical source shows particularly poor transfer of cross-study shrinkage relative to its observed pilot.',
'c, Added-state minus pilot-ridge core-program RMSE under ten fixed randomized cell partitions. Dots are computational split results, and black ticks mark their medians; they are not biological replicates or confidence intervals. Partition 1 was prespecified as primary, and no split was selected by its performance. Replogle partition effects include both signs. The results concern an observed pilot and contemporaneous control state; they do not establish that biological pretreatment state is irrelevant.'),file.path(archive,'Fig_cross_study_transfer_legend.txt'))
writeLines(c('Backend: R for all drawing, SVG/PDF export and PNG preview.',
'Dimensions: design183x156mm; transfer183x150mm.',
'Fonts: Helvetica7pt body,8pt bold lowercase panel labels; standard vector text.',
'No raster source images, image manipulation or rescaled biological imagery.',
'Figure1 is a design schematic; arrows denote data-processing relationships only.',
'Figure2 source data are exported TSVs; primary core3 model differences use identical units and fixed crossed-cluster intervals.',
'Quantitative QA: all plotted values read directly from completed analysis outputs; three studies, not unit counts, define source holdouts.',
'Visual QA pending direct inspection of R-rendered full-size previews.'),file.path(archive,'figure_quality_scope.txt'))
cat('Created two editable183mm vector figures and R previews.\n')
