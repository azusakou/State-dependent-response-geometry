#!/usr/bin/env Rscript
# Fixed generic essentiality comparator from an independently released DepMap matrix.
suppressPackageStartupMessages({library(data.table);library(digest)})
a<-commandArgs(TRUE);if(length(a)<1)stop('Usage: Rscript extract_generic_essentiality_baseline.R CRISPRGeneEffect.csv [output_directory]')
p<-normalizePath(a[1],mustWork=TRUE);out<-if(length(a)>1)a[2]else'.';dir.create(out,recursive=TRUE,showWarnings=FALSE)
r<-fread(file.path(out,'harmonized_discovery_ranking.tsv'));g<-r$gene
mapping<-file.path(out,'discovery_gene_identifier_mapping.tsv')
if(file.exists(mapping)){
 map<-fread(mapping,colClasses=c(gene='character',ENTREZID='character'));map<-map[gene%in%g]
}else{
 suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db)})
 map<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=g,keytype='SYMBOL',columns='ENTREZID'));setnames(map,'SYMBOL','gene');missing<-setdiff(g,map[!is.na(ENTREZID),gene]);if(length(missing)){q<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=missing,keytype='ALIAS',columns='ENTREZID'));setnames(q,'ALIAS','gene');map<-rbind(map,q,fill=TRUE)}
 map<-unique(map[!is.na(ENTREZID)]);single<-map[,.(ids=uniqueN(ENTREZID)),by=gene][ids==1,gene];map<-map[gene%in%single];fwrite(map,mapping,sep='\t')
}
header<-names(fread(p,nrows=0));ids<-sub('^.*\\(([0-9]+)\\)$','\\1',header);take<-which(ids%in%map$ENTREZID);d<-fread(p,select=take)
b<-data.table(ENTREZID=ids[take],depmap_generic_essentiality=-vapply(d,mean,numeric(1),na.rm=TRUE),depmap_lines=vapply(d,function(x)sum(is.finite(x)),integer(1)));b<-merge(map,b,by='ENTREZID');b<-unique(b,by='gene');fwrite(b,file.path(out,'discovery_generic_essentiality_baseline.tsv'),sep='\t')
fwrite(data.table(source=p,bytes=file.info(p)$size,sha256=digest(p,file=TRUE,algo='sha256'),release='DepMap24Q4',source_dataset_url='https://doi.org/10.25452/figshare.plus.27993248.v1',method_doi='10.1186/s13059-021-02540-7',statistic='Negative mean over all finite source cellline effects',role='Comparator only, not a response-validation endpoint'),file.path(out,'generic_essentiality_source.tsv'),sep='\t')
