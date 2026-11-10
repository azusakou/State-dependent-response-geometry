suppressPackageStartupMessages({library(data.table);library(Matrix);library(hdf5r);library(digest)})

verify_locks<-function(out){
 for(nm in c('membership_lock_manifest.tsv','scoring_lock_manifest.tsv')){
  x<-fread(file.path(out,nm));stopifnot(all(vapply(file.path(out,x$file),digest,character(1),file=TRUE,algo='sha256')==x$sha256))
 }
}
annotation<-function(f,group,name){
 o<-f[[group]][[name]]
 if(inherits(o,'H5Group')){code<-as.integer(o[['codes']][]);cats<-o[['categories']][];v<-rep(NA_character_,length(code));v[code>=0L]<-cats[code[code>=0L]+1L];return(v)}
 v<-o[];if('__categories'%in%names(f[[group]])&&name%in%names(f[[group]][['__categories']]))v<-f[[group]][['__categories']][[name]][][as.integer(v)+1L]
 v
}
make_coeff<-function(genes,spec){
 make<-function(direction){i<-j<-integer();v<-numeric();for(k in seq_along(spec)){gs<-spec[[k]][[direction]];if(!length(gs))next;stopifnot(all(gs%in%genes));i<-c(i,match(gs,genes));j<-c(j,rep(k,length(gs)));v<-c(v,rep(1/length(gs),length(gs)))};sparseMatrix(i=i,j=j,x=v,dims=c(length(genes),length(spec)),dimnames=list(genes,names(spec)))}
 list(positive=make('positive'),negative=make('negative'),np=lengths(lapply(spec,`[[`,'positive')),nn=lengths(lapply(spec,`[[`,'negative')))
}
score_logs<-function(lg,genes,spec,coeff,assigned_gene=NULL){
 p<-as.matrix(t(lg)%*%coeff$positive);n<-as.matrix(t(lg)%*%coeff$negative)
 if(!is.null(assigned_gene))for(g in unique(assigned_gene)){
  if(is.na(g)||!g%in%genes)next
  ii<-which(assigned_gene==g);self<-lg[match(g,genes),ii]
  pos<-which(vapply(spec,function(s)g%in%s$positive,logical(1)))
  neg<-which(vapply(spec,function(s)g%in%s$negative,logical(1)))
  if(length(pos))for(k in pos)p[ii,k]<-if(coeff$np[k]>1L)(coeff$np[k]*p[ii,k]-self)/(coeff$np[k]-1L)else NA_real_
  if(length(neg))for(k in neg)n[ii,k]<-if(coeff$nn[k]>1L)(coeff$nn[k]*n[ii,k]-self)/(coeff$nn[k]-1L)else 0
 }
 p-n
}

get_metadata<-function(study,project,out,prep){
 study_name<-study;roster<-fread(file.path(out,'fixed_study_gene_rosters.tsv'))[study==study_name,gene]
 if(study=='K562'){
  md<-copy(prep$source_cells);md[,context:='K562'];md[,gem:=as.character(gem_group)];md[,technical_group_type:='GEM'];md[,genuine_GEM_metadata:=TRUE]
 }else if(study=='RPE1'){
  md<-fread(file.path(project,'paper19opp-revision2/reconstruction/cross_dataset/phase2_upr_rpe1/RPE1_author_roster_metadata_only.tsv.gz'))
  md[,context:='RPE1'];md[,gem:=as.character(gem_group)];md[,technical_group_type:='GEM'];md[,genuine_GEM_metadata:=TRUE]
 }else{
  f<-hdf5r::H5File$new(prep$paths[[study]],'r');on.exit(f$close_all())
  if(study=='Adamson'){
   md<-data.table(source_row=seq_along(annotation(f,'obs','_index')),cell_id=annotation(f,'obs','_index'),gene=annotation(f,'obs','perturbation'),guide=annotation(f,'obs','guide_identity'),UMI_count=as.numeric(annotation(f,'obs','UMI count')))
   md[,`:=`(context='Adamson_K562_UPR',is_control=gene=='63(mod)',gem=sub('.*-','',cell_id),technical_group_type='capture_lane',genuine_GEM_metadata=FALSE)]
   stopifnot(uniqueN(md$gem)==10L)
  }else{
   md<-data.table(source_row=seq_along(annotation(f,'obs','cell_name')),cell_id=annotation(f,'obs','cell_name'),gene=annotation(f,'obs','perturbation'),context=annotation(f,'obs','perturbation_2'),guide=annotation(f,'obs','guide_id'),sgRNA=annotation(f,'obs','sgRNA'),UMI_count=as.numeric(annotation(f,'obs','ncounts')),n_genes_by_counts=as.numeric(annotation(f,'obs','ngenes')),mitopercent=as.numeric(annotation(f,'obs','percent_mito')))
   md[,`:=`(is_control=gene=='control',gem='pooled',technical_group_type='pooled_context',genuine_GEM_metadata=FALSE)]
  }
 }
 md<-md[is_control==TRUE|gene%in%roster];setorder(md,source_row);md[,study:=study_name]
 stopifnot(!anyDuplicated(md$cell_id),all(is.finite(md$UMI_count)&md$UMI_count>0))
 md
}
