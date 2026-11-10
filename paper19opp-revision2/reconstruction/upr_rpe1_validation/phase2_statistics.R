# Locked phase-2 inference. Source-only ranks, memberships and matched sets are inputs.
# Technical intervals use shared-GEM Rademacher wild weights and fixed control SDs.
# A missing gene/GEM contributes no residual term; it is never an observed zero effect.
suppressPackageStartupMessages(library(data.table))
phase2_quantile <- function(x) {
  y <- x[is.finite(x)]
  if(length(y)<100L) return(c(NA_real_,NA_real_))
  as.numeric(quantile(y,c(.025,.975)))
}
phase2_stratified_permutation <- function(strata) {
  ii <- seq_along(strata)
  for(jj in split(ii,as.character(strata),drop=TRUE)) {
    if(length(jj)>1L) ii[jj] <- jj[sample.int(length(jj),length(jj),replace=FALSE)]
  }
  ii
}
phase2_spearman <- function(x,y) {
  if(any(!is.finite(x))||any(!is.finite(y))||length(x)<3L) return(NA_real_)
  suppressWarnings(cor(x,y,method='spearman'))
}
phase2_make_sets <- function(matches,ranking,k=20L) {
  m <- as.data.table(matches);r <- as.data.table(ranking)[order(rank)]
  top <- r[rank<=k,gene]
  stopifnot(length(top)==k,uniqueN(top)==k)
  sets <- lapply(top,function(g) c(g,sort(m[selected_gene==g,background_gene])))
  stopifnot(all(lengths(sets)==6L),uniqueN(unlist(sets))==6L*k)
  names(sets) <- top
  sets
}
phase2_set_contrasts <- function(y,sets) {
  vapply(sets,function(g) {
    z <- y[g]
    if(length(z)!=6L||any(!is.finite(z))) return(NA_real_)
    z[1]-mean(z[-1])
  },numeric(1))
}
phase2_set_permutation <- function(y,sets,roles,k) {
  vals <- vapply(sets, function(g) as.numeric(y[g]), numeric(6))
  if(any(!is.finite(vals))) return(rep(NA_real_,nrow(roles)))
  null <- matrix(0,nrow(roles),k)
  for(j in seq_len(k)) null[,j] <- 1.2*vals[cbind(roles[,j],j)]-sum(vals[,j])/5
  rowMeans(null)
}
phase2_wild_effects <- function(gem_table,genes,point,B=2000L) {
  d <- as.data.table(gem_table)[gene%in%genes&target_cells>0]
  stopifnot(all(is.finite(d$contrast_z)),all(is.finite(d$target_cells)),all(d$target_cells>0))
  # A GEM table may have repeated guide rows; collapse with observed-cell weights.
  d <- d[,.(target_cells=sum(target_cells),contrast_z=weighted.mean(contrast_z,target_cells)),by=.(gene,gem)]
  chk <- d[,.(delta_z=weighted.mean(contrast_z,target_cells),n=sum(target_cells)),by=gene]
  stopifnot(setequal(chk$gene,genes))
  p <- setNames(chk$delta_z,chk$gene)[genes]
  stopifnot(max(abs(p-as.numeric(point[genes])))<1e-8)
  total <- setNames(chk$n,chk$gene)
  d[,weighted_residual:=target_cells/total[gene]*(contrast_z-p[gene])]
  kk <- sort(unique(as.character(d$gem)))
  R <- matrix(0,length(genes),length(kk),dimnames=list(genes,kk))
  R[cbind(match(d$gene,genes),match(as.character(d$gem),kk))] <- d$weighted_residual
  # Every technical GEM receives one shared sign for all its genes and references.
  signs <- matrix(sample(c(-1,1),length(kk)*B,replace=TRUE),nrow=length(kk),ncol=B)
  values <- sweep(R%*%signs,1,p,'+')
  dimnames(values) <- list(genes,NULL)
  list(effects=values,gem_labels=kk,wild_signs=signs,point=p)
}
phase2_primary <- function(source,validation,ranking,matches,source_gem,validation_gem,
                           B=2000L,P=10000L,seed=20261106L) {
  set.seed(seed)
  s <- as.data.table(source);v <- as.data.table(validation);r <- as.data.table(ranking)
  stopifnot(!anyDuplicated(s$gene),!anyDuplicated(v$gene),!anyDuplicated(r$gene))
  genes <- s$gene
  stopifnot(setequal(genes,r$gene),all(genes%in%v$gene),all(!is.na(s$stratum)))
  x <- setNames(s$source_delta_z,s$gene)[genes]
  score <- setNames(s$priority_score,s$gene)[genes]
  y <- setNames(v$validation_delta_z,v$gene)[genes]
  stopifnot(sum(is.finite(score)&score>0)>=20L)
  sets <- phase2_make_sets(matches,r,20L)
  stopifnot(all(unlist(sets)%in%genes))
  observed_rho <- phase2_spearman(x,y)
  observed_gap <- phase2_set_contrasts(y,sets)
  rho_null <- vapply(seq_len(P),function(b)phase2_spearman(x,y[phase2_stratified_permutation(s$stratum)]),numeric(1))
  roles <- matrix(sample.int(6L,P*20L,replace=TRUE),nrow=P,ncol=20L)
  null20 <- phase2_set_permutation(y,sets,roles,20L)
  null10 <- phase2_set_permutation(y,sets,roles,10L)
  obs <- c(observed_rho,mean(observed_gap),mean(observed_gap[1:10]))
  nulls <- list(rho_null,null20,null10)
  pv <- vapply(seq_len(3),function(i) {
    if(!is.finite(obs[i])||any(!is.finite(nulls[[i]])))return(NA_real_)
    (1+sum(abs(nulls[[i]])>=abs(obs[i])-1e-12))/(P+1)
  },numeric(1))
  holm <- p.adjust(ifelse(is.finite(pv),pv,1),'holm',n=3);holm[!is.finite(pv)]<-NA_real_
  gene_ids <- replicate(B,sample.int(length(genes),length(genes),replace=TRUE),simplify=FALSE)
  set20_ids <- replicate(B,sample.int(20L,20L,replace=TRUE),simplify=FALSE)
  set10_ids <- replicate(B,sample.int(10L,10L,replace=TRUE),simplify=FALSE)
  brho <- vapply(gene_ids,function(ii)phase2_spearman(x[ii],y[ii]),numeric(1))
  b20 <- vapply(set20_ids,function(ii)mean(observed_gap[ii]),numeric(1))
  b10 <- vapply(set10_ids,function(ii)mean(observed_gap[ii]),numeric(1))
  sw <- phase2_wild_effects(source_gem,genes,x,B)
  vw <- phase2_wild_effects(validation_gem,genes,y,B)
  wrho <- xrho <- w20 <- x20 <- w10 <- x10 <- rep(NA_real_,B)
  for(b in seq_len(B)) {
    xb <- sw$effects[,b];yb <- setNames(vw$effects[,b],genes)
    wrho[b] <- phase2_spearman(xb,yb)
    ii <- gene_ids[[b]];xrho[b] <- phase2_spearman(xb[ii],yb[ii])
    gaps <- phase2_set_contrasts(yb,sets)
    w20[b] <- mean(gaps);w10[b]<-mean(gaps[1:10])
    x20[b] <- mean(gaps[set20_ids[[b]]]);x10[b]<-mean(gaps[set10_ids[[b]]])
  }
  qb <- lapply(list(brho,b20,b10),phase2_quantile)
  qw <- lapply(list(wrho,w20,w10),phase2_quantile)
  qx <- lapply(list(xrho,x20,x10),phase2_quantile)
  result <- data.table(endpoint=c('Signed source–RPE1 Spearman','Top20 matched mean gap','Top10 matched mean gap'),
    genes=c(length(genes),120L,60L),estimate=obs,two_sided_permutation_p=pv,Holm_3_p=holm,
    positive_support=is.finite(obs)&obs>0&is.finite(holm)&holm<.05,
    gene_bootstrap_low=vapply(qb,'[',numeric(1),1),gene_bootstrap_high=vapply(qb,'[',numeric(1),2),
    technical_GEM_wild_low=vapply(qw,'[',numeric(1),1),technical_GEM_wild_high=vapply(qw,'[',numeric(1),2),
    crossed_gene_GEM_low=vapply(qx,'[',numeric(1),1),crossed_gene_GEM_high=vapply(qx,'[',numeric(1),2),
    permutation_draws=P,bootstrap_draws=B)
  diag <- data.table(genes=length(genes),source_positive=sum(x>0),source_score_positive=sum(score>0),
    score_spearman=phase2_spearman(score,y),source_GEMs=length(sw$gem_labels),validation_GEMs=length(vw$gem_labels),
    Top20_mean=mean(y[r[order(rank)][1:20,gene]]),Top20_background_mean=mean(y[unlist(sets)[rep(c(FALSE,rep(TRUE,5)),20)]]),
    Top10_mean=mean(y[r[order(rank)][1:10,gene]]),Top10_background_mean=mean(y[unlist(sets[1:10])[rep(c(FALSE,rep(TRUE,5)),10)]]))
  list(endpoints=result,diagnostics=diag,sets=sets,point_source=x,point_validation=y,
       source_wild=sw,validation_wild=vw,gene_ids=gene_ids,
       set20_ids=set20_ids,set10_ids=set10_ids,permutation_roles=roles,
       permutation_null=data.table(draw=seq_len(P),signed_rho=rho_null,Top20_gap=null20,Top10_gap=null10),
       bootstrap=data.table(draw=seq_len(B),gene_rho=brho,gene_Top20=b20,gene_Top10=b10,
          wild_rho=wrho,wild_Top20=w20,wild_Top10=w10,crossed_rho=xrho,crossed_Top20=x20,crossed_Top10=x10))
}
phase2_baselines <- function(primary,source,comparators,common_genes,raw_upr_method,
                             B=2000L,P=10000L,seed=20261107L) {
  set.seed(seed)
  s<-as.data.table(source);cc<-as.character(common_genes)
  cmp<-as.data.table(comparators)[gene%in%cc];methods<-unique(cmp$method)
  stopifnot(raw_upr_method%in%methods,all(cc%in%names(primary$point_validation)))
  for(m in methods)stopifnot(setequal(cmp[method==m,gene],cc))
  y<-primary$point_validation[cc];strata<-setNames(s$stratum,s$gene)[cc]
  ids<-replicate(B,sample.int(length(cc),length(cc),replace=TRUE),simplify=FALSE)
  perms<-replicate(P,phase2_stratified_permutation(strata),simplify=FALSE)
  rows<-list()
  for(k in c(20L,10L)) {
    if(length(cc)<=k)next
    rawset<-cmp[method==raw_upr_method][order(rank)][1:k,gene]
    rawsel<-cc%in%rawset
    for(m in methods) {
      selset<-cmp[method==m][order(rank)][1:k,gene];sel<-cc%in%selset
      diff<-mean(y[sel])-mean(y[rawsel])
      boot<-vapply(ids,function(ii){u<-sel[ii];v<-rawsel[ii];if(!any(u)||!any(v))return(NA_real_);mean(y[ii][u])-mean(y[ii][v])},numeric(1))
      null<-vapply(perms,function(ii)mean(y[ii][sel])-mean(y[ii][rawsel]),numeric(1))
      p<-if(m==raw_upr_method)NA_real_ else(1+sum(abs(null)>=abs(diff)-1e-12))/(P+1)
      wild<-cross<-rep(NA_real_,B)
      for(b in seq_len(B)) {
        yy<-primary$validation_wild$effects[cc,b]
        wild[b]<-mean(yy[sel])-mean(yy[rawsel]);ii<-ids[[b]];u<-sel[ii];v<-rawsel[ii]
        if(any(u)&&any(v))cross[b]<-mean(yy[ii][u])-mean(yy[ii][v])
      }
      qb<-phase2_quantile(boot);qw<-phase2_quantile(wild);qx<-phase2_quantile(cross)
      rows[[length(rows)+1L]]<-data.table(method=m,TopK=k,complete_case_genes=length(cc),mean_RPE1_UPR=mean(y[sel]),
        minus_raw_UPR_mean=diff,gene_bootstrap_low=qb[1],gene_bootstrap_high=qb[2],technical_GEM_wild_low=qw[1],technical_GEM_wild_high=qw[2],
        crossed_gene_GEM_low=qx[1],crossed_gene_GEM_high=qx[2],two_sided_permutation_p=p,selection_overlap_raw=sum(sel&rawsel),selected_genes=paste(selset,collapse=';'))
    }
  }
  r<-rbindlist(rows);r[,exploratory_Holm_p:=NA_real_];r[is.finite(two_sided_permutation_p),exploratory_Holm_p:=p.adjust(two_sided_permutation_p,'holm')]
  r
}
phase2_module_specificity <- function(primary,source,validation,null_source,null_validation,null_matches,null_rankings) {
  s<-as.data.table(source);ns<-as.data.table(null_source);nv<-as.data.table(null_validation)
  nr<-as.data.table(null_rankings);nm<-as.data.table(null_matches);genes<-s$gene;y<-primary$point_validation[genes]
  ids<-sort(unique(as.character(ns$module_id)));stopifnot(length(ids)==100L)
  rows<-lapply(ids,function(id) {
    a<-ns[as.character(module_id)==id];b<-nv[as.character(module_id)==id]
    x<-setNames(a$source_delta_z,a$gene)[genes];v<-setNames(b$validation_delta_z,b$gene)[genes]
    r<-nr[as.character(module_id)==id];m<-nm[as.character(module_id)==id]
    stopifnot('priority_score'%in%names(r));valid_top<-nrow(r)>=20L&&all(1:20%in%r$rank)&&sum(is.finite(r$priority_score)&r$priority_score>0)>=20L&&uniqueN(m$selected_gene)==20L
    gap<-topmean<-NA_real_
    if(valid_top) {
      sets<-phase2_make_sets(m,r,20L);gap<-mean(phase2_set_contrasts(y,sets));topmean<-mean(y[r[order(rank)][1:20,gene]])
    }
    data.table(module_id=id,positive_source_genes=sum(is.finite(r$priority_score)&r$priority_score>0),genes=sum(is.finite(x)&is.finite(v)),signed_rho=phase2_spearman(x,v),Top20_on_real_UPR_mean=topmean,
               Top20_on_real_UPR_matched_gap=gap,valid_Top20=valid_top)
  })
  r<-rbindlist(rows)
  observed<-c(primary$endpoints$estimate[1],primary$endpoints$estimate[2])
  fields<-c('signed_rho','Top20_on_real_UPR_matched_gap')
  sums<-rbindlist(lapply(seq_along(fields),function(i){z<-r[[fields[i]]];z<-z[is.finite(z)];data.table(calibration=fields[i],observed=observed[i],valid_null_modules=length(z),
    upper_tail_fraction=(1+sum(z>=observed[i]-1e-12))/(1+length(z)),empirical_percentile=(sum(z<observed[i])+0.5*sum(abs(z-observed[i])<1e-12))/length(z),
    null_median=median(z),null_95_low=quantile(z,.025),null_95_high=quantile(z,.975))}))
  sums[,specificity_Holm2:=p.adjust(upper_tail_fraction,'holm')]
  list(all_modules=r,calibration=sums)
}
