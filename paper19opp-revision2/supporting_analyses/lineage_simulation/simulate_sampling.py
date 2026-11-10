#!/usr/bin/env python3
"""Known-null lineage geometry experiment; conditional on recorded baseline design.
No independent biological replication is asserted. See README.txt for assumptions.
"""
from pathlib import Path
import csv, json, time
import numpy as np
OUT = Path(__file__).resolve().parent
SEED=20261106
REPLICATES=200
SIGMAS=(0.5,1.0,2.0)

def read(path):
    with path.open() as f: return list(csv.DictReader(f,delimiter='\t'))
def groups_neighbors(path):
    rows=read(path); n=max(int(x['row']) for x in rows)
    a=[[] for _ in range(n)]
    for x in rows:a[int(x['row'])-1].append(int(x['neighbor'])-1)
    groups=[]
    for k in sorted(set(map(len,a))):
        q=np.array([i for i,x in enumerate(a) if len(x)==k])
        groups.append((q,np.array([a[i] for i in q])))
    return groups

def rank(y,groups):
    result=np.empty(len(y))
    for q,idx in groups:
        yy=y[idx];centered=yy-yy.mean(axis=1,keepdims=True)
        cov=np.einsum('nkb,nkc->nbc',centered,centered)/(idx.shape[1]-1)
        ev=np.maximum(np.linalg.eigvalsh(cov),0)
        total=ev.sum(axis=1,keepdims=True)
        ev/=np.maximum(total,1e-300)
        result[q]=np.where(total[:,0]>0,np.exp(-np.sum(ev*np.log(np.maximum(ev,1e-300)),axis=1)),0)
    return result

def multinomial(p,n,rng):
    # Sequential binomial construction supports row-specific sample size.
    remain=n.copy();den=np.ones(len(p));ans=np.zeros_like(p)
    for j in range(p.shape[1]-1):
        pp=np.clip(p[:,j]/np.maximum(den,1e-300),0,1)
        ans[:,j]=rng.binomial(remain,pp)
        remain-=ans[:,j].astype(int);den-=p[:,j]
    ans[:,-1]=remain
    return ans/n[:,None]

def write(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

start=time.time();results=[];calibration=[]
for di,nm in enumerate(('larry_cytokine','larry_in_vitro')):
    d=read(OUT/f'{nm}_design.tsv');branches=read(OUT/f'{nm}_branches.tsv');columns=[k for k in branches[0] if k!='row']
    y=np.array([[float(x[k]) for k in columns] for x in branches])
    n=np.array([int(x['descendants']) for x in d]);sr=np.array([float(x['stress_residual']) for x in d]);den=(sr*sr).sum()
    condition=np.array([x['condition'] for x in d]);clones=np.array([x['clone_id'] for x in d]);ids,cloneindex=np.unique(clones,return_inverse=True)
    nn=groups_neighbors(OUT/f'{nm}_neighbors.tsv');deep=read(OUT/f'{nm}_deep_design.tsv');sel=np.array([int(x['original_row'])-1 for x in deep]);dsr=np.array([float(x['stress_residual']) for x in deep]);dden=(dsr*dsr).sum();dnn=groups_neighbors(OUT/f'{nm}_deep_neighbors.tsv')
    # Condition-specific offsets use descendant-weighted empirical branch frequencies
    # with a half-count pseudocount. They do not estimate unobserved clone truth.
    offset=np.empty_like(y);strata=[]
    for co in np.unique(condition):
        ix=np.where(condition==co)[0];strata.append(ix)
        counts=(y[ix]*n[ix,None]).sum(axis=0)+0.5;pi=counts/counts.sum()
        offset[ix]=np.log(pi)
        calibration += [dict(dataset=nm,condition=co,branch=b,offset_frequency=float(v),units=len(ix),total_descendants=int(n[ix].sum())) for b,v in zip(columns,pi)]
    reconstructed=rank(y,nn)
    archived=np.array([float(x["observed_rank"]) for x in d])
    assert np.max(np.abs(reconstructed-archived))<1e-8, (nm,np.max(np.abs(reconstructed-archived)))
    obs=float(np.sum(sr*reconstructed)/den)
    for sigma in SIGMAS:
      for rep in range(REPLICATES):
        rng=np.random.default_rng(np.random.SeedSequence([SEED,di,int(sigma*10),rep]))
        # One Gaussian vector per clone is shared across its recorded conditions.
        z=rng.normal(size=(len(ids),y.shape[1]))[cloneindex]
        logits=offset+sigma*z;logits-=logits.max(axis=1,keepdims=True)
        p=np.exp(logits);p/=p.sum(axis=1,keepdims=True)
        observed=multinomial(p,n,rng);equal=multinomial(p,np.repeat(10,len(n)),rng)
        pn=n.copy()
        for ix in strata:pn[ix]=rng.permutation(n[ix])
        permuted=multinomial(p,pn,rng)
        # Random donor neighborhoods preserve condition, k and exclusion of self,
        # while discarding baseline-neighbor structure. Same random graph for oracle
        # and noisy probabilities in this replicate.
        randomnn=[]
        for q,idx in nn:
            randomidx=[]
            for row in q:
                pool=np.where(condition==condition[row])[0];pool=pool[pool!=row]
                randomidx.append(rng.choice(pool,idx.shape[1],replace=False))
            randomnn.append((q,np.array(randomidx)))
        scenarios={'oracle_probability':rank(p,nn),
                   'recorded_depth':rank(observed,nn),
                   'equal_depth_10':rank(equal,nn),
                   'permuted_depth':rank(permuted,nn),
                   'random_neighbors_oracle':rank(p,randomnn),
                   'random_neighbors_recorded_depth':rank(observed,randomnn)}
        for scenario,ranks in scenarios.items():
            results.append(dict(dataset=nm,latent_logit_sd=sigma,replicate=rep+1,scenario=scenario,units=len(n),slope=float(np.sum(sr*ranks)/den),median_rank=float(np.median(ranks)),observed_slope=obs))
        for scenario,yy in [('deep_oracle_probability',p[sel]),('deep_recorded_depth',observed[sel]),('deep_equal_depth_10',equal[sel])]:
            ranks=rank(yy,dnn)
            results.append(dict(dataset=nm,latent_logit_sd=sigma,replicate=rep+1,scenario=scenario,units=len(sel),slope=float(np.sum(dsr*ranks)/dden),median_rank=float(np.median(ranks)),observed_slope=obs))
      print(nm,sigma,'complete',round(time.time()-start,1),'seconds',flush=True)
write(OUT/'simulation_replicates.tsv',results);write(OUT/'simulation_calibration.tsv',calibration)
summary=[]
for nm in ('larry_cytokine','larry_in_vitro'):
 for sigma in SIGMAS:
  for scenario in sorted(set(x['scenario'] for x in results)):
   v=[x for x in results if x['dataset']==nm and x['latent_logit_sd']==sigma and x['scenario']==scenario];b=np.array([x['slope'] for x in v]);r=np.array([x['median_rank'] for x in v]);obs=v[0]['observed_slope']
   summary.append(dict(dataset=nm,latent_logit_sd=sigma,scenario=scenario,replicates=len(b),units=v[0]['units'],median_slope=float(np.median(b)),slope_q025=float(np.quantile(b,.025)),slope_q975=float(np.quantile(b,.975)),negative_fraction=float(np.mean(b<0)),median_rank=float(np.median(r)),observed_slope=obs,null_slopes_at_or_below_observed=int(np.sum(b<=obs))))
write(OUT/'simulation_summary.tsv',summary)
# Paired sampling-minus-oracle contrasts isolate observation bias on identical latent draws.
paired=[]
for nm in ('larry_cytokine','larry_in_vitro'):
 for sigma in SIGMAS:
  vv=[x for x in results if x['dataset']==nm and x['latent_logit_sd']==sigma];lookup={(x['replicate'],x['scenario']):x['slope'] for x in vv}
  for target,ref in [('recorded_depth','oracle_probability'),('equal_depth_10','oracle_probability'),('permuted_depth','oracle_probability'),('random_neighbors_recorded_depth','random_neighbors_oracle'),('deep_recorded_depth','deep_oracle_probability'),('deep_equal_depth_10','deep_oracle_probability')]:
   b=np.array([lookup[(i,target)]-lookup[(i,ref)] for i in range(1,REPLICATES+1)])
   paired.append(dict(dataset=nm,latent_logit_sd=sigma,comparison=f'{target}_minus_{ref}',replicates=len(b),median_difference=float(np.median(b)),q025=float(np.quantile(b,.025)),q975=float(np.quantile(b,.975))))
write(OUT/'simulation_paired_bias.tsv',paired)
(OUT/'metadata.json').write_text(json.dumps(dict(seed=SEED,replicates=REPLICATES,latent_logit_sd=SIGMAS,numpy=np.__version__,elapsed_seconds=time.time()-start),indent=2)+'\n')
print('All complete',round(time.time()-start,1),'seconds',flush=True)
