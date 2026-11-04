#!/usr/bin/env python3
"""Train genuine nested scVI references; query encoding never updates weights."""
from __future__ import annotations
import argparse, hashlib, importlib.metadata, json, os, time, sys, subprocess, tempfile, urllib.request, zipfile
from pathlib import Path
TASK = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'paper19_scvi_cache/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir())/'paper19_scvi_cache'))
os.environ.setdefault('NUMBA_CACHE_DIR', str(Path(tempfile.gettempdir())/'paper19_scvi_cache/numba'))
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
import anndata as ad
import scanpy as sc
import scvi
import torch
from filelock import FileLock
from lightning.pytorch.callbacks import Callback
from lightning.pytorch.plugins.environments import LightningEnvironment

SPEC = dict(n_hidden=128,n_latent=30,n_layers=2,dropout_rate=0.1,dispersion='gene',gene_likelihood='nb',use_observed_lib_size=True,latent_distribution='normal')
SOURCE_NAMES=['DMSO_24hr_expt1.qc.h5ad','DMSO_24hr_expt3.qc.h5ad','DMSO_6hr_expt1.qc.h5ad','DMSO_6hr_expt3.qc.h5ad','DMSO_expt10.qc.h5ad']

def save_json(v,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+f'.{os.getpid()}.tmp');t.write_text(json.dumps(v,indent=2,default=str));t.replace(p)

def key(ex,tm,ln):
    def f(x):return 'NA' if pd.isna(x) else format(float(x),'g')
    return [f'{f(e)}|{f(t)}|{l}' for e,t,l in zip(ex,tm,ln)]

def matrix_hash(x):
    x=sparse.csr_matrix(x);h=hashlib.sha256();h.update(np.asarray(x.shape,dtype='<i8').tobytes())
    for values in [x.indptr,x.indices,x.data]:
        for start in range(0,len(values),1048576):h.update(values[start:start+1048576].astype('<i8').tobytes())
    return h.hexdigest()

def weights_hash(model):
    h=hashlib.sha256()
    for n,t in sorted(model.module.state_dict().items()):
        h.update(n.encode());h.update(str(tuple(t.shape)).encode());h.update(t.detach().cpu().contiguous().numpy().tobytes())
    return h.hexdigest()

def versions():
    return {p:importlib.metadata.version(p) for p in ['scvi-tools','scanpy','anndata','torch','lightning','numpy','scipy','pandas','h5py']}

EXPECTED_VERSIONS={'scvi-tools':'1.4.3','scanpy':'1.11.5','anndata':'0.12.19','torch':'2.5.1','lightning':'2.6.6','numpy':'1.26.4','scipy':'1.13.1','pandas':'2.2.2','h5py':'3.16.0'}
COUNT_SHA='518bf5f91c9250b00ad115c0d4c3baccb8324a707665894a3a6320b70c6fe1ea'

def compatible_versions(actual):
    normalized={k:(v.split('+')[0] if k=='torch' else v) for k,v in actual.items()}
    return normalized==EXPECTED_VERSIONS

def check_runtime():
    actual=versions()
    assert compatible_versions(actual), f'Use the pinned clean environment: {actual}'
    return actual

def data_fingerprints(a,manifest_path=None):
    assert a.shape==(21936,32738) and a.obs_names.is_unique and a.var_names.is_unique
    assert a.obs.cell_line.nunique()==169 and a.obs.technical_batch.nunique()==5
    assert np.isfinite(a.X.data).all() and (a.X.data>=0).all() and np.equal(a.X.data,np.floor(a.X.data)).all()
    digest=matrix_hash(a.X);assert digest==COUNT_SHA
    columns=['source_file','original_cell_id','cell_line','technical_batch','baseline_context']
    metadata=a.obs[columns].astype(str).to_csv(sep='\t',index=True)
    config=dict(architecture=SPEC,seed=20260713,hvg_flavor='cell_ranger',hvg_n=4000,hvg_batch='technical_batch',normalization_target=10000,max_epochs=200,early_stopping_patience=20,validation_every=5,batch_size=256)
    return dict(count_matrix_sha256=digest,cell_metadata_sha256=hashlib.sha256(metadata.encode()).hexdigest(),feature_order_sha256=hashlib.sha256('\n'.join(a.var_names).encode()).hexdigest(),fold_manifest_sha256=hashlib.sha256(Path(manifest_path).read_bytes()).hexdigest() if manifest_path is not None else None,configuration_sha256=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest())

def validate_completed_fit(out,m,fingerprints,allow_legacy=False):
    out=Path(out);audit=json.loads((out/'complete.json').read_text());adopted_fingerprints='fingerprints' not in audit
    assert audit['architecture']==SPEC and audit['seed']==20260713 and audit['max_epochs']==200 and audit['early_stopping_patience']==20 and audit['validation_frequency_epochs']==5 and audit['batch_size']==256
    assert compatible_versions(audit['versions'])
    assert audit['weight_hash_before_query']==audit['weight_hash_after_query'] and audit['query_updates_weights'] is False
    assert audit['training_batch_count']==5 and audit['query_unseen_batches']==0
    assert 0<audit['actual_epochs']<=200
    required=['epoch_log.jsonl','run_fingerprints.json','task_latent_means.tsv','selected_genes.tsv','training_hvg_statistics.tsv.gz','cell_membership.tsv.gz','history_elbo_train.tsv','history_elbo_validation.tsv','model/model.pt']
    if int(m.inner_fold.iloc[0])==0:required+=['all_baseline_context_latent_means.tsv']
    for name in required:assert (out/name).is_file(), f'Missing completed-fit output: {out/name}'
    provenance=json.loads((out/'run_fingerprints.json').read_text());assert provenance['fingerprints']==fingerprints and provenance['fit_id']==audit['fit_id'] and compatible_versions(provenance['versions'])
    part=pd.read_csv(out/'task_latent_means.tsv',sep='\t',float_precision='round_trip');want=m[m.task_role.isin(['train','validation','test'])].sort_values('row_id')
    for col in ['row_id','pair_id','cell_line','baseline_context','fit_id','task_role']:assert list(part[col])==list(want[col]),col
    assert np.isfinite(part[[f'z{i}_mean' for i in range(30)]].to_numpy()).all()
    cells=pd.read_csv(out/'cell_membership.tsv.gz',sep='\t',index_col=0);held=set(m.loc[m.task_role.ne('train'),'cell_line']);query=set(m.loc[m.task_role.isin(['validation','test']),'cell_line'])
    expected=np.where(~cells.cell_line.isin(held),'reference_fit',np.where(cells.cell_line.isin(query),'frozen_query','outer_test_excluded'))
    assert len(cells)==21936 and cells.index.is_unique and np.array_equal(expected,cells.cell_role)
    assert int((cells.cell_role=='reference_fit').sum())==audit['training_cells'] and int((cells.cell_role=='frozen_query').sum())==audit['query_cells']
    internal=set();reference=set(cells.index[cells.cell_role.eq('reference_fit')])
    for name,count in audit['internal_split_counts'].items():
        assert count is not None
        file=out/f'internal_{name}.tsv.gz';assert file.is_file()
        assignments=pd.read_csv(file,sep='\t');ids=set(assignments.cell_uid)
        assert len(assignments)==count and len(ids)==count and ids<=reference and not ids&internal
        internal.update(ids)
    assert internal==reference
    for name in ['elbo_train','elbo_validation']:
        history=pd.read_csv(out/f'history_{name}.tsv',sep='\t',index_col=0)
        assert len(history)>0 and np.isfinite(history.to_numpy(dtype=float)).all()
    epoch_records=[json.loads(x) for x in (out/'epoch_log.jsonl').read_text().splitlines()]
    assert epoch_records and max(x['epoch'] for x in epoch_records)==audit['actual_epochs']
    genes=pd.read_csv(out/'selected_genes.tsv',sep='\t');assert len(genes)==4000 and genes.gene.is_unique
    hvg=pd.read_csv(out/'training_hvg_statistics.tsv.gz',sep='\t',index_col=0);assert list(hvg.index[hvg.highly_variable])==list(genes.gene)
    if 'fingerprints' in audit:assert audit['fingerprints']==fingerprints,'Stale input/configuration fingerprint'
    else:
        assert allow_legacy, 'Missing recorded fingerprints; explicit independent audit is required'
        # Adoption requires the source, membership, configuration, runtime
        # and required-output checks above.
        audit['fingerprints']=fingerprints;audit['fingerprints_independently_verified_after_training']=True
    if int(m.inner_fold.iloc[0])==0:
        contexts=pd.read_csv(out/'all_baseline_context_latent_means.tsv',sep='\t',float_precision='round_trip')
        assert contexts.baseline_context.is_unique and set(contexts.baseline_context)==set(cells.baseline_context)
        assert np.isfinite(contexts[[f'z{i}_mean' for i in range(30)]].to_numpy()).all()
        expected_n=cells.groupby('baseline_context',observed=True).size()
        assert np.array_equal(contexts.baseline_ncells,expected_n.loc[contexts.baseline_context])
        aligned=contexts.set_index('baseline_context').loc[part.baseline_context,[f'z{i}_mean' for i in range(30)]].to_numpy()
        assert np.array_equal(aligned,part[[f'z{i}_mean' for i in range(30)]].to_numpy())
    if adopted_fingerprints or not audit.get('completed_outputs_validated',False):
        audit['completed_outputs_validated']=True;save_json(audit,out/'complete.json')
    return audit


def export_context_means(means,a,out):
    meta=a.obs[['baseline_context','cell_line','experiment','time_hours']].drop_duplicates()
    assert not meta.baseline_context.duplicated().any()
    meta=meta.merge(a.obs.groupby('baseline_context',observed=True).size().rename('baseline_ncells'),on='baseline_context',validate='one_to_one')
    table=means.reset_index().merge(meta,on='baseline_context',how='left',validate='one_to_one')
    assert table.cell_line.notna().all()
    table.to_csv(Path(out)/'all_baseline_context_latent_means.tsv',sep='\t',index=False,float_format='%.17g')

class EpochLog(Callback):
    def __init__(self,path,checkpoint):
        self.path=Path(path);self.started=time.time();self.checkpoint=Path(checkpoint);self.checkpoint.parent.mkdir(parents=True,exist_ok=True)
    def on_train_epoch_end(self,trainer,pl_module):
        e=int(trainer.current_epoch)+1
        if e==1 or e%5==0:
            vals={k:float(v.detach().cpu()) if torch.is_tensor(v) and v.numel()==1 else str(v) for k,v in trainer.callback_metrics.items()}
            with self.path.open('a') as h:h.write(json.dumps(dict(epoch=e,elapsed_seconds=time.time()-self.started,metrics=vals))+'\n')
            print(f'epoch={e} elapsed_seconds={time.time()-self.started:.1f}',flush=True)
        if e%10==0:trainer.save_checkpoint(str(self.checkpoint))

def run_fit(input_path,manifest_path,fit_id,out,threads=4):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    # Concurrent resume queues may target a completed or still-active fit.
    # The OS-backed lock is released on process exit, including interruption.
    with FileLock(str(out/'.fit.lock')):
        return _run_fit(input_path,manifest_path,fit_id,out,threads)

def _run_fit(input_path,manifest_path,fit_id,out,threads=4):
    started=time.time();out=Path(out);out.mkdir(parents=True,exist_ok=True)
    check_runtime()
    torch.set_num_threads(threads)
    try:torch.set_num_interop_threads(1)
    except RuntimeError:pass
    scvi.settings.seed=20260713
    scvi.settings.dl_num_workers=0
    a=ad.read_h5ad(input_path);m=pd.read_csv(manifest_path,sep='\t',float_precision='round_trip');m=m[m.fit_id.eq(fit_id)].copy()
    assert len(m)==120 and not m.row_id.duplicated().any()
    fingerprints=data_fingerprints(a,manifest_path)
    ckpt=out/'checkpoints/last.ckpt';provenance=out/'run_fingerprints.json'
    if ckpt.exists() or provenance.exists() or (out/'complete.json').exists():
        assert provenance.exists(), 'Checkpoint has no input/configuration provenance; refuse automatic recovery'
        previous=json.loads(provenance.read_text())
        assert previous['fingerprints']==fingerprints and previous['fit_id']==fit_id and compatible_versions(previous['versions'])
    else:save_json(dict(fit_id=fit_id,fingerprints=fingerprints,versions=versions(),recorded_before_training=True),provenance)
    kw={'ckpt_path':str(ckpt)} if ckpt.exists() else {}
    if (out/'complete.json').exists():
        validate_completed_fit(out,m,fingerprints);print('Validated complete:',fit_id,flush=True);return
    test_lines=set(m.loc[m.task_role.isin(['test','validation','outer_test_excluded']),'cell_line'])
    train_lines=set(m.loc[m.task_role.eq('train'),'cell_line']);query_lines=set(m.loc[m.task_role.isin(['test','validation']),'cell_line'])
    assert not train_lines&test_lines and not query_lines&train_lines
    train_mask=~a.obs.cell_line.isin(test_lines).to_numpy();query_mask=a.obs.cell_line.isin(query_lines).to_numpy()
    assert not np.any(train_mask&query_mask)
    assert set(a.obs.loc[train_mask,'cell_line']).isdisjoint(test_lines)
    assert set(m.loc[m.task_role.eq('train'),'baseline_context'])<=set(a.obs.loc[train_mask,'baseline_context'])
    assert set(a.obs.loc[train_mask,'technical_batch'].astype(str))==set(a.obs.technical_batch.astype(str))
    assert len(set(a.obs.technical_batch))==5
    cellrole=np.where(train_mask,'reference_fit',np.where(query_mask,'frozen_query','outer_test_excluded'))
    cellmem=a.obs.copy();cellmem['fit_id']=fit_id;cellmem['cell_role']=cellrole;cellmem.to_csv(out/'cell_membership.tsv.gz',sep='\t',index=True)
    tr=a[train_mask].copy();norm=tr.copy();sc.pp.normalize_total(norm,target_sum=1e4);sc.pp.log1p(norm)
    sc.pp.highly_variable_genes(norm,n_top_genes=min(4000,norm.n_vars),flavor='cell_ranger',batch_key='technical_batch')
    genes=norm.var_names[norm.var.highly_variable].copy();assert len(genes)==4000
    norm.var[['highly_variable','highly_variable_nbatches','means','dispersions','dispersions_norm']].to_csv(out/'training_hvg_statistics.tsv.gz',sep='\t',index=True)
    pd.DataFrame({'gene':genes}).to_csv(out/'selected_genes.tsv',sep='\t',index=False)
    del norm
    tr=tr[:,genes].copy();tr.layers['counts']=tr.X
    assert np.isfinite(tr.X.data).all() and (tr.X.data>=0).all() and np.equal(tr.X.data,np.floor(tr.X.data)).all()
    scvi.model.SCVI.setup_anndata(tr,layer='counts',batch_key='technical_batch')
    model=scvi.model.SCVI(tr,**SPEC)
    callbacks=[EpochLog(out/'epoch_log.jsonl',out/'checkpoints/last.ckpt')]
    model.train(max_epochs=200,accelerator='cpu',devices=1,early_stopping=True,early_stopping_patience=20,check_val_every_n_epoch=5,batch_size=256,enable_progress_bar=False,default_root_dir=str(out),callbacks=callbacks,logger=None,log_save_dir=str(out/'training_log'),plugins=[LightningEnvironment()],**kw)
    model.save(out/'model',overwrite=True,save_anndata=False)
    model.module.eval()
    before=weights_hash(model)
    # Reference posterior means and query posterior means use the same frozen encoder.
    zref=model.get_latent_representation(adata=tr,give_mean=True,batch_size=256)
    q=a[query_mask,genes].copy();q.layers['counts']=q.X
    assert q.var_names.equals(tr.var_names)
    assert set(q.obs.technical_batch.astype(str))<=set(tr.obs.technical_batch.astype(str))
    zquery=model.get_latent_representation(adata=q,give_mean=True,batch_size=256)
    after=weights_hash(model);assert before==after
    zr=pd.DataFrame(zref,index=tr.obs_names,columns=[f'z{i}_mean' for i in range(30)]);zq=pd.DataFrame(zquery,index=q.obs_names,columns=zr.columns)
    allz=pd.concat([zr,zq]);allobs=pd.concat([tr.obs[['baseline_context']],q.obs[['baseline_context']]])
    allz['baseline_context']=allobs.loc[allz.index,'baseline_context'].to_numpy();means=allz.groupby('baseline_context',sort=True).mean()
    if int(m.inner_fold.iloc[0])==0:export_context_means(means,a,out)
    used=m[m.task_role.isin(['train','validation','test'])].copy();assert set(used.baseline_context)<=set(means.index)
    result=used.merge(means,on='baseline_context',how='left',validate='many_to_one').sort_values('row_id')
    assert np.isfinite(result[[f'z{i}_mean' for i in range(30)]].to_numpy()).all()
    assert set(result.loc[result.task_role.eq('train'),'cell_line']).isdisjoint(set(result.loc[result.task_role.ne('train'),'cell_line']))
    result.to_csv(out/'task_latent_means.tsv',sep='\t',index=False,float_format='%.17g')
    for name,h in model.history.items():h.to_csv(out/f'history_{name}.tsv',sep='\t',index=True)
    indices={}
    for name in ['train_indices','validation_indices','test_indices']:
        try:
            idx=np.asarray(getattr(model,name),dtype=int);indices[name]=len(idx)
            pd.DataFrame({'cell_uid':tr.obs_names[idx],'internal_role':name}).to_csv(out/f'internal_{name}.tsv.gz',sep='\t',index=False)
            assert set(tr.obs.iloc[idx].cell_line).isdisjoint(test_lines)
        except AttributeError:indices[name]=None
    audit=dict(fit_id=fit_id,outer_fold=int(m.outer_fold.iloc[0]),inner_fold=int(m.inner_fold.iloc[0]),training_cells=tr.n_obs,query_cells=q.n_obs,excluded_outer_cells=int((cellrole=='outer_test_excluded').sum()),reference_lines=tr.obs.cell_line.nunique(),query_lines=q.obs.cell_line.nunique(),fixed_auxiliary_lines=len(set(a.obs.cell_line)-set(m.cell_line)),training_batch_count=tr.obs.technical_batch.nunique(),query_unseen_batches=0,selected_genes=len(genes),architecture=SPEC,seed=20260713,max_epochs=200,actual_epochs=int(model.trainer.current_epoch),history_epochs_current_session=len(model.history.get('elbo_train',[])),recovered_from_checkpoint=bool(kw),early_stopping_patience=20,validation_frequency_epochs=5,batch_size=256,internal_split_counts=indices,device='cpu',threads=threads,versions=versions(),weight_hash_before_query=before,weight_hash_after_query=after,query_updates_weights=False,query_method='get_latent_representation give_mean=True; no load_query_data or query training',elapsed_seconds=time.time()-started)
    audit['fingerprints']=fingerprints
    save_json(audit,out/'complete.json');validate_completed_fit(out,m,fingerprints);print(json.dumps(audit,default=str),flush=True)

def _md5(path):
    h = hashlib.md5()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _member(archive, suffix):
    matches = [n for n in archive.namelist()
               if n.endswith(suffix) and not n.startswith("__MACOSX/")]
    if len(matches) != 1:
        raise ValueError(f"Expected one {suffix}: {matches}")
    return matches[0]


def iter_count_blocks(source_dir, cell_manifest, gene_table, source_manifest,
                      download=False):
    """Yield integer cell-by-gene CSR blocks, retaining the archived QC mask.

    Gene selection and model fitting are performed downstream. This function
    recovers source counts and the recorded retained cells; it does not rerun QC.
    """
    source_dir = Path(source_dir)
    obs = pd.read_csv(cell_manifest, sep="\t", index_col="cell_uid")
    genes = pd.read_csv(gene_table, sep="\t", index_col="gene", keep_default_na=False)
    sources = pd.read_csv(source_manifest, sep="\t", keep_default_na=False)
    if not obs.index.is_unique or not genes.index.is_unique:
        raise ValueError("Cell and feature identifiers must be unique")
    covered = set()
    for row in sources.itertuples(index=False):
        path = source_dir / row.archive
        if not path.exists():
            if not download:
                raise FileNotFoundError(f"Missing source archive: {path}")
            source_dir.mkdir(parents=True, exist_ok=True)
            partial = path.with_suffix(path.suffix + ".part")
            with urllib.request.urlopen(row.download_url, timeout=120) as response, partial.open("wb") as output:
                for chunk in iter(lambda: response.read(8 * 1024 * 1024), b""):
                    output.write(chunk)
            if _md5(partial) != row.md5:
                partial.unlink()
                raise ValueError(f"Downloaded source checksum mismatch: {row.archive}")
            partial.replace(path)
        if path.stat().st_size != row.bytes or _md5(path) != row.md5:
            raise ValueError(f"Source checksum/size mismatch: {row.archive}")
        selected = obs.loc[obs.source_file.eq(row.source_file)].copy()
        if len(selected) != row.retained_cells:
            raise ValueError(f"Retained-cell count mismatch: {row.archive}")
        with zipfile.ZipFile(path) as archive:
            raw_genes = pd.read_csv(archive.open(_member(archive, "genes.tsv")),
                                    sep="\t", header=None, names=["gene_id", "gene_symbol"],
                                    keep_default_na=False)
            if not raw_genes.gene_id.equals(genes.gene_id.reset_index(drop=True)):
                raise ValueError(f"Source gene-ID order mismatch: {row.archive}")
            if not raw_genes.gene_symbol.equals(genes.gene_symbol.reset_index(drop=True)):
                raise ValueError(f"Source gene-symbol order mismatch: {row.archive}")
            barcodes = pd.read_csv(archive.open(_member(archive, "barcodes.tsv")),
                                   sep="\t", header=None)[0].astype(str)
            positions = pd.Index(barcodes).get_indexer(selected.original_cell_id)
            if np.any(positions < 0):
                raise ValueError(f"Retained barcodes absent from source: {row.archive}")
            metadata = pd.read_csv(archive.open(_member(archive, "classifications.csv")))
            source_lines = metadata.set_index("barcode").loc[selected.original_cell_id, "singlet_ID"]
            if not np.array_equal(source_lines.astype(str), selected.cell_line.astype(str)):
                raise ValueError(f"Cell-line identity mismatch: {row.archive}")
            matrix = mmread(archive.open(_member(archive, "matrix.mtx"))).tocsc()
            counts = matrix[:, positions].T.tocsr()
            del matrix
        if counts.shape != (len(selected), len(genes)) or np.any(counts.data < 0):
            raise ValueError(f"Count matrix shape/value mismatch: {row.archive}")
        counts.sum_duplicates()
        counts.sort_indices()
        if not np.issubdtype(counts.dtype, np.integer) and not np.all(counts.data == np.floor(counts.data)):
            raise ValueError(f"Noninteger source counts: {row.archive}")
        if counts.data.size and counts.data.max() > np.iinfo(np.int32).max:
            raise ValueError(f"Counts exceed the int32 range: {row.archive}")
        counts = counts.astype(np.int32)
        covered.update(selected.index)
        yield counts, selected, genes
    if covered != set(obs.index):
        raise ValueError("Source manifest does not cover the retained cell manifest")


def load_counts(source_dir, cell_manifest, gene_table, source_manifest,
                download=False):
    """Return AnnData with source integer counts and the original retained-cell order."""
    import anndata as ad

    blocks, annotations = [], []
    for counts, obs, genes in iter_count_blocks(source_dir, cell_manifest, gene_table,
                                               source_manifest, download=download):
        blocks.append(counts)
        annotations.append(obs)
    obs = pd.concat(annotations)
    counts = sparse.vstack(blocks, format="csr")
    del blocks
    expected = pd.read_csv(cell_manifest, sep="\t", usecols=["cell_uid"])["cell_uid"]
    positions = obs.index.get_indexer(expected)
    if not np.array_equal(positions, np.arange(len(obs))):
        counts, obs = counts[positions].copy(), obs.iloc[positions].copy()
    data = ad.AnnData(X=counts, obs=obs.copy(), var=genes.copy())
    data.uns["source_dataset_doi"] = "10.6084/m9.figshare.10298696.v3"
    data.uns["source_license"] = "CC BY 4.0; https://creativecommons.org/licenses/by/4.0/"
    data.uns["counts_status"] = "Original UMI counts; recorded baseline QC mask, without rerunning QC"
    return data


def reconstruct(bundle,source_dir,out,download=False):
    bundle=Path(bundle);out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    a=load_counts(source_dir,bundle/'inputs/mixseq_baseline_cells.tsv.gz',bundle/'inputs/mixseq_baseline_genes.tsv.gz',bundle/'upstream/mixseq_baseline_sources.tsv',download=download)
    assert a.shape==(21936,32738)
    assert a.obs.cell_line.nunique()==169 and a.obs.technical_batch.nunique()==5
    assert np.isfinite(a.X.data).all() and (a.X.data>=0).all() and np.equal(a.X.data,np.floor(a.X.data)).all()
    digest=matrix_hash(a.X)
    assert digest=='518bf5f91c9250b00ad115c0d4c3baccb8324a707665894a3a6320b70c6fe1ea'
    check_runtime();data_fingerprints(a)
    a.write_h5ad(out,compression='gzip')
    save_json(dict(shape=a.shape,nonzero=a.X.nnz,matrix_sha256=digest,all_counts_integer_nonnegative=True,source_dataset_doi='10.6084/m9.figshare.10298696.v3',source_license='CC BY 4.0',retained_cells='Archived QC mask; no new QC fitting',download_requested=download),out.with_suffix('.audit.json'))
    return out


def consolidate(manifest_path,root,input_path):
    check_runtime();fingerprints=data_fingerprints(ad.read_h5ad(input_path),manifest_path)
    root=Path(root);manifest=pd.read_csv(manifest_path,sep='\t',float_precision='round_trip');fits=list(manifest.fit_id.unique());latent=[];audits=[];gene_lists=[];epochs=[]
    for fit in fits:
        folder=root/fit
        audit=validate_completed_fit(folder,manifest[manifest.fit_id.eq(fit)],fingerprints)
        assert audit['weight_hash_before_query']==audit['weight_hash_after_query']
        assert audit['training_batch_count']==5 and audit['query_unseen_batches']==0
        assert audit['seed']==20260713 and audit['selected_genes']==4000 and audit['actual_epochs']<=200
        part=pd.read_csv(folder/'task_latent_means.tsv',sep='\t',float_precision='round_trip')
        want=manifest.loc[manifest.fit_id.eq(fit)&manifest.task_role.isin(['train','validation','test'])].sort_values('row_id')
        assert list(part.row_id)==list(want.row_id) and list(part.pair_id)==list(want.pair_id) and list(part.task_role)==list(want.task_role)
        assert np.isfinite(part[[f'z{i}_mean' for i in range(30)]].to_numpy()).all()
        latent.append(part)
        record={k:v for k,v in audit.items() if not isinstance(v,(dict,list))}
        record.update(audit['fingerprints'])
        membership=pd.read_csv(folder/'cell_membership.tsv.gz',sep='\t',index_col=0)
        training_ids=membership.index[membership.cell_role.eq('reference_fit')]
        record['training_cell_order_sha256']=hashlib.sha256('\n'.join(training_ids).encode()).hexdigest()
        record['software_versions']=json.dumps(audit['versions'],sort_keys=True)
        record['internal_split_counts']=json.dumps(audit['internal_split_counts'],sort_keys=True)
        audits.append(record)
        genes=pd.read_csv(folder/'selected_genes.tsv',sep='\t');genes.insert(0,'fit_id',fit);gene_lists.append(genes)
        for line in (folder/'epoch_log.jsonl').read_text().splitlines():
            e=json.loads(line);epochs.append(dict(fit_id=fit,epoch=e['epoch'],elapsed_seconds=e['elapsed_seconds'],**e['metrics']))
    distinct_states=len({x['weight_hash_before_query'] for x in audits})
    # Symmetric inner folds legitimately use the same reference cells, genes,
    # seed and configuration; independent training can yield identical states.
    expected_sets={tuple(sorted(set(d.loc[d.task_role.ne('train'),'cell_line']))) for _,d in manifest.groupby('fit_id')}
    training_sets=len({x['training_cell_order_sha256'] for x in audits})
    assert training_sets==len(expected_sets)
    state_memberships={}
    for row in audits:state_memberships.setdefault(row['weight_hash_before_query'],set()).add(row['training_cell_order_sha256'])
    assert all(len(group)==1 for group in state_memberships.values()), 'Identical final encoder state across different reference-cell memberships'
    pd.concat(latent,ignore_index=True).to_csv(root/'nested_scvi_task_latents.tsv.gz',sep='\t',index=False,float_format='%.17g')
    pd.DataFrame(audits).to_csv(root/'nested_scvi_fit_audit.tsv',sep='\t',index=False)
    pd.concat(gene_lists,ignore_index=True).to_csv(root/'nested_scvi_selected_genes.tsv.gz',sep='\t',index=False)
    pd.DataFrame(epochs).to_csv(root/'nested_scvi_epoch_metrics.tsv.gz',sep='\t',index=False)
    save_json(dict(fits=len(fits),task_latent_rows=sum(len(x) for x in latent),complete=True,distinct_final_encoder_states=distinct_states,distinct_training_cell_sets=training_sets,expected_distinct_training_line_sets=len(expected_sets),identical_state_only_with_identical_training_membership=True,all_frozen_query_weight_hashes_identical=True,all_training_batches_present=True,all_scvi_fits_fresh=True,scope='Five outer cell-line folds with four independently refit inner representations each. Fixed original QC cell mask; auxiliary nonbenchmark lines retained. Other holdout axes and prospective task construction are outside this sensitivity.'),root/'nested_scvi_verification.json')


def batch(input_path,manifest_path,out,threads=3,parallel=2):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    check_runtime();fingerprints=data_fingerprints(ad.read_h5ad(input_path),manifest_path);manifest=pd.read_csv(manifest_path,sep='\t');queue=[]
    for fit in manifest.fit_id.unique():
        if (out/fit/'complete.json').exists():validate_completed_fit(out/fit,manifest[manifest.fit_id.eq(fit)],fingerprints)
        else:queue.append(fit)
    active={};failed=[];done=[]
    env=os.environ.copy()
    for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:env[v]=str(threads)
    while queue or active:
        while queue and len(active)<parallel:
            fit=queue.pop(0);folder=out/fit;folder.mkdir(parents=True,exist_ok=True);log=(folder/'run.log').open('a')
            cmd=[sys.executable,str(Path(__file__).resolve()),'fit','--input',str(input_path),'--manifest',str(manifest_path),'--fit-id',fit,'--out',str(folder),'--threads',str(threads)]
            job=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
            (folder/'pid').write_text(str(job.pid));active[fit]=(job,log,time.time());print(json.dumps({'event':'launch','fit_id':fit,'pid':job.pid}),flush=True)
        for fit,(job,log,started) in list(active.items()):
            code=job.poll()
            if code is None:continue
            log.close();del active[fit];ok=code==0 and (out/fit/'complete.json').exists()
            (done if ok else failed).append(fit);print(json.dumps({'event':'complete' if ok else 'failed','fit_id':fit,'exit_code':code,'elapsed_seconds':time.time()-started}),flush=True)
            if not ok:queue=[]
        save_json(dict(active=list(active),remaining=queue,completed_this_run=done,failed=failed,updated_unix=time.time()),out/'batch_status.json')
        if active:time.sleep(5)
    if failed:raise RuntimeError('Fit failed: '+', '.join(failed))
    consolidate(manifest_path,out,input_path)


def main():
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='command',required=True)
    x=sp.add_parser('reconstruct');x.add_argument('--bundle',default='.');x.add_argument('--source-cache',required=True);x.add_argument('--out',required=True);x.add_argument('--download',action='store_true')
    x=sp.add_parser('fit');x.add_argument('--input',required=True);x.add_argument('--manifest',required=True);x.add_argument('--fit-id',required=True);x.add_argument('--out',required=True);x.add_argument('--threads',type=int,default=3)
    x=sp.add_parser('batch');x.add_argument('--input',required=True);x.add_argument('--manifest',required=True);x.add_argument('--out',required=True);x.add_argument('--threads',type=int,default=3);x.add_argument('--parallel',type=int,default=2)
    x=sp.add_parser('consolidate');x.add_argument('--manifest',required=True);x.add_argument('--out',required=True);x.add_argument('--input',required=True)
    x=sp.add_parser('all');x.add_argument('--bundle',default='.');x.add_argument('--source-cache',default='results/source_cache/mixseq');x.add_argument('--counts-cache',default='results/source_cache/mixseq_baseline_counts.h5ad');x.add_argument('--manifest',default='inputs/mixseq_nested_folds.tsv');x.add_argument('--out',default='results/nested_scvi/fits');x.add_argument('--download',action='store_true');x.add_argument('--threads',type=int,default=3);x.add_argument('--parallel',type=int,default=2)
    a=p.parse_args()
    if a.command=='reconstruct':reconstruct(a.bundle,a.source_cache,a.out,a.download)
    elif a.command=='fit':run_fit(a.input,a.manifest,a.fit_id,a.out,a.threads)
    elif a.command=='batch':batch(a.input,a.manifest,a.out,a.threads,a.parallel)
    elif a.command=='consolidate':consolidate(a.manifest,a.out,a.input)
    else:
        check_runtime();b=Path(a.bundle).resolve();counts=b/a.counts_cache
        if not counts.exists():reconstruct(b,b/a.source_cache,counts,a.download)
        data_fingerprints(ad.read_h5ad(counts),b/a.manifest)
        batch(counts,b/a.manifest,b/a.out,a.threads,a.parallel)
if __name__=='__main__':main()
