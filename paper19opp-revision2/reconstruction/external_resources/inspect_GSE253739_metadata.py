#!/usr/bin/env python3
"""Inspect sample/support metadata only; no response model or outcome test."""
import csv,gzip,re,collections
from pathlib import Path
out=Path(__file__).resolve().parent
source=Path('/Volumes/wang/paper8ikk/data/raw/GSE253739_FIGSHARE/all_data_final_lineages.RData')
with gzip.open(source,'rb')as f:header=f.read(16384)
names=[x.decode()for x in re.findall(rb'[A-Za-z0-9_][A-Za-z0-9_.:-]{3,80}',header)[:8]]
root=Path('/Volumes/wang/paper8ikk/results/GSE253739/history_conditioned_fate')
rows=list(csv.DictReader(gzip.open(root/'longitudinal_lineage_exposure_pairs.tsv.gz','rt'),delimiter='\t'))
base={}
for r in rows:
 key=(r['lineage'],r['first']);n=int(r['parent_cells'])
 assert key not in base or base[key]==n
 base[key]=n
summary=[]
for k in [1,2,3,5,10]:
 eligible={key for key,n in base.items()if n>=k}
 for depth in [1,2,5]:
  a=[r for r in rows if (r['lineage'],r['first'])in eligible and int(r['desc_cells'])>=depth]
  summary.append(dict(minimum_parent_cells=k,eligible_parent_history_units=len(eligible),eligible_parent_unique_lineages=len({x[0]for x in eligible}),minimum_descendant_cells=depth,observed_paired_tasks=len(a),paired_unique_lineages=len({r['lineage']for r in a})))
with(out/'GSE253739_metadata_support.tsv').open('w')as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
(out/'GSE253739_metadata_assessment.txt').write_text(f'''Purpose: feasibility of independent measured clone-state validation only. No temporal policy, treatment-duration optimization, response model or new biological test was run.

Source RData: {source.stat().st_size:,} compressed bytes. A 16-KiB decompressed-header inspection found RDX3 serialization, an object named all_data and assays/counts slots, consistent with the author-described Seurat object. The full object was not loaded. Raw matrices and cell barcodes for twelve samples are also present locally.

Primary study: Schaff et al. Clonal differences underlie variable responses to sequential and prolonged treatment. Cell Systems (2024), GSE253739. Primary paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11003565/ . The study measures barcoded WM989 melanoma clones following cisplatin, cobalt chloride or dabrafenib/trametinib, then continued or switched treatment. A measured clone state after the first treatment can precede measured relatives after the second treatment; it is not untreated pretreatment state and is not the same individual cell followed longitudinally.

Compact pre-existing metadata are available in lineage_condition_centroids.tsv.gz (776,888 bytes) and longitudinal_lineage_exposure_pairs.tsv.gz under the prior local analysis. The former contains lineage, first, second, n_cells and expression-PC columns; the latter contains lineage, first, second, parent_cells, desc_cells and derived outcome/PC fields. Only lineage/treatment/cell-support metadata were inspected; no derived outcome tests or PC values were interpreted. These counts have not been independently reconstructed from the author RData in this task.

The compact ledger contains 2,561 baseline lineage-by-first-treatment units representing 2,441 unique lineages, crossed with three second treatments to give 7,683 exposure rows. These are correlated repeated exposures, not 7,683 independent biological units. Of these, 727 tasks from 510 lineages have at least one sampled descendant. At >=3 measured parent cells and >=2 measured descendant cells, support is 257 tasks from 166 lineages. At >=10 parent cells and >=5 descendant cells, support is 85 tasks from 52 lineages. The depth restriction materially changes the population.

This resource is usable for a separate clone-grouped assessment of whether directly measured pre-second-treatment clone states predict subsequent transcriptomic response or observed clonal abundance, conditional on first-treatment history, parent abundance and cell support. All branches from the same lineage must remain together across splits and bootstrap draws; states and dimensions must be fit using baseline/training samples only. Absent sampled descendants do not establish clone death. Original culture-level biological replication and complete barcode mapping must be checked from author metadata before causal or generalizable claims. No such analysis has been launched.
''')
print('Metadata feasibility report written; no response effects analyzed.')
