import csv,gzip,json,collections,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[1];w=r/'final_adjustment'
rank=list(csv.DictReader((r/'independent_validation/frozen/all_source_comparator_rankings.tsv').open(),delimiter='\t'))
rank=sorted([x for x in rank if x['method']=='UPR_priority'],key=lambda x:int(x['rank']))[:20]
for cell in ['HCT116','HEK293T']:
 p=w/f'{cell}_filtered_guide_calls_per_cell.csv.gz'
 if not p.exists():continue
 try:
  cnt=collections.Counter();total=0
  with gzip.open(p,'rt') as f:
   for x in csv.DictReader(f):
    if x['pass_guide_filter']=='True':cnt[x['gene_target']]+=1;total+=1
 except EOFError:continue
 out=[{'gene':x['gene'],'rank':x['rank'],'cells':cnt['GARS1' if x['gene']=='GARS' else x['gene']],'metadata_symbol':'GARS1' if x['gene']=='GARS' else x['gene']} for x in rank]
 with (w/f'{cell}_candidate_cell_counts.tsv').open('w') as f:
  wr=csv.DictWriter(f,fieldnames=['gene','rank','cells','metadata_symbol'],delimiter='\t');wr.writeheader();wr.writerows(out)
 print(cell,total,[(x['gene'],x['cells']) for x in out if x['cells']<30],flush=True)
 (w/f'{cell}_metadata_QA.json').write_text(json.dumps({'source':str(p.name),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'pass_guide_cells':total,'cell_counts_are_guide_metadata_not_RNA_QC':True,'minimum_frozen_top20':min(x['cells'] for x in out)},indent=2))
