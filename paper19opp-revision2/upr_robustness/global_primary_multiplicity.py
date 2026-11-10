from pathlib import Path
import csv,hashlib,json,datetime,platform
root=Path(__file__).resolve().parent
inputs=[root/'ranking/joint_Holm7_tests.tsv',root/'specificity/primary80_rank_associations.tsv']
rows=[]
for p in inputs:
 with p.open() as f:
  data=list(csv.DictReader(f,delimiter='\t'))
 for d in data:
  if p.parent.name=='ranking':
   rows.append(dict(family='quintile_and_ranking7',hypothesis=d['endpoint'],cohort=d['cohort'],axis='Hallmark90',p=float(d['two_sided_P']),declared_family_q=float(d['Holm7_P'])))
  else:
   rows.append(dict(family='curated_signature80',hypothesis='signed_Spearman',cohort=d['comparison'],axis=d['axis'],p=float(d['raw_permutation_p']),declared_family_q=float(d['BH80'])))
assert len(rows)==87
order=sorted(range(len(rows)),key=lambda i:rows[i]['p'])
m=len(order);last=1.0
for rankpos in range(m-1,-1,-1):
 i=order[rankpos];q=min(1.,rows[i]['p']*m/(rankpos+1));last=min(last,q);rows[i]['global_BH87']=last
out=root/'all87_additional_primary_tests.tsv'
with out.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
manifest={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':platform.python_version(),'inputs':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'tests':87,'scope':'Additional globalBH sensitivity; does not replace prespecified Holm7 or BH80, and does not convert reused-outcome analyses to independent validation.'}
(root/'global_primary_multiplicity_manifest.json').write_text(json.dumps(manifest,indent=2))
print('Completed globalBH87 sensitivity; all original family corrections retained.')
for d in rows:
 if d['family']=='quintile_and_ranking7' or (d['cohort']=='RPE1' and d['axis'] in ['Reactome_UPR_full','Reactome_UPR_no_Hallmark','GO_ER_UPR_full','GO_ER_UPR_no_Hallmark']):print(d)
