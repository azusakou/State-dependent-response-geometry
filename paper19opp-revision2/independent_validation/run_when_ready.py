from pathlib import Path
import time,subprocess,json,hashlib,datetime
w=Path(__file__).resolve().parent;root=w.parent
while not ((w/'DOWNLOAD_COMPLETE').exists() and (w/'HepG2_complete.h5ad').exists()):time.sleep(30)
h=hashlib.sha256()
with (w/'HepG2_complete.h5ad').open('rb') as f:
 for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)

if (w/'complete_input_sha256.json').exists():
 expected=json.loads((w/'complete_input_sha256.json').read_text())['sha256']
 assert h.hexdigest()==expected,'Public matrix differs from pinned original input'
(w/'complete_input_sha256.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':(w/'HepG2_complete.h5ad').stat().st_size,'sha256':h.hexdigest(),'source':'GSE264667; HepG2; official GEO raw matrix'},indent=2))
for stem in ['extract_effects','analyze_holdout','shared_GEM_sensitivity']:
 with (w/(stem+'.stdout.log')).open('wb') as out,(w/(stem+'.stderr.log')).open('wb') as err:
  p=subprocess.run(['Rscript',str(w/(stem+'.R'))],cwd=root,stdout=out,stderr=err)
 if p.returncode:
  (w/'PIPELINE_FAILED.json').write_text(json.dumps({'stage':stem,'exit_code':p.returncode}));raise SystemExit(p.returncode)
(w/'ANALYSIS_COMPLETE').write_text(datetime.datetime.now(datetime.timezone.utc).isoformat())
