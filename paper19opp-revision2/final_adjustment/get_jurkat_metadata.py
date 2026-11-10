from pathlib import Path
import urllib.request,os,json,time,hashlib
w=Path(__file__).resolve().parent;N=9366490264;url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad';fd=os.open(w/'Jurkat_metadata_sparse.h5ad',os.O_RDWR|os.O_CREAT);os.ftruncate(fd,N);records=[]
for a,b in [(0,1024**2-1),(N-64*1024**2,N-32*1024**2-1),(N-32*1024**2,N-1)]:
 for retry in range(6):
  try:
   req=urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}'})
   with urllib.request.urlopen(req,timeout=180) as r:
    assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/');v=r.read();assert len(v)==b-a+1;os.pwrite(fd,v,a)
   records.append({'start':a,'end':b,'sha256':hashlib.sha256(v).hexdigest()});print(a,b,'done',flush=True);break
  except Exception as e:
   print(retry,str(e),flush=True);time.sleep(15*(retry+1))
 else:raise RuntimeError('metadata download failed')
os.close(fd);(w/'Jurkat_metadata_ranges.json').write_text(json.dumps(records,indent=2))
