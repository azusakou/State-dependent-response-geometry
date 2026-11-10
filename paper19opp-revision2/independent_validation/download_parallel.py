import urllib.request,os,concurrent.futures,json,hashlib,time
from pathlib import Path
w=Path(__file__).resolve().parent;target=w/'HepG2_complete.h5ad';size=5614460941;chunk=32*1024**2
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad'

(w/'DOWNLOAD_COMPLETE').unlink(missing_ok=True)
fd=os.open(target,os.O_CREAT|os.O_RDWR);os.ftruncate(fd,size);d=w/'download_chunks';d.mkdir(exist_ok=True)
def fetch(i):
 a=i*chunk;b=min(size,(i+1)*chunk)-1;p=d/str(i)
 if p.exists():
  rec=json.loads(p.read_text())
  if hashlib.sha256(os.pread(fd,b-a+1,a)).hexdigest()==rec['sha256']:return
 for retry in range(10):
  try:
   req=urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}'})
   with urllib.request.urlopen(req,timeout=120) as r:
    assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/'),str(r.headers)
    pos=a;h=hashlib.sha256()
    while True:
     buf=r.read(1024**2)
     if not buf:break
     os.pwrite(fd,buf,pos);pos+=len(buf);h.update(buf)
    assert pos==b+1
   p.write_text(json.dumps({'start':a,'end':b,'sha256':h.hexdigest()}));print(i,'done',flush=True);return
  except Exception as e:
   print(i,'retry',retry,str(e),flush=True);time.sleep(min(300,30*(retry+1)))
 raise RuntimeError(i)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(fetch,range(int(os.environ.get('UPR_START_CHUNK','0')),(size+chunk-1)//chunk)))
os.fsync(fd);os.close(fd);
if all((d/str(i)).exists() for i in range((size+chunk-1)//chunk)):(w/'DOWNLOAD_COMPLETE').write_text(str(size));print('COMPLETE',flush=True)
