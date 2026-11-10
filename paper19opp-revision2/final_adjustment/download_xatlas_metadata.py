"""Acquire guide metadata only; never access expression outcomes."""
import urllib.request,json,hashlib
from pathlib import Path
w=Path(__file__).resolve().parent
url='https://api.figshare.com/v2/articles/29190726'
data=urllib.request.urlopen(url).read();(w/'XAtlas_figshare.json').write_bytes(data)
for entry in json.loads(data)['files']:
 if entry['name'].endswith(('.csv','.csv.gz')):
  p=w/entry['name']
  if not p.exists() or p.stat().st_size!=entry['size']:
   urllib.request.urlretrieve(entry['download_url'],p)
  assert p.stat().st_size==entry['size']
  if entry.get('computed_md5'):
   assert hashlib.md5(p.read_bytes()).hexdigest()==entry['computed_md5']
  print(entry['name'],p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest(),flush=True)
