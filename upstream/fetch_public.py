#!/usr/bin/env python3
"""Fetch the two compact original public files with fixed direct URLs."""
import argparse
from pathlib import Path
from urllib.request import urlopen
import shutil
FILES={
 'replogle_breadth':('https://ndownloader.figshare.com/files/35774443','K562_gwps_raw_bulk_01.h5ad'),
 'multiperturb':('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277747/suppl/GSE277747_RAW.tar','GSE277747_RAW.tar')}
p=argparse.ArgumentParser();p.add_argument('dataset',choices=FILES);p.add_argument('--output-dir',type=Path,default=Path('downloaded_sources'));a=p.parse_args()
u,name=FILES[a.dataset];a.output_dir.mkdir(parents=True,exist_ok=True);target=a.output_dir/name
if target.exists():raise SystemExit('Destination exists; preserve or remove it explicitly before downloading again.')
with urlopen(u,timeout=120) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
print(target)
