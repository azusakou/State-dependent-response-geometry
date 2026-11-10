"""Verify the published compact revision snapshot without downloading data."""
from pathlib import Path
import hashlib
r=Path(__file__).resolve().parents[1]/"paper19opp-revision2"
n=0
for line in (r/"GITHUB_SNAPSHOT.sha256").read_text().splitlines():
    digest,name=line.split("  ",1)
    p=r/name
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    if h.hexdigest()!=digest:
        raise SystemExit("Checksum mismatch: "+name)
    n+=1
print(f"Verified {n} revision files; all SHA256 checks passed.")
