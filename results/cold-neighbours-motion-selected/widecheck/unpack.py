"""Recreate the exact merged W2-neighbour table from review-sized frozen CSV chunks."""
import gzip
import hashlib
import json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/"observations.json").read_text())
pieces=[]
for i,p in enumerate(manifest["parts"]):
    data=(root/p["path"]).read_bytes()
    if hashlib.sha256(data).hexdigest()!=p["sha256"]:raise ValueError(p["path"])
    rows=gzip.decompress(data)
    pieces.append(rows if i==0 else rows.partition(b"\n")[2])
merged=b"".join(pieces)
if hashlib.sha256(merged).hexdigest()!=manifest["sha256"]:raise ValueError("merged digest mismatch")
(root/"detections.csv").write_bytes(merged)
