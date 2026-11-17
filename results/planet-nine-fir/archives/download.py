"""Refetch the public FITS used by decam_difference.py and verify their recorded SHA256."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

root = Path(__file__).resolve().parent
cache = root / 'receipts'
cache.mkdir(exist_ok=True)
for row in json.loads((root/'decam_inputs.json').read_text()):
    path = cache / (row['label']+'.response')
    if path.exists() and hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==row['sha256']:
        continue
    with urlopen(row['url'],timeout=180) as remote:
        data=remote.read()
    if hashlib.sha256(data).hexdigest()!=row['sha256']:
        raise ValueError(f"Unexpected archive bytes for {row['label']}")
    path.write_bytes(data)
