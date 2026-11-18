"""Restore exact upstream inputs when still served, checking SHA-256 before use.
Usage: python fetch_receipts.py /path/to/cache [receipt_key ...]
"""
import hashlib
import json
import shutil
import sys
import urllib.request
from pathlib import Path

base=Path(__file__).resolve().parent
cache=Path(sys.argv[1]);cache.mkdir(exist_ok=True,parents=True)
wanted=set(sys.argv[2:])
for r in json.loads((base/"receipts.json").read_text()):
    if wanted and r["key"] not in wanted:continue
    path=cache/r["filename"]
    if path.exists() and hashlib.file_digest(path.open("rb"),"sha256").hexdigest()==r["sha256"]:continue
    temp=path.with_suffix(path.suffix+".download")
    try:
        with urllib.request.urlopen(r["url"],timeout=120) as src,temp.open("wb") as dst:
            shutil.copyfileobj(src,dst)
        with temp.open("rb") as dst: digest=hashlib.file_digest(dst,"sha256").hexdigest()
        if digest!=r["sha256"]:raise ValueError(f"{r['key']} upstream snapshot changed")
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)
