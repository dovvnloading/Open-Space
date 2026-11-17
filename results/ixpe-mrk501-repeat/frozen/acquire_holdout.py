"""Explicit future acquisition command; do not run before root authorization."""
import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--receipt',required=True);a=ap.parse_args()
    root=Path(a.root);receipt=Path(a.receipt)
    if root.exists() or receipt.exists():raise FileExistsError('root and receipt must be new')
    root.mkdir(parents=True);m=json.loads((Path(__file__).parent/'holdout_inputs.json').read_text(encoding='utf-8'));results=[]
    for x in m['files']:
        path=root/x['relative'];path.parent.mkdir(parents=True,exist_ok=True);digest=hashlib.sha256();size=0
        with urllib.request.urlopen(x['url'],timeout=60) as r,path.open('xb') as f:
            while chunk:=r.read(1024*1024):f.write(chunk);digest.update(chunk);size+=len(chunk)
        if size!=x['bytes']:raise ValueError('archive size changed: '+x['relative'])
        results.append({'relative':x['relative'],'url':x['url'],'bytes':size,'sha256':digest.hexdigest()})
    receipt.write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
