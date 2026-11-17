"""Post-unblind adapter: run the untouched frozen pipeline with supplied paths."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(HERE/'frozen'))
import development
import conditional_repeat
import centroid_rule
from astropy.coordinates import SkyCoord
from astropy import units as u

def compare(a,b,path='result'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise ValueError('key mismatch: '+path)
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise ValueError('list mismatch: '+path)
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        if not isinstance(b,(int,float)) or not math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-9):raise ValueError('numeric drift: '+path)
    elif a!=b:raise ValueError('value mismatch: '+path)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--holdout',type=Path,required=True);ap.add_argument('--caldb',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise FileExistsError('output must be a new file')
    development.CAL=args.caldb;conditional_repeat.CAL=args.caldb;conditional_repeat.ROOT=args.baseline
    base_receipt=json.loads((HERE/'evidence'/'baseline-receipt.json').read_text(encoding='utf-8'))
    manifest=json.loads((HERE/'frozen'/'baseline_inputs.json').read_text(encoding='utf-8'))
    expected_files={x['relative']:x['url'] for x in manifest['files']}
    if len(base_receipt)!=len(expected_files) or {x['relative']:x['url'] for x in base_receipt}!=expected_files:raise ValueError('baseline receipt/manifest mismatch')
    for x in base_receipt:
        path=args.baseline/x['relative'];sha=hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):sha.update(block)
        if path.stat().st_size!=x['bytes'] or sha.hexdigest()!=x['sha256']:raise ValueError('baseline receipt mismatch: '+x['relative'])
    pos=centroid_rule.baseline(root=args.baseline)
    center=SkyCoord(float(pos['ra_deg'])*u.deg,float(pos['dec_deg'])*u.deg)
    baseline_actual={'centroid':{'ra_deg':float(pos['ra_deg']),'dec_deg':float(pos['dec_deg']),'offset_arcsec':pos['offset_arcsec']},'nominal':development.primary(development.CENTER,args.baseline),'xy_centroid':development.primary(center,args.baseline)}
    baseline_expected=json.loads((HERE/'frozen'/'frozen_baseline.json').read_text(encoding='utf-8'))
    compare(baseline_expected,baseline_actual,'baseline')
    receipt=json.loads((HERE/'evidence'/'receipt.json').read_text(encoding='utf-8'))
    actual=conditional_repeat.run(args.holdout,receipt)
    expected=json.loads((HERE/'evidence'/'comparison.json').read_text(encoding='utf-8'))
    compare(expected,actual)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(actual,indent=2)+'\n',encoding='utf-8')
    print('Raw replay matched deposited result within 1e-9; wrote',args.output)
