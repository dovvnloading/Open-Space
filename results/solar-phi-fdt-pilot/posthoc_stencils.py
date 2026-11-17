"""Post-holdout compaction of frozen g42 bilinear stencils, not preblind code."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

HERE=Path(__file__).resolve().parent

def run(raw_root,out):
    if out.exists():raise FileExistsError('output must be new')
    pred=json.loads((HERE/'evidence'/'prediction.json').read_text(encoding='utf-8'))
    dev=json.loads((HERE/'evidence'/'development.json').read_text(encoding='utf-8'))
    names=dev['development_obs_ids']+dev['holdout_obs_ids'];hashes=dev['development_sha256']+dev['holdout_sha256']
    if pred['patch_id']!='g42' or len(pred['tracked_footprint'])!=4:raise ValueError('unexpected frozen selection')
    pixels=[];fx=[];fy=[];mus=[]
    for name,expected,m in zip(names,hashes,pred['tracked_footprint']):
        path=raw_root/name
        if m['filename']!=name or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('input mismatch: '+name)
        x=np.asarray(m['x']);y=np.asarray(m['y']);mu=np.asarray(m['mu']);ix=np.floor(x).astype(int);iy=np.floor(y).astype(int)
        with fits.open(path,memmap=True) as h:
            image=h[0].data
            st=np.stack((image[iy,ix],image[iy,ix+1],image[iy+1,ix],image[iy+1,ix+1]),axis=1).astype(float)
        if not np.isfinite(st).all() or np.any(np.abs(st)>=4900) or not np.isfinite(mu).all():raise ValueError('invalid science stencil')
        pixels.append(st);fx.append(x-ix);fy.append(y-iy);mus.append(mu)
    out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,stencils=np.array(pixels),fraction_x=np.array(fx),fraction_y=np.array(fy),mu=np.array(mus),obs_ids=np.array(names))
    print('Saved 4 x 317 four-pixel stencils to',out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.raw_root,a.out)
