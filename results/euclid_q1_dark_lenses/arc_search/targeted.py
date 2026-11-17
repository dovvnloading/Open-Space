"""Apply the frozen VIS ring morphology test to catalogue centres lacking a valid radius."""
import argparse
import csv
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

from search import HERE,find_rings,world
from followup import cut,ring_profile,render

sys.path.insert(0,str(HERE.parent/'imaging'))
from cutouts import ARCHIVES, CutoutArchive, image_extensions


def process(i,archives,where,stamps):
    out={'id_str':i,'archive':where.get(i,''),'status':'no_cutout','sha256':'','eligible':0,'triggers':0,'pass_vis':0,'positions':[]}
    if i not in where:
        return out
    archive=archives[where[i]]
    part=archive.members(i,[i+'.fits'])[i+'.fits']
    out['sha256']=hashlib.sha256(part).hexdigest()
    images=image_extensions(part)
    h,a=images['VIS_FLUX']
    rms=images['VIS_RMS'][1]
    assert abs(abs(float(h['CD1_1']))*3600-.1)<.002
    ny,nx=a.shape[0]//2,a.shape[1]//2
    a=a[:2*ny,:2*nx].astype('f4')
    rms=rms[:2*ny,:2*nx]
    mask=np.isfinite(a)&(a!=0)&np.isfinite(rms)&(rms>0)&(rms<1e8)
    a[~mask]=0
    a=a.reshape(ny,2,nx,2).mean(axis=(1,3))
    mask=mask.reshape(ny,2,nx,2).all(axis=(1,3))
    peaks,out['eligible']=find_rings(a,mask)
    out['status']='examined' if out['eligible'] else 'insufficient_coverage'
    out['triggers']=len(peaks)
    for _,x,y,r,_,_ in peaks:
        binned=cut(a,x,y)
        measure=ring_profile(binned,r)
        ra,dec=world(h,2*x+.5,2*y+.5)
        pos=dict(ra=ra,dec=dec,r=r*.2,**measure)
        out['positions'].append(pos)
        if measure['pass_vis']:
            out['pass_vis']+=1
            render(binned,stamps/f'{i}_{x}_{y}.png')
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,default=HERE/'unmodelled.json')
    p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args()
    rows=list(csv.DictReader((HERE.parent/'catalog/slde_screen.csv').open()))
    targets=[r['id_str'] for r in rows if r['model_state']!='positive_radius']
    archives={key:CutoutArchive(key,args.cache) for key in ARCHIVES}
    where={}
    for name,arc in archives.items():
        for key in arc.entries:
            bits=key.split('/')
            if len(bits)==3 and bits[2]==bits[1]+'.fits':
                if bits[1] in where:
                    # The current release may carry a recenter copy. Prefer the corrected
                    # cutout if one is ever supplied as a successfully modelled system.
                    if where[bits[1]]=='lens':continue
                where[bits[1]]=name
    targetdir=HERE/'stamps'
    targetdir.mkdir(exist_ok=True)
    def run(i):
        try:
            return process(i,archives,where,targetdir)
        except (OSError,ValueError,KeyError,AssertionError) as exc:
            return dict(id_str=i,archive=where.get(i,''),status='error',error=repr(exc))
    with ThreadPoolExecutor(6) as pool:
        result=list(pool.map(run,targets))
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    from collections import Counter
    print(Counter(r['status'] for r in result))
if __name__=='__main__':
    main()
