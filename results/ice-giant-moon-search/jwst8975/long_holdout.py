"""Test six frozen short-exposure tracks in independent long F140M visits 005 and 007."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'jwst2739_uranus'))
from screen import setup,projector,search


def run(cube,ee,short,out):
    rows=list(csv.DictReader(short.open()))
    c,cuts,maps=setup(cube,ee)
    rows=[r for r in rows if r['planet']=='Uranus' and r['filter']=='F140M']
    epochs={float(r['epoch']) for r in rows};epoch,=epochs
    _,project=projector(c,epoch)
    aa=np.array([float(r['radius_km']) for r in rows]);pp=np.array([float(r['phase_deg']) for r in rows])
    x,y=project(aa,pp,slice(None))
    n,w=search.sample_paths(maps,x,y,cuts)
    z=search.score(n,w,c['visit'],cuts)
    out.mkdir(exist_ok=True)
    obs=[];res=[]
    for j,(a,p) in enumerate(zip(aa,pp,strict=True)):
        res.append(dict(a_km=a,phase_deg=p,z_005=z['visit_snr'][j,0],z_007=z['visit_snr'][j,1],n_005=z['n'][j,0],n_007=z['n'][j,1],rank=z['rank'][j]))
        for k,frame in enumerate(c['exposure_id']):
            known=c['known_path'][:,k]
            d=np.hypot(known[:,:,:,None][:,:,0,:]-x[j,k][None,None,:],known[:,:,:,None][:,:,1,:]-y[j,k][None,None,:]).min((1,2))
            idx=np.argmin(d)
            obs.append(dict(a_km=a,phase_deg=p,exposure=frame,x=x[j,k,x.shape[-1]//2],y=y[j,k,x.shape[-1]//2],
                has_support=bool(w[j,k]>0),frame_score=search.snr(n,w)[j,k],
                nearest_known=int(c['known_id'][idx]),nearest_known_trail_px=float(d[idx]),
                perdita_trail_px=float(d[list(c['known_id']).index(725)]),belinda_trail_px=float(d[list(c['known_id']).index(714)])))
    for fn,rr in [('tracks.csv',res),('exposures.csv',obs)]:
        with (out/fn).open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=rr[0]);wr.writeheader();wr.writerows(rr)
    print(res,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','ee','short','out']:p.add_argument(k,type=Path)
    a=p.parse_args();run(a.cube,a.ee,a.short,a.out)
