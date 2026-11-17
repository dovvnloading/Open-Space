"""Measure the Jan 2024 source's longitude relative to predicted Perdita and Belinda."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import AU,ARCSEC,basis,tangent_basis
from search import load_npz


def angle(positions,cube,ix,j):
    u,v,h=basis('Uranus',2451545.)
    s,e,n=tangent_basis(cube['path_ra'][ix,j],cube['path_dec'][ix,j]);d=cube['path_delta'][ix,j]*AU
    sc=float(cube['scale_arcsec'])/ARCSEC
    dx=-(positions[:,0]-cube['center_xy'][0])*sc;dy=(positions[:,1]-cube['center_xy'][1])*sc
    ray=s+dx[:,None]*e+dy[:,None]*n;rg=d*(s@h)/(ray@h)
    r=rg[:,None]*ray-d[:,None]*s
    return np.rad2deg(np.arctan2(r@v,r@u))


def run(cubes,folder):
    root=Path(folder);res=[]
    for path in cubes:
        c=load_npz(path);band=str(c['band']).split('/')[1]
        for body,suf in [(725,'_native.csv'),(714,'_belinda_native.csv')]:
            rows=list(csv.DictReader((root/(band+suf)).open()))
            ix=np.array([list(c['exposure_id']).index(r['frame'].rsplit('_',2)[0]) for r in rows])
            obs=np.array([[float(r['sky_x']),float(r['sky_y'])] for r in rows])
            ref=c['known_xy'][ix,list(c['known_id']).index(body)];j=c['path_mjd'].shape[1]//2
            ph=(angle(obs,c,ix,j)-angle(ref,c,ix,j)+180)%360-180
            for i,r in enumerate(rows):
                res.append(dict(body=body,band=band,visit=r['visit'],frame=r['frame'],mjd_utc=r['mjd_utc'],phase_delta_deg=ph[i],
                                center_residual_px=np.linalg.norm(obs[i]-ref[i])))
    with (root/'frame_longitudes.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=res[0]);w.writeheader();w.writerows(res)
    stats={}
    for body in [725,714]:
        selected=[r for r in res if r['body']==body];phi=np.array([r['phase_delta_deg'] for r in selected])
        stats[str(body)]=dict(n=len(selected),median_deg=float(np.median(phi)),mad_deg=float(np.median(abs(phi-np.median(phi)))),
            by_visit={v:np.median([r['phase_delta_deg'] for r in selected if r['visit']==v]) for v in sorted({r['visit'] for r in selected})})
    (root/'frame_longitudes.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(stats)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder');p.add_argument('cubes',nargs='+')
    a=p.parse_args();run(a.cubes,a.folder)
