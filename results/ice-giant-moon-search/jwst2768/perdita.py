"""Test nominal, phase-shifted, and independently predicted Perdita-like trajectories."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import ARCSEC,AU,pole,tangent_basis,generate
from search import Cuts,load_npz,prepare,sample_paths,score,validate_cube


def rotate_tracks(cube,offsets,body=725):
    p=cube['known_path'][list(cube['known_id']).index(body)]
    los,east,north=tangent_basis(cube['path_ra'],cube['path_dec'])
    h=pole(str(cube['planet']),2451545)
    scale=cube['scale_arcsec']/ARCSEC
    delta=cube['path_delta']*AU
    dx=-(p[:,:,0]-cube['center_xy'][0])*scale
    dy= (p[:,:,1]-cube['center_xy'][1])*scale
    # Intersect the observed ray (los + dx*east + dy*north) with
    # the planet's equatorial plane, then rotate that displacement.
    ray=los+dx[:,:,None]*east+dy[:,:,None]*north
    distance=delta*(los@h)/(ray@h)
    r=distance[:,:,None]*ray-delta[:,:,None]*los
    an=np.deg2rad(offsets)[:,None,None,None]
    rot=r*np.cos(an)+np.cross(h,r)*np.sin(an)
    dist=delta+(rot*los).sum(-1)
    x=cube['center_xy'][0]-(rot*east).sum(-1)/dist/scale
    y=cube['center_xy'][1]+(rot*north).sum(-1)/dist/scale
    return x,y


def run(paths,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    angles=np.arange(-180.,180.,.25); cuts=Cuts(psf_support=.8)
    with TemporaryDirectory() as tmp:
        el=Path(tmp)/'orbit.csv'
        el.write_text('a_km,e,i_deg,node_deg,peri_deg,mean_deg\n76337.26,0,0,0,0,254.74603\n')
        for fname in paths:
            c=load_npz(fname)
            if str(c['unmasked_control_body'])!='725':
                raise ValueError('require Perdita-only known-body control cube')
            band=str(c['band']).split('/')[1]
            m=prepare(c,validate_cube(c),cuts)
            x,y=rotate_tracks(c,angles)
            p=c['known_path'][list(c['known_id']).index(725)]
            if np.hypot(x[angles==0]-p[:,:,0],y[angles==0]-p[:,:,1]).max()>1e-4:
                raise ValueError('zero-offset trajectory fails original ephemeris check')
            path=Path(tmp)/'t.npz'
            generate(fname,path,elements=el,model='circular',epoch=2460955.765061979)
            pred=load_npz(path)
            x=np.concatenate([x,pred['x_path']]);y=np.concatenate([y,pred['y_path']])
            n,w=sample_paths(m,x,y,cuts);fit=score(n,w,c['visit'],cuts)
            neg=score(-n,w,c['visit'],cuts)
            visits=np.unique(c['visit'])
            with (out/(band+'_longitude_sweep.csv')).open('w',newline='') as f:
                writer=csv.writer(f);writer.writerow(['model','longitude_offset_deg','rank','pooled_score','pooled_flux_uJy',
                               'negative_rank',*[f'score_{v}' for v in visits],*[f'dropout_{v}' for v in visits],*[f'n_{v}' for v in visits]])
                for i in range(len(x)):
                    writer.writerow(['2025 circular orbit' if i==len(angles) else 'rotated URA184',
                                '' if i==len(angles) else angles[i],fit['rank'][i],fit['snr'][i],fit['flux'][i]*1e6,neg['rank'][i],
                                *fit['visit_snr'][i],*fit['dropout'][i],*fit['n'][i]])
            k=fit['rank'][:-1].argmax(); zero=int(np.flatnonzero(angles==0)[0])
            valid=fit['eligible'][:-1]
            print(band,'precomputed 2025',fit['visit_snr'][-1],fit['n'][-1],'nominal',fit['visit_snr'][zero],fit['n'][zero],
                  'top',angles[k],fit['rank'][k],fit['visit_snr'][k], 'passing',sum(fit['rank'][:-1]>=1),'eligible',valid.sum(),flush=True)
            np.savez_compressed(Path(fname).parent/(band+'-perdita-phase-tracks.npz'),angle=angles,x=x,y=y)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output');p.add_argument('cubes',nargs='+')
    a=p.parse_args();run(a.cubes,a.output)
