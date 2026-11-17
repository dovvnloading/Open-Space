"""Two-epoch circular search over the Perdita–Puck gap in Uranus program 2739."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from astropy.time import Time
import numpy as np

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
sys.path.insert(0, str(root/'jwst8975'))
import geometry
import search
from pilot_inputs import radial_psf
from ring_check import subtract_ring_profile


def setup(path, ee):
    cube = search.load_npz(path)
    band = next(f for f in str(cube['band']).split('/') if f.startswith('F'))
    psf = radial_psf(ee, float(cube['scale_arcsec']), band)
    cube.update(psf=psf, psf_fraction=psf.sum())
    cuts = search.Cuts(psf_support=.80)
    maps = search.prepare(subtract_ring_profile(cube), search.validate_cube(cube), cuts)
    return cube, cuts, maps


def projector(cube, epoch=None):
    delta = cube['path_delta']
    jd = Time(cube['path_mjd'],format='mjd',scale='utc').tdb.jd-delta*geometry.AU/geometry.C/86400
    if epoch is None:
        first = cube['visit']==np.unique(cube['visit'])[0]
        epoch = float(np.median(jd[first]))
    dt = (jd-epoch)*86400
    u,v,_ = geometry.basis('Uranus',epoch)
    geometry_axes = geometry.tangent_basis(cube['path_ra'],cube['path_dec'])
    component = [(axis@u,axis@v) for axis in geometry_axes]
    scale = geometry.ARCSEC/cube['scale_arcsec']
    xref,yref = cube['center_xy']
    def project(a, phase, samples):
        theta = np.deg2rad(phase)[:,None,None]+geometry.frequency('Uranus',a)[:,None,None]*dt[:,samples]
        x = a[:,None,None]*np.cos(theta)
        y = a[:,None,None]*np.sin(theta)
        los,east,north = [x*cu[:,samples]+y*cv[:,samples] for cu,cv in component]
        dist = delta[:,samples]*geometry.AU+los
        return xref-east/dist*scale,yref+north/dist*scale
    return epoch,project


def scan(path, ee, out):
    cube, cuts, maps=setup(path,ee)
    epoch,project=projector(cube)
    middle=cube['path_mjd'].shape[1]//2
    radii=np.arange(79500.,80500.01,.5)
    phases=np.arange(0,360.,.75)
    result=[];eligible=0;screened=[0,0];maximum_spacing=0.
    for r in np.array_split(radii,(len(radii)+15)//16):
        a,ph=[v.ravel() for v in np.meshgrid(r,phases,indexing='ij')]
        x,y=project(a,ph,[middle])
        num,w=search.sample(maps,x[:,:,0],y[:,:,0])
        eligible+=int(search.score(num,w,cube['visit'],cuts)['eligible'].sum())
        for j,sign in enumerate([1,-1]):
            pre=search.score(sign*num,w,cube['visit'],cuts)
            ids=np.flatnonzero(pre['rank']>=.8)
            screened[j]+=len(ids)
            if len(ids)==0:continue
            x,y=project(a[ids],ph[ids],slice(None))
            maximum_spacing=max(maximum_spacing,float(np.hypot(np.diff(x),np.diff(y)).max()))
            n,i=search.sample_paths(maps,x,y,cuts)
            z=search.score(sign*n,i,cube['visit'],cuts)
            for k in np.flatnonzero(z['rank']>=1):
                result.append(dict(sign=sign,a_km=a[ids[k]],phase_deg=ph[ids[k]],rank=z['rank'][k],
                    snr1=z['visit_snr'][k,0],snr2=z['visit_snr'][k,1],n1=z['n'][k,0],n2=z['n'][k,1],
                    flux_uJy=sign*1e6*z['flux'][k]))
    columns=['sign','a_km','phase_deg','rank','snr1','snr2','n1','n2','flux_uJy']
    with out.with_suffix('.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=columns);wr.writeheader();wr.writerows(result)
    rpt=dict(cube_sha256=search.file_hash(path),ee_sha256=search.file_hash(ee),epoch_jd_tdb=epoch,
        band=str(cube['band']),dates_mjd=cube['mjd'].tolist(),grid=[79500,80500,.5,.75],
        templates=int(len(radii)*len(phases)),eligible_midpoint=eligible,midpoint_prefilter_rank=.8,
        positive_and_negative_to_full_trail=screened,positive_and_negative_final=[sum(r['sign']==s for r in result) for s in (1,-1)],
        maximum_spacing_on_scored_full_trails_px=maximum_spacing)
    out.with_suffix('.json').write_text(json.dumps(rpt,indent=2)+'\n')
    print(rpt['templates'],eligible,screened,rpt['positive_and_negative_final'],flush=True)


def control(path,ee,out):
    c,cuts,maps=setup(path,ee)
    body=int(c['unmasked_control_body'])
    pathxy=c['known_path'][np.flatnonzero(c['known_id']==body).item()]
    ofs=[(0,0)]+[(r*np.cos(t),r*np.sin(t)) for r in [7,11,15] for t in np.arange(12)*np.pi/6]
    x=np.array([pathxy[:,:,0]+d for d,_ in ofs]); y=np.array([pathxy[:,:,1]+d for _,d in ofs])
    n,w=search.sample_paths(maps,x,y,cuts)
    z=search.score(n,w,c['visit'],cuts)
    out.write_text(json.dumps(dict(body=body,visit_snr=z['visit_snr'][0].tolist(),visit_counts=z['n'][0].tolist(),
        offsets=ofs,control_ranks=z['rank'].tolist(),flux_uJy=1e6*z['flux'][0],cube_sha256=search.file_hash(path)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube',type=Path);p.add_argument('ee',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--known',action='store_true')
    a=p.parse_args();(control if a.known else scan)(a.cube,a.ee,a.out)
