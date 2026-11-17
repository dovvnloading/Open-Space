"""Test the October source against known satellites at each measured centroid epoch."""
import argparse
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.optimize import minimize_scalar
from geometry import ARCSEC,AU,direction,tangent_basis,pole,query


def run(table,cube,cache,out):
    with table.open() as f:rows=list(csv.DictReader(f))
    xcol,ycol=('x_sky','y_sky') if 'x_sky' in rows[0] else ('x_measured','y_measured')
    rows=[r for r in rows if np.isfinite(float(r[xcol])+float(r[ycol]))]
    times=np.array([float(r['mjd_utc']) for r in rows]);jd=times+2400000.5
    xy=np.array([[float(r[xcol]),float(r[ycol])] for r in rows])
    labels=np.array([r['frame'] if 'frame' in r else r['exposure_integration'].rsplit('_',1)[0] for r in rows])
    groups=np.unique(labels)
    with np.load(cube) as z:scale=float(z['scale_arcsec']);center=z['center_xy']
    planet=query(799,'JWST',jd,cache)
    l,e,n=tangent_basis(planet[:,0],planet[:,1]);dist=planet[:,4,None]*AU
    h=pole('Uranus',np.mean(jd))
    result=[]
    for body,name in [(725,'Perdita'),(714,'Belinda'),(727,'Cupid'),(713,'Rosalind'),(715,'Puck')]:
        sat=query(body,'JWST',jd,cache)
        d=direction(sat[:,0],sat[:,1])*sat[:,4,None]*AU-l*dist
        # Redistribute the known planetocentric XYZ only in orbital longitude.
        # This is a phenomenological phase shift, not a JPL uncertainty model.
        par=(d@h)[:,None]*h; perp=d-par; cross=np.cross(h,perp)
        def positions(theta):
            t=np.deg2rad(theta)
            v=perp*np.cos(t)+cross*np.sin(t)+par
            denom=dist[:,0]+(v*l).sum(axis=1)
            return center+np.column_stack((-(v*e).sum(axis=1),(v*n).sum(axis=1)))*ARCSEC/(scale*denom[:,None])
        def metrics(theta):
            err=positions(theta)-xy
            return np.hypot(err[:,0],err[:,1])
        def fit(idx):
            fn=lambda t:np.mean(np.minimum(metrics(t)[idx],1000.)**2)
            grid=np.arange(-180,180,2)
            v=grid[np.argmin([fn(t) for t in grid])]
            return minimize_scalar(fn,bounds=(v-2.,v+2.),method='bounded').x
        set1=np.arange(len(rows)) < len(rows)//2
        a=fit(np.ones(len(rows),bool));a1=fit(set1);a2=fit(~set1)
        native=metrics(0.);shifted=metrics(a)
        # Delete one complete exposure at a time. This diagnoses repeatability
        # under correlated integrations rather than assuming independent centroids.
        loo=np.array([fit(labels!=g) for g in groups])
        jackknife_se=np.sqrt((len(loo)-1)/len(loo)*np.sum((loo-loo.mean())**2))
        result.append(dict(moon=name,id=body,epochs=len(rows),nominal_median_px=np.median(native),nominal_min_px=native.min(),
                           fitted_phase_offset_deg=a,phase_first_half_deg=a1,phase_last_half_deg=a2,
                           n_exposure_groups=len(loo),leave_one_exposure_min_deg=loo.min(),
                           leave_one_exposure_max_deg=loo.max(),cluster_jackknife_se_deg=jackknife_se,
                           shifted_rms_px=np.sqrt(np.mean(shifted**2)),shifted_median_px=np.median(shifted),
                           last_half_rms_using_first_half_phase_px=np.sqrt(np.mean(metrics(a1)[~set1]**2))))
    with (out/'identity_phase_tests.csv').open('w') as f:
        w=csv.DictWriter(f,list(result[0]),lineterminator='\n');w.writeheader();w.writerows(result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['astrometry','cube','cache','out']:p.add_argument(name,type=Path)
    a=p.parse_args();print(json.dumps(run(a.astrometry,a.cube,a.cache,a.out),indent=2))
