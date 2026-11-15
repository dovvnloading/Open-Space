"""Reproduce the retrospective eclipse geometry from archived state vectors.

Usage: python replay.py [states.csv.gz]
Requires numpy and scipy. All vectors are barycentric geometric ICRF states.
"""
import csv
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

DAY = 86400.0
C = 299792.458
RT, RN, RS = 1352.6, 170.0, 695700.0
MODELS = {
    'Horizons': ('horizons','horizons'),
    'Horizons Triton / NEP105 Nereid': ('horizons','nep105'),
    'NEP097': ('nep097','nep097'),
    'SAI-2015 Triton / SAI-2011 Nereid': ('sai801','sai802'),
}


def solar_fraction(separation, solar, triton):
    """Fraction of an unobscured, uniformly bright solar disk blocked."""
    d = np.maximum(np.asarray(separation),1e-12)
    a = np.clip((d*d+solar**2-triton**2)/(2*d*solar),-1,1)
    b = np.clip((d*d+triton**2-solar**2)/(2*d*triton),-1,1)
    area = solar**2*np.arccos(a)+triton**2*np.arccos(b)
    area -= 0.5*np.sqrt(np.maximum(0.,(-d+solar+triton)*(d+solar-triton)*
                                           (d-solar+triton)*(d+solar+triton)))
    return np.clip(area/(np.pi*solar**2),0,1)


def load_tables(path):
    with gzip.open(path,'rt') as f:
        by_key = {}
        for r in csv.reader(f):
            if r[0]=='event':
                continue
            by_key.setdefault(tuple(r[:3]),[]).append(tuple(map(float,r[3:])))
    return {k:np.array(v) for k,v in by_key.items()}


def measure(tables, year, model):
    triton, nereid = MODELS[model]
    bodies = {}
    for body, solution in [('801',triton),('802',nereid),('399','horizons'),
                           ('10','horizons'),('899','horizons')]:
        data=tables[year,solution,body]
        bodies[body]=CubicHermiteSpline(data[:,0],data[:,1:4],data[:,4:]*DAY,
                                       extrapolate=False)
    origin=tables[year,'horizons','801'][0,0]
    def received(at,source,t):
        receiver=bodies[at](t)
        emission=t
        for _ in range(4):
            delta=bodies[source](emission)-receiver
            emission=t-np.linalg.norm(delta)/(C*DAY)
        return bodies[source](emission)-receiver
    def geometry(offset):
        # Sun crosses about 4.2 light-hours from the barycentre to Nereid.
        # For the local, short state table only the solar state needs extrapolation.
        t=origin+offset
        a=received('802','801',t)
        pos=bodies['802'](t)
        sun=bodies['10'](t)
        step=np.linalg.norm(sun-pos)/(C*DAY)
        # The Sun's barycentric acceleration changes its apparent displacement
        # by far below 1 km in 4.2h, and is further scaled by <0.001 here.
        solar=bodies['10'](t)-bodies['10'](t,1)*step-pos
        n=solar/np.linalg.norm(solar)
        length=np.dot(a,n)
        assert length>0
        d=np.linalg.norm(np.cross(a,n))
        spread=RS*length/np.linalg.norm(solar)
        return d,spread
    def earth_time(h):
        t=origin+h
        pos=bodies['802'](t)
        recv=t+0.175
        for _ in range(4):
            recv=t+np.linalg.norm(bodies['399'](recv)-pos)/(C*DAY)
        return recv
    def utc(t):
        return datetime.fromtimestamp(round((t-2440587.5)*DAY-69.184),tz=timezone.utc).isoformat(timespec='seconds')
    # Center the bounded search on the local evening, retaining both limbs.
    h=minimize_scalar(lambda z:geometry(z)[0],bounds=(0.015,0.14),
                      method='bounded',options={'xatol':1e-12}).x
    d,s=geometry(h)
    def contact(z):
        dist,spread=geometry(z)
        return dist-RT-RN-spread
    ingress=egress=None
    if contact(h)<0:
        ingress=utc(earth_time(brentq(contact,0.001,h)))
        egress=utc(earth_time(brentq(contact,h,.19)))
    rec=earth_time(h)
    earth=bodies['399'](rec)
    def angle(a,b):
        return np.degrees(np.arctan2(np.linalg.norm(np.cross(a,b)),np.dot(a,b)))
    seen=received('399','802',rec)
    # Check geometry with Neptune as possible solar occulter.
    nep=received('802','899',origin+h)
    to_sun=bodies['10'](origin+h)-bodies['802'](origin+h)
    nep_miss=np.linalg.norm(np.cross(nep,to_sun/np.linalg.norm(to_sun)))
    # Uniform projected Nereid disk, uniform solar intensity. This is a
    # sensitivity model, not a measured optical flux or a surface law.
    nodes,wt=np.polynomial.legendre.leggauss(32)
    r=np.sqrt((nodes+1)/2)*RN
    phi=np.arange(192)*2*np.pi/192
    local=np.sqrt(d*d+r[:,None]**2+2*d*r[:,None]*np.cos(phi))
    depth=np.dot(wt/2,solar_fraction(local,s,RT).mean(axis=1))
    return dict(year=year,model=model,peak_utc=utc(rec),
                impact_km=round(float(d),2),solar_radius_km=round(float(s),2),
                ingress_utc=ingress,egress_utc=egress,
                uniform_disk_flux_loss=round(float(depth),4),
                center_solar_coverage=round(float(solar_fraction(d,s,RT)),4),
                neptune_shadow_axis_miss_km=round(float(nep_miss)),
                elongation_deg=round(angle(seen,bodies['10'](rec)-earth),3),
                neptune_sep_arcsec=round(angle(seen,received('399','899',rec))*3600,3),
                triton_sep_arcsec=round(angle(seen,received('399','801',rec))*3600,3))


if __name__=='__main__':
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('states.csv.gz')
    tables=load_tables(path)
    results=[measure(tables,year,model) for year in ('2024','2026') for model in MODELS]
    print(json.dumps(results,indent=2,allow_nan=False))
