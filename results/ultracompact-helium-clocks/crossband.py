"""Compare independent ZTF g exposures against simultaneous TESS modulation.

Use only groups of at least two g-band exposures in one UTC night. Demean
within each night and use an exact within-night permutation as a sensitivity
test. The TESS phase is fitted without any ZTF measurements.
"""
import argparse
from itertools import permutations, product
import json
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers
from scipy.ndimage import median_filter
from measure import load_cube, FREQ


def compare(cube_path,alert_path):
    t,cube,_=load_cube(cube_path)
    aperture=cube[:,2:5,2:6].sum(axis=(1,2))-12*np.median(cube[:,:,:2],axis=(1,2))
    resid=aperture-median_filter(aperture,size=81)
    with gzip.open(alert_path,'rt') as f: d=pd.DataFrame(json.load(f))
    d=d[(d.mjd>=60290)&(d.mjd<60309)&(d.magpsf<16)&(d.isdiffpos==1)&(d.fid==1)].copy()
    d=d.sort_values('sigmapsf').drop_duplicates('pid').sort_values('mjd')
    d['night']=np.floor(d.mjd).astype(int)
    d=d.groupby('night').filter(lambda g:len(g)>=2).copy()
    iers.conf.auto_download=False
    site=EarthLocation.from_geodetic(-116.8648*u.deg,33.3563*u.deg,1706*u.m)
    utc=Time(d.mjd.values,format='mjd',scale='utc',location=site)
    btjd=(utc.tdb+utc.light_travel_time(SkyCoord(89.6895*u.deg,39.2592*u.deg))).jd-2457000
    predictions=[]
    for epoch in btjd:
        near=np.abs(t-epoch)<.23
        tau=t[near]-epoch
        if np.min(abs(tau))>.0105 or len(tau)<60:raise ValueError("No contemporaneous TESS coverage")
        m=np.column_stack([np.ones(len(tau)),tau,np.sin(2*np.pi*FREQ*tau),np.cos(2*np.pi*FREQ*tau)])
        predictions.append(np.linalg.lstsq(m,resid[near],rcond=None)[0][-1])
    d['tess_e_s']=predictions
    groups=[g for _,g in d.groupby('night')]
    xs=[g.tess_e_s.to_numpy()-g.tess_e_s.mean() for g in groups]
    ys=[g.magpsf.mean()-g.magpsf.to_numpy() for g in groups]
    x,y=np.concatenate(xs),np.concatenate(ys)
    corr=float(np.corrcoef(x,y)[0,1]); alternatives=[]
    for perm in product(*(list(permutations(z)) for z in ys)):
        alternatives.append(np.corrcoef(x,np.concatenate(perm))[0,1])
    tau=np.concatenate([g.mjd-g.mjd.mean() for g in groups])
    def detrend(z):
        return z-tau*np.dot(z,tau)/np.dot(tau,tau)
    return dict(correlation_after_common_linear_trend=float(np.corrcoef(detrend(x),detrend(y))[0,1]),
                utc_nights=[int(g.night.iloc[0]) for g in groups],exposures=len(d),
                correlation=corr,slope_brightness_mag_per_e_s=float(np.dot(x,y)/np.dot(x,x)),
                permutations=len(alternatives),one_sided_permutation_p=float(np.mean(np.asarray(alternatives)>=corr-1e-12)),
                exposures_used=d[['mjd','pid','night','magpsf','sigmapsf','tess_e_s']].to_dict('records'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fits');p.add_argument('alerts')
    p.add_argument('--output',default='crossband.json');args=p.parse_args()
    Path(args.output).write_text(json.dumps(compare(args.fits,args.alerts),indent=2)+'\n')
