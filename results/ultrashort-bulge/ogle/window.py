"""Exact idealized five-exposure window fraction on real EWS alert-star timestamps.

Threshold crossings are computed for a noiseless injected light curve relative to the
per-epoch flux error. This measures cadence/noise opportunity, conditional on alerts.
"""
import argparse
import csv
from datetime import datetime,timezone
from pathlib import Path
import sys

import numpy as np
import VBMicrolensing

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'ogle-dark-lenses'))
import lens

# Uniform source. Uniform event centers, for one fixed impact parameter.
SCENARIOS=[("2h_PSPL",2/24,0.,1.),("6h_PSPL",.25,0.,1.),("12h_PSPL",.5,0.,1.),
           ("6h_FSPL",.25,1.5,1.),("6h_FSPL_blend",.25,1.5,.3)]
U0=.5
R=np.geomspace(.0001,100,3500)


def union_fraction(t,half,start,end):
    # For epochs i..i+4 each peak center must lie inside all five individual
    # detectability intervals simultaneously. Union across all consecutive groups.
    if len(t)<5:return 0.
    l=np.maximum(np.lib.stride_tricks.sliding_window_view(t-half,5).max(axis=1),start)
    r=np.minimum(np.lib.stride_tricks.sliding_window_view(t+half,5).min(axis=1),end)
    good=r>l
    order=np.argsort(l[good]); l,r=l[good][order],r[good][order]
    before=np.r_[start,np.maximum.accumulate(r)[:-1]] if len(r) else np.array([])
    return np.maximum(0,r-np.maximum(l,before)).sum()/(end-start)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);cache=ap.parse_args().cache
    rows=list(csv.DictReader((HERE/'catalog.csv').open()))
    vb=VBMicrolensing.VBMicrolensing()
    vb.LoadESPLTable(str(Path(VBMicrolensing.__file__).parent/'data/ESPL.tbl'))
    curves={0.:lens.magnification(R**2),1.5:np.minimum.accumulate([vb.ESPLMag2(v,1.5) for v in R])}
    out=[]
    for row in rows:
        event=row['event'];year,_,n=event.split('-');p=cache/year/f'blg-{n}.dat'
        if not p.exists():continue
        t,mag,err,_,_=lens.read_phot(p)
        if not len(t):continue
        # HJD - 2450000. Feb 1 to Nov 1 00h UT, identically in both seasons.
        def date(month):
            return datetime(int(year),month,1,tzinfo=timezone.utc).timestamp()/86400+2440587.5-2450000.
        start,end=date(2),date(11)
        in_season=(t>=start-2)&(t<end+2)
        ind=np.argsort(t[in_season]);season_t,season_err=t[in_season][ind],err[in_season][ind]
        # Inflate per-epoch precision if the typical baseline scatters more than the
        # supplied formal errors with the same error floor as the fits.
        flux,sig=lens.to_flux(mag,err)
        base=np.median(flux)
        scale=max(1.,np.median(np.abs((flux-base)/sig))/.67449)
        point=dict(event=event,year=year,field=row['field'].split('.')[0],n_season=len(season_t))
        for label,te,rho,fs in SCENARIOS:
            # F_source(A-1) >= 3 sigma; quoted noise applied as a baseline fractional error.
            min_amp=1+3*scale*np.hypot(season_err,.004)*(.4*np.log(10))/fs
            rp=np.interp(min_amp,curves[rho][::-1],R[::-1],left=R[-1],right=0)
            half=te*np.sqrt(np.maximum(0,rp**2-U0**2))
            point[label]=float(union_fraction(season_t,half,start,end))
        out.append(point)
    if len(out)!=len(rows):raise ValueError(f'only {len(out)} of {len(rows)} light curves available')
    groups={}
    for row in out:
        for k in [(row['year'],row['field']),(row['year'],'ALL')]:
            groups.setdefault(k,[]).append(row)
    fields=['year','field','n_curves','median_n_season']+[v[0] for v in SCENARIOS]
    with (HERE/'window_by_field.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fields,lineterminator='\n');w.writeheader()
        for (year,field),g in sorted(groups.items()):
            w.writerow(dict(year=year,field=field,n_curves=len(g),median_n_season=np.median([x['n_season'] for x in g]),
                            **{label:round(float(np.mean([x[label] for x in g])),6) for label,*_ in SCENARIOS}))

if __name__=='__main__':main()
