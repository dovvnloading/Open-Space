"""Replay the ground-archive checks. Requires Python 3, NumPy, SciPy and Astropy."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'inputs'
SOURCE = SkyCoord(21.06275*u.deg, -73.3259444444*u.deg)
PERIOD = 1.1696719
EPOCH = 2460000.0369
SITES = {'DECam':(-70.81489,-30.16606,2200), 'SkyMapper':(149.061194,-31.27336,1165),
         'VISTA':(-70.3975,-24.6158,2518)}
iers.conf.auto_download = False


def table(name):
    with (RAW/name).open() as f:
        return list(csv.DictReader(f))


def hjd(site, mjd, duration):
    lon,lat,height=SITES[site]
    obs=EarthLocation.from_geodetic(lon*u.deg,lat*u.deg,height*u.m)
    mid=Time(np.asarray(mjd)+np.asarray(duration)/172800,format='mjd',scale='utc',location=obs)
    return mid.utc.jd + mid.light_travel_time(SOURCE,kind='heliocentric').to_value(u.day)


def wrap(x):
    return (x+.5)%1-.5


def phases(t):
    return ((t-EPOCH)/PERIOD)%1


def archive_epochs():
    exposures={r['expnum']:float(r['exptime']) for r in table('delve-exposures.csv')}
    groups={}
    for site, file, tkey, durkey, mkey, ekey, bandkey in [
        ('DECam','delve-meas.csv','mjd',None,'mag','err','filter'),
        ('DECam','nsc-exposures.csv','mjd','exptime','mag_auto','magerr_auto','filter'),
        ('SkyMapper','skymapper-epochs.csv','date','exp_time','mag_psf','e_mag_psf','filter')]:
        data=table(file)
        duration=[float(r[durkey]) if durkey else exposures[r['id'].split('_')[0]] for r in data]
        times=hjd(site,[float(r[tkey]) for r in data],duration)
        groups[file]=[dict(hjd=float(t),phase=float(phases(t)),band=r[bandkey],mag=float(r[mkey]),
                           err=float(r[ekey]),flags=r.get('flags'),id=r.get('id',r.get('image_id',r.get('measid'))))
                      for r,t in zip(data,times,strict=True)]
    for key in ('ogle3-V.dat','ogle4-V.dat'):
        t,m,e=np.loadtxt(RAW/key).T
        if key.startswith('ogle3'):
            t+=2450000
        groups[key]=[dict(hjd=float(h),phase=float(phases(h)),band='V',mag=float(tmag),err=float(terr))
                     for h,tmag,terr in zip(t,m,e)]
    return groups


def nightly_baseline(rows):
    by_day={}
    for r in rows:
        d=min(abs(wrap(r['phase'])),abs(wrap(r['phase']-.5)))
        if d>.11:
            by_day.setdefault(int(r['hjd']),[]).append(r['mag'])
    return (float(np.median([np.median(m) for m in by_day.values()])),len(by_day)) if by_day else (None,0)


def photometry(rows):
    usable=[r for r in rows if .0<r['err']<=.1 and 0<r['mag']<40]
    with_excess=[]
    bands={}
    for band in sorted({r['band'] for r in usable}):
        filtered=[r for r in usable if r['band']==band]
        base,n=nightly_baseline(filtered)
        bands[band]={'off_eclipse_median':base,'off_eclipse_nights':n}
        if base is None:
            continue
        for r in filtered:
            with_excess.append(dict(**r,excess_mag=r['mag']-base,baseline_nights=n))
    return {'raw':len(rows),'error_cut':len(usable),'first_hjd':min(r['hjd'] for r in rows),
            'last_hjd':max(r['hjd'] for r in rows),'baselines':bands,'epochs':with_excess}


def fit(data):
    # Limit DELVE to the 2017-20 cluster and g/r/i. Template centers share phase across bands.
    # The duration and nominal dip magnitude are imported from the earlier OGLE analysis.
    # The baseline is estimated anew per passband, with each observing night weighted once.
    rows=[r for r in data if 2457754.5<=r['hjd']<2459215.5 and r['band'] in 'gri']
    ph=np.array([r['phase'] for r in rows])
    exc=np.array([r['excess_mag'] for r in rows])
    sig=np.hypot([r['err'] for r in rows],.025)
    def optimize(width,lo_depth,hi_depth):
        def residual(x):
            c1,c2,d1,d2=x
            y=d1*np.exp(-.5*(wrap(ph-c1)/width)**2)+d2*np.exp(-.5*(wrap(ph-c2)/width)**2)
            return (y-exc)/sig
        result=least_squares(residual,[0,.52,.4,.4],bounds=([-.06,.44,lo_depth,lo_depth],[.06,.58,hi_depth,hi_depth]),loss='soft_l1')
        return dict(phase_1=result.x[0],phase_2=result.x[1],spacing=result.x[1]-result.x[0],
                    depth_1=result.x[2],depth_2=result.x[3],width=width)
    return {'exposures':len(rows),'distinct_nights':len({int(r['hjd']) for r in rows}),
            'near_primary_nights':sorted({int(r['hjd']) for r in rows if abs(wrap(r['phase']))<.045}),
            'near_secondary_nights':sorted({int(r['hjd']) for r in rows if abs(wrap(r['phase']-.52))<.045}),
            'fits':[optimize(width,.25,.6) for width in (.017,.022,.028,.035)],
            'warning':'Only one primary-core night, pooled epoch, sparse secondary. Centers depend on assumed width and allowed depth; do not use as a self-standing period or drift measurement.'}


def main():
    inputs={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in RAW.iterdir() if f.is_file()}
    groups=archive_epochs()
    out={'input_hashes_sha256':inputs,'photometry':{k:photometry(v) for k,v in groups.items()}}
    out['DELVE_spacing_diagnostic']=fit(out['photometry']['delve-meas.csv']['epochs'])
    ids=lambda name,key:{int(r[key].split('_')[0]) for r in table(name)}
    nsc={int(r['expnum']) for r in table('nsc-exposures.csv')}
    delve=ids('delve-meas.csv','id')
    smash={int(r['expnum']) for r in table('smash-meas.csv')}
    out['DECam_overlap']={'DELVE':len(delve),'NSC':len(nsc),'SMASH':len(smash),'DELVE_and_NSC':len(delve&nsc),
                          'SMASH_and_DELVE':len(smash&delve),'SMASH_and_NSC':len(smash&nsc)}
    objects=table('delve-cone.csv')
    target=next(r for r in objects if r['objid']=='0210m732.7445')
    neigh=next(r for r in objects if r['objid']=='0210m732.48525')
    out['neighbor']={'id':neigh['objid'],'separation_from_OGLE_arcsec':float(SOURCE.separation(SkyCoord(float(neigh['ra'])*u.deg,float(neigh['dec'])*u.deg)).arcsec),'bands':{}}
    for b in 'gri':
        d=float(neigh[b+'mag'])-float(target[b+'mag'])
        out['neighbor']['bands'][b]={'mag_difference':d,'max_unresolved_depth_if_neighbor_disappears_mag':float(2.5*np.log10(1+10**(-.4*d)))}
    out['VISTA']={}
    for band in ['Y','J','Ks']:
        x=json.loads((RAW/f'vmc-{band}.json').read_text())
        r=[dict(zip([m['name'] for m in x['metadata']],line)) for line in x['data']]
        out['VISTA'][band]={'n':len(r),'date_range_MJD':[min(z['MJD'] for z in r),max(z['MJD'] for z in r)],
                           'zero_quality_bits':sum(int(z[band.upper()+'PPERRBITS'])==0 for z in r)}
    out['late_OGLE_V']=[r for r in groups['ogle4-V.dat'] if r['hjd']>=2459700.]
    (ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Spacing diagnostic',out['DELVE_spacing_diagnostic'])
    print('Overlap',out['DECam_overlap'],'Neighbor',out['neighbor'],'VISTA',out['VISTA'])

if __name__=='__main__':
    main()
