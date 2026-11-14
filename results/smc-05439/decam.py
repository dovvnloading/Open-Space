"""Test the two OGLE templates on unforced DELVE DR3 Magellanic Clouds sources."""
import json
import gzip
import csv
import numpy as np
from replay import DATA, HERE, records, training, templates, fit


def read(name):
    with gzip.open(DATA/name,'rt') as f:
        return list(csv.DictReader(f))


def hjd(mjd,seconds):
    # Low-order heliocentric correction, from the smc-01293 replay, at this target.
    jd=mjd+2400000.5+seconds/172800
    d=jd-2451545
    anomaly=np.radians((357.528+.9856003*d)%360)
    lon=np.radians((280.460+.9856474*d+1.915*np.sin(anomaly)+.020*np.sin(2*anomaly))%360)
    obl=np.radians(23.439-.0000004*d)
    radius=1.00014-.01671*np.cos(anomaly)-.00014*np.cos(2*anomaly)
    ra,dec=np.radians([17.2197967,-72.3103774])
    proj=np.cos(lon)*np.cos(ra)*np.cos(dec)+np.sin(lon)*(np.cos(obl)*np.sin(ra)*np.cos(dec)+np.sin(obl)*np.sin(dec))
    return jd-radius*proj*499.005/86400-2450000


def usable(name):
    keep={}
    for r in read(name):
        if int(r['forced']) or r['filter'] not in ('g','r','i') or not 0<float(r['err'])<.12 or float(r['mag'])>=40:
            continue
        key=r['id'].partition('_')[0]
        if key not in keep or float(r['err'])<float(keep[key]['err']):
            keep[key]=r
    return keep


def project(source,fs,b,zp=None,seconds=None):
    rows=list(source.items())
    if zp is not None:
        rows=[(k,v) for k,v in rows if k in zp]
    if len(rows)<20:
        return {'usable':len(rows)}
    exp,rr=zip(*rows)
    labels,band=np.unique([r['filter'] for r in rr],return_inverse=True)
    bandcols=np.eye(len(labels))[band]
    y=np.array([float(r['mag']) for r in rr])
    e=np.hypot([float(r['err']) for r in rr],.008)
    if zp is not None:
        y-=np.array([zp[k] for k in exp])
    sec=np.array([seconds.get(k,0) for k in exp]) if seconds else np.zeros(len(rr))
    t=hjd(np.array([float(r['mjd']) for r in rr]),sec)
    x=np.column_stack((bandcols,templates(t,fs,b)))
    coeff,cov=fit(x,y,e,t)
    return dict(usable=len(t),nights=len(np.unique(np.floor(t))),
                filters=dict(zip(labels.tolist(),np.bincount(band).tolist())),
                scale_orbit_even_odd=coeff[-3:].tolist(),se_orbit_even_odd=np.sqrt(np.diag(cov)[-3:]).tolist())


def main():
    fs,b=training(list(records()))
    target=usable('dc-46872.csv.gz')
    sec={}
    for p in DATA.glob('durations*.csv.gz'):
        for r in read(p.name):
            sec[str(int(r['expnum'])).zfill(8)]=float(r['exptime'])
    result={'target':project(target,fs,b,seconds=sec), 'matched_exposure_times':sum(k in sec for k in target)}
    by_exp={}
    for p in sorted(DATA.glob('dc-*.csv.gz')):
        n=p.name.split('-')[1].split('.')[0]
        if n=='46872':
            continue
        rows=usable(p.name)
        if n=='59764':
            result['nearest_neighbor']=project(rows,fs,b,seconds=sec)
            continue
        for band in ('g','r','i'):
            mags=[float(r['mag']) for r in rows.values() if r['filter']==band]
            if len(mags)<3:
                continue
            med=np.median(mags)
            for k,r in rows.items():
                if r['filter']==band:
                    by_exp.setdefault(k,[]).append(float(r['mag'])-med)
    zp={k:np.median(v) for k,v in by_exp.items() if len(v)>=3}
    result['target_differential']=project(target,fs,b,zp,sec)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
