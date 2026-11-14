"""Resolve the daily alias of a 4.25-hour brightness signal toward SMC-ECL-01947."""
from pathlib import Path
import csv
import gzip
import json

import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

DATA = Path(__file__).resolve().parent / 'data'
PERIOD = 2.4013454
EPOCH = 10000.0436  # HJD minus 2450000. A common numerical reference.
HARMONICS = 11


def orbit(t, order):
    angle = (t - EPOCH)*2*np.pi/PERIOD
    return np.column_stack([op(k*angle) for k in range(1, order+1) for op in (np.cos,np.sin)])


def fit(t, y, err, freq=None):
    x = np.column_stack((np.ones(len(t)), orbit(t,HARMONICS)))
    if freq is not None:
        phi = 2*np.pi*(t-EPOCH)*freq
        x = np.column_stack((x,np.cos(phi),np.sin(phi)))
    w = 1/err
    for step in range(4):
        coef = np.linalg.lstsq(x*w[:,None], y*w, rcond=None)[0]
        residual = y-x@coef
        scatter = max(.003, 1.483*np.median(abs(residual-np.median(residual))))
        if step < 3:
            w = np.sqrt(np.minimum(1.,3*scatter/np.maximum(abs(residual),1e-12)))/err
    bread = np.linalg.pinv((x*w[:,None]).T@(x*w[:,None]))
    ng,groups = np.unique(np.floor(t),return_inverse=True)
    scores = np.zeros((len(ng),x.shape[1]))
    np.add.at(scores,groups,x*(residual*w*w)[:,None])
    cov = bread @ scores.T @ scores @ bread
    z = abs(residual/scatter)
    loss = np.sum(np.where(z<3,z*z,6*z-9)*scatter**2/err**2)
    return coef, cov, loss, residual


def series():
    for stage in [4,3]:
        for band in 'IV':
            t,y,e = np.loadtxt(DATA/f'ogle{stage}-{band}.dat.gz').T
            ok = np.isfinite(t+y+e) & (e>0) & (e<min(.13,max(.03,3*np.median(e))))
            t,y,e = t[ok],y[ok],e[ok]
            yield f'OGLE-{stage}-{band}',t-(2450000 if stage==4 else 0),y-np.median(y),np.maximum(e,.7*np.median(e))
    with gzip.open(DATA/'xtime.tsv.gz','rt') as source:
        lines = [line for line in source if line.strip() and not line.startswith('#')]
    rows = list(csv.DictReader(lines[:1]+lines[3:], delimiter='\t'))
    for obj in sorted({r['EROS'].strip() for r in rows}):
        for band in 'RB':
            seq = []
            for r in rows:
                if r['EROS'].strip()!=obj or not r[band+'mag'].strip() or not r['e_'+band+'mag'].strip():
                    continue
                t,y,e = map(float,(r['HJD'],r[band+'mag'],r['e_'+band+'mag']))
                if np.isfinite(t+y+e) and 0 < e < .2 and y < 40:
                    seq.append((t,y,e))
            if len(seq)>30:
                t,y,e = np.array(seq).T
                yield f'EROS-{obj}-{band}',t,y-np.median(y),np.maximum(e,.6*np.median(e))


def tess(path, freq):
    t,y = np.loadtxt(path,delimiter=',',skiprows=1).T
    t = t+7000-.0008  # approximate TDB–UTC conversion, <10 s residual HJD–BJD
    idx = np.unique(np.floor((t-t.min())/.8),return_inverse=True)[1]
    blocks = np.eye(idx.max()+1)[idx]
    X = np.column_stack((blocks, blocks*((t-t.min())%.8-.4)[:,None],orbit(t,9)))
    residual = y-X@np.linalg.lstsq(X,y,rcond=1e-9)[0]
    freq_grid = np.arange(.1,15,.0033)
    power = LombScargle(t,np.clip(residual,-5*np.std(residual),5*np.std(residual))).power(freq_grid)
    peaks = find_peaks(power,distance=8)[0]
    local = peaks[abs(freq_grid[peaks]-freq)<.025]
    result = {'cadences': len(t), 'highest_per_day': float(freq_grid[peaks[np.argmax(power[peaks])]])}
    if len(local):
        k = local[np.argmax(power[local])]
        result.update(local_peak_per_day=float(freq_grid[k]),rank=int(1+sum(power[peaks]>power[k])))
    return result


def analyze():
    datasets = list(series())
    _,t,y,e = datasets[0]
    res = fit(t,y,e)[3]
    freq_grid = np.arange(.1,12,1/(4*np.ptp(t)))
    ls = LombScargle(t, np.clip(res,-5*np.std(res),5*np.std(res)),e)
    power = ls.power(freq_grid)
    peak = freq_grid[np.argmax(power)]
    def loss(f):
        return fit(t,y,e,f)[2]
    def refine(center):
        return minimize_scalar(loss,bounds=(center-.00015,center+.00015),method='bounded',
                               options={'xatol':1.e-10})
    selected = refine(peak)
    f = float(selected.x)
    coeff = fit(t,y,e,f)[0][-2:]
    unit = coeff/np.linalg.norm(coeff)
    projections = {}
    for label,t,y,e in datasets:
        coeff,cov,_,_ = fit(t,y,e,f)
        projections[label] = {'rows':len(t),'nights':len(np.unique(np.floor(t))),
                              'semi_mag': float(unit@coeff[-2:]),
                              'cluster_se_mag': float(np.sqrt(unit@cov[-2:,-2:]@unit))}
    # Re-evaluate alternatives only on the discovery series.
    _,t,y,e = datasets[0]
    alias = {str(d): float(refine(f+d).fun-selected.fun) for d in [-1.0027379,-1,1,1.0027379]}
    return {'frequency_per_day':f,'period_hours':24/f,'projections':projections,'alias_delta_loss':alias,
            'tess':{path.name.split('.')[0]:tess(path,f) for path in sorted(DATA.glob('tess-*.csv.gz'))}}


if __name__ == '__main__':
    print(json.dumps(analyze(),indent=2))
