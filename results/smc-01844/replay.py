"""Measure the short-period component and project its fixed phase into other surveys."""
from pathlib import Path
import csv
import gzip
import json

import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
EPOCH = 10000.0436  # relative to HJD 2450000
BINARY_PERIOD = 1.8398765
SEARCH_BOUNDS = (4.91210, 4.91235)  # preliminary OGLE periodogram peak, cycles/day


def design(time, freq, harmonics):
    phase = 2*np.pi*(time - EPOCH)
    orbit = [column for n in range(1, harmonics+1)
             for column in (np.cos(n*phase/BINARY_PERIOD), np.sin(n*phase/BINARY_PERIOD))]
    return np.column_stack([np.ones(len(time)), *orbit, np.cos(phase*freq), np.sin(phase*freq)])


def fit(time, magnitude, error, freq, harmonics=11):
    x = design(time, freq, harmonics)
    weights = 1/error
    for _ in range(4):
        coeff = np.linalg.lstsq(x*weights[:, None], magnitude*weights, rcond=None)[0]
        residual = magnitude - x@coeff
        scatter = max(.003, 1.483*np.median(abs(residual-np.median(residual))))
        weights = 1/error/np.sqrt(np.maximum(1., abs(residual)/(3*scatter)))
    score = x * (residual * weights**2)[:, None]
    labels = np.unique(np.floor(time), return_inverse=True)[1]
    grouped = np.stack([np.bincount(labels, weights=col) for col in score.T])
    bread = np.linalg.pinv((x*weights[:, None]).T @ (x*weights[:, None]))
    cov = bread @ grouped @ grouped.T @ bread
    scaled = abs(residual/scatter)
    loss = np.sum(np.where(scaled<3, scaled**2, 6*scaled-9)*scatter**2/error**2)
    return coeff, cov, loss, residual


def query_rows(name):
    with gzip.open(DATA/name, 'rt') as handle:
        lines = [s for s in handle if s.strip() and not s.startswith('#')]
    return list(csv.DictReader(lines[:1]+lines[3:], delimiter='\t', skipinitialspace=True))


def series():
    for era in ('iv', 'iii'):
        for band in ('I', 'V'):
            t, y, e = np.loadtxt(DATA/f'{era}-{band}.dat.gz', unpack=True)
            keep = np.isfinite(t+y+e) & (e>0) & (e<min(.13, max(.03, 3*np.median(e))))
            t,y,e=t[keep],y[keep],e[keep]
            if era=='iv':
                t-=2450000
            yield f'OGLE-{era.upper()} {band}', t, y-np.median(y), np.maximum(e,.7*np.median(e)), 11

    rows = query_rows('xtime.tsv.gz')
    for obj in sorted({row['EROS'].strip() for row in rows}):
        for band in ('R','B'):
            values=[]
            for row in rows:
                if row['EROS'].strip()!=obj:
                    continue
                try:
                    t,y,e = (float(row[k]) for k in ('HJD',band+'mag','e_'+band+'mag'))
                except ValueError:
                    continue
                if np.isfinite(t+y+e) and 0<e<.20 and y<40:
                    values.append((t,y,e))
            if len(values)>=30:
                t,y,e=np.asarray(values).T
                yield f'EROS {obj} {band}', t, y-np.median(y), np.maximum(e,.6*np.median(e)), 11

    with gzip.open(DATA/'mj-table2.tsv.gz','rt') as handle:
        rows=list(csv.DictReader(handle, delimiter='\t',skipinitialspace=True))
    for band in ('u','VJ','IC'):
        values=[(float(r['HJD'])-2450000,-2.5*np.log10(float(r['light'])))
                for r in rows if r['Filt'].strip()==band]
        t,y=np.asarray(values).T
        yield f'Mt John {band}', t, y, np.ones(len(t)), 8


def analyze():
    seq=list(series())
    _,t,y,e,h=seq[0]
    opt=minimize_scalar(lambda f:fit(t,y,e,f,h)[2],bounds=SEARCH_BOUNDS,
                        method='bounded',options={'xatol':1.e-10})
    freq=float(opt.x)
    coeff=fit(t,y,e,freq,h)[0][-2:]
    unit=coeff/np.linalg.norm(coeff)
    out=dict(binary_period_d=BINARY_PERIOD,epoch_hjd=EPOCH+2450000,
             frequency_per_d=freq,period_hours=24/freq,
             training_cos_sin_mag=coeff.tolist(),results={})
    for name,t,y,e,h in seq:
        n=len(t)
        b,c,loss,residual=fit(t,y,e,freq,h)
        amp=unit@b[-2:];se=np.sqrt(unit@c[-2:,-2:]@unit)
        out['results'][name]=dict(n=n,nights=len(np.unique(np.floor(t))),
                min_hjd=float(min(t)+2450000),max_hjd=float(max(t)+2450000),
                signed_semiamplitude_mag=float(amp),cluster_se_mag=float(se),
                cosine_sine_mag=b[-2:].tolist(),model_harmonics=h)
    return out


if __name__=='__main__':
    print(json.dumps(analyze(),indent=2))
