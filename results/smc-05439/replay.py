"""Test the full, alternating secondary period in independent survey photometry."""
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.signal import find_peaks

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
EPOCH=9000. # HJD 2459000.0
START=np.array([1/1.4511308,1/(2*.582732)])


def harmonics(time, freq, count):
    phase=2*np.pi*(time-EPOCH)*freq
    return np.column_stack([v for n in range(1,count+1) for v in (np.cos(n*phase),np.sin(n*phase))])


def design(t, fs):
    return np.column_stack((np.ones(len(t)),harmonics(t,fs[0],9),harmonics(t,fs[1],8)))


def fit(x,y,e,covariance_time=None):
    w=1/e
    for _ in range(4):
        b=np.linalg.lstsq(x*w[:,None],y*w,rcond=None)[0]
        residual=y-x@b
        scale=max(.003,1.483*np.median(abs(residual-np.median(residual))))
        w=1/e/np.sqrt(np.maximum(1.,abs(residual)/(3*scale)))
    loss=np.sum(np.where(abs(residual)<3*scale,residual**2,6*scale*abs(residual)-9*scale**2)/e**2)
    if covariance_time is None:
        return b,loss
    labels=np.unique(np.floor(covariance_time),return_inverse=True)[1]
    score=x*(residual*w**2)[:,None]
    meat=np.stack([np.bincount(labels,weights=c) for c in score.T])
    bread=np.linalg.pinv((x*w[:,None]).T@(x*w[:,None]))
    return b,bread@meat@meat.T@bread


def records():
    for era in (4,3):
        for band in ('I','V'):
            t,y,e=np.loadtxt(DATA/f'ogle{era}-{band}.dat.gz',unpack=True)
            ok=np.isfinite(t+y+e)&(e>0)&(e<min(.13,max(.03,3*np.median(e))))
            t,y,e=t[ok],y[ok],e[ok]
            yield f'OGLE-{era} {band}',t-2450000*(era==4),y-np.median(y),np.maximum(e,.7*np.median(e))
    with gzip.open(DATA/'xtime.tsv.gz','rt') as fp:
        lines=[l for l in fp if l.strip() and not l.startswith('#')]
    rows=list(csv.DictReader(lines[:1]+lines[3:],delimiter='\t',skipinitialspace=True))
    for obj in sorted({r['EROS'].strip() for r in rows}):
        for band in ('R','B'):
            vals=[]
            for r in rows:
                if r['EROS'].strip()!=obj:
                    continue
                try:
                    t,y,e=(float(r[i]) for i in ('HJD',band+'mag','e_'+band+'mag'))
                except ValueError:
                    continue
                if np.isfinite(t+y+e) and 0<e<.2 and y<40:
                    vals.append([t,y,e])
            if len(vals)>30:
                t,y,e=np.array(vals).T
                yield f'EROS {obj} {band}',t,y-np.median(y),np.maximum(e,.6*np.median(e))


def training(rows):
    _,t,y,e=rows[0]
    def obj(delta):
        return fit(design(t,START+delta*1.e-6),y,e)[1]
    opt=minimize(obj,[0.,0.],method='Nelder-Mead',
                 options={'initial_simplex':np.array([[0.,0.],[1.,0.],[0.,4.]]),'xatol':.001})
    fs=START+opt.x*1.e-6
    b,loss=fit(design(t,fs),y,e)
    return fs,b


def templates(t,fs,b):
    x=design(t,fs)
    even=np.repeat(np.arange(1,9)%2==0,2)
    h=x[:,-16:]
    return np.column_stack((x[:,1:19]@b[1:19],h[:,even]@b[-16:][even],h[:,~even]@b[-16:][~even]))


def main():
    rows=list(records())
    fs,b=training(rows)
    phi=np.arange(10000)/10000
    g=harmonics(EPOCH+phi/fs[1],fs[1],8)@b[-16:]
    peaks=find_peaks(np.tile(g,3),prominence=.004)[0]
    centers=(peaks[(peaks>=10000)&(peaks<20000)]-10000)
    out=dict(epoch_hjd=EPOCH+2450000,periods_days=(1/fs).tolist(),
             secondary_relative_dim_peaks=[{'phase':float(phi[i]),'mag_above_brightest':float(g[i]-min(g))} for i in centers],
             rows={})
    for name,t,y,e in rows:
        x=design(t,fs)
        a=np.column_stack((x[:,:19],templates(t,fs,b)[:,1:]))
        coeff,cov=fit(a,y,e,t)
        out['rows'][name]=dict(n=len(t),nights=len(np.unique(np.floor(t))),scale_even_odd=coeff[-2:].tolist(),se_even_odd=np.sqrt(cov.diagonal()[-2:]).tolist())
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
