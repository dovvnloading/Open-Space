"""Reproduce eclipse spacings and the independent preselected-frequency tests."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.ndimage import percentile_filter
from scipy.optimize import least_squares
from astropy.timeseries import LombScargle

HERE = Path(__file__).resolve().parent
PERIOD = 7.281
REF = np.array([.032, .436])
RNG = np.random.default_rng(196989952)

def center(profile, phase, guess=None):
    time, flux = profile
    x = ((time / PERIOD - phase + .5) % 1 - .5) * PERIOD
    observed = np.abs(x) < .8
    x, flux = x[observed], flux[observed]
    def diff(b):
        mid, depth, width, shape, level, slope = b
        return flux - (level + slope*x - depth*np.exp(-.5*abs((x-mid)/width)**shape))
    initial = [0, .45, .08, 3, 1, 0] if guess is None else guess
    return least_squares(diff, initial, bounds=([-.6,.03,.015,.8,.8,-.15],[.6,.95,.3,9,1.2,.15]),
                         loss="soft_l1", f_scale=.012, max_nfev=90).x

def binary_phases(t):
    # Linear interpolation of measured eclipse phases used only to exclude entire eclipse windows.
    drift = np.array([-.0000042, -.0000013])
    return np.array([.03198,.4371]) + (t[:,None]-2950)*drift

def outside(t):
    delta = np.abs((t[:,None]/PERIOD-binary_phases(t)+.5)%1-.5)
    return np.all(delta > np.array([.043,.025]), axis=1)

def frequency_fit(t, y, freq, errors=None):
    time = t-t.mean()
    columns = [np.ones(len(t))]
    for f in freq:
        columns.extend((np.sin(2*np.pi*f*time), np.cos(2*np.pi*f*time)))
    A = np.stack(columns,axis=1)
    weight = np.ones(len(t)) if errors is None else 1/errors
    beta = np.linalg.lstsq(A*weight[:,None], y*weight,rcond=None)[0]
    return A@beta, np.linalg.norm(beta[1:].reshape(-1,2),axis=1)

def inputs():
    tess = pd.read_csv(HERE/'tess.csv.gz')
    wasp = pd.read_csv(HERE/'wasp.csv.gz')
    wasp = wasp[(wasp['magnitude error']<.035)&wasp.magnitude.between(10.5,13.7)].copy()
    wasp['btjd'] = wasp.HJD-2457000
    wasp['flux'] = 10**(-.4*(wasp.magnitude-wasp.groupby('camera').magnitude.transform('median')))
    gaia = pd.read_csv(HERE/'gaia.csv')
    return tess, wasp, gaia

def normalize(sector):
    t=sector.btjd.to_numpy();y=sector.flux.to_numpy()
    return t,y/np.median(y)

def timing(tess, wasp):
    groups=[('SuperWASP',wasp)]
    groups.extend((f'TESS {s}',g) for s,g in tess.groupby('sector'))
    rows=[]
    for tag,group in groups:
        t,y=normalize(group);coefs=[center((t,y),ph) for ph in REF]
        mid=REF+np.array([b[0] for b in coefs])/PERIOD
        # Resample complete nights for WASP and complete binary cycles for TESS.
        # The small TESS event count makes this a sensitivity interval.
        keys=np.floor(t+0.5) if tag=='SuperWASP' else np.floor(t/PERIOD+.2)
        blocks=[np.where(keys==k)[0] for k in np.unique(keys)]
        replicates=[]
        for draw in range(250):
            k=np.concatenate([blocks[i] for i in RNG.integers(len(blocks),size=len(blocks))])
            if any(np.sum(np.abs((t[k]/PERIOD-m+.5)%1-.5)<.01)<5 for m in mid):
                continue
            samples=[center((t[k],y[k]),ph,b)[0] for ph,b in zip(REF,coefs)]
            replicates.append((REF[1]-REF[0])*PERIOD+samples[1]-samples[0])
        rows.append(dict(group=tag,mean_btjd=float(t.mean()),cadences=len(t),
                         centers_phase=mid.tolist(),separation_days=float((mid[1]-mid[0])*PERIOD),
                         bootstrap_p16_p84=np.quantile(replicates,[.16,.84]).tolist(),
                         resamples_used=len(replicates)))
    return rows

def modes(tess):
    out=[]
    for _,g in tess.groupby('sector'):
        t,y=normalize(g)
        # Only remove drifts on timescales longer than the candidate waves.
        cadence=np.median(np.diff(t));window=round(3/cadence)//2*2+1
        base=percentile_filter(y,percentile=65,size=window)
        use=outside(t)&(abs(y/base-1)<.1)
        out.append((t[use],y[use]/base[use]-1))
    t,y=np.hstack([a[0] for a in out]),np.hstack([a[1] for a in out]);residual=y.copy()
    grid=np.arange(.4,4,.00002);freq=[]
    for _ in range(2):
        power=LombScargle(t,residual).power(grid)
        freq.append(float(grid[np.argmax(power)]))
        model,amp=frequency_fit(t,y,freq)
        residual=y-model
    fitted={str(s):frequency_fit(t,y,freq)[1].tolist() for s,(t,y) in zip(sorted(tess.sector.unique()),out)}
    return freq,amp.tolist(),fitted

def ground(wasp, freq):
    t=wasp.btjd.to_numpy();y=wasp.flux.to_numpy()
    use=outside(t)&(abs(y-1)<.09);t,y=t[use],y[use]
    ls=LombScargle(t,y);scan=np.arange(.5,3,.00005);power=ls.power(scan)
    model,amp=frequency_fit(t,y,freq)
    return dict(retained=len(t),peak_cpd=float(scan[np.argmax(power)]),preselected_amplitude=amp.tolist(),
                fixed_power=[float(ls.power(f)) for f in freq])

def gaia_test(gaia,freq):
    t=gaia.TimeG.to_numpy()-1802.5
    use=(gaia.GrVFlag==0).to_numpy()&outside(t)&np.isfinite(gaia.Gmag)
    t=t[use]; y=gaia.Gmag.to_numpy()[use]
    # Gaia visits minutes apart are correlated. Average each UTC-like observing day
    # in this BJD-based series, then randomize the daily magnitudes across actual days.
    day=np.floor(t+.5)
    df=pd.DataFrame({'t':t,'mag':y,'day':day}).groupby('day',as_index=False).mean()
    t=df.t.to_numpy();y=df.mag.to_numpy();y-=y.mean()
    fitted,amp=frequency_fit(t,y,freq)
    r2=1-np.sum((y-fitted)**2)/np.sum(y*y)
    design=np.column_stack((np.ones(len(t)),*[v for f in freq for v in
                    (np.sin(2*np.pi*f*t),np.cos(2*np.pi*f*t))]))
    Q=np.linalg.qr(design)[0];exceed=0;draws=30000
    for _ in range(draws//1000):
        null=np.array([RNG.permutation(y) for _ in range(1000)])
        ratios=np.sum((null@Q)**2,axis=1)/np.sum(y*y)
        exceed+=int(np.sum(ratios>=r2))
    return dict(transits=int(np.sum(use)),daily_groups=len(t),amplitude_mag=amp.tolist(),
                joint_r2=float(r2),permutation_p=(1+exceed)/(1+draws))
if __name__=='__main__':
    tess,wasp,gaia=inputs()
    freq,amp,persector=modes(tess)
    result=dict(tess_frequencies_cpd=freq,tess_semiamplitudes_flux=amp,tess_sector_amplitudes_flux=persector,
                superwasp=ground(wasp,freq),gaia=gaia_test(gaia,freq),spacings=timing(tess,wasp))
    (HERE/'measurements.json').write_text(json.dumps(result,indent=2)+'\n')
