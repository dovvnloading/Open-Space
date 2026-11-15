"""Reproduce the Gaia19dyr nightly photometry analysis. Requires numpy, pandas."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
N_DRAW = 20000
JD_REF = 2459000.0


def center_by_year(values, years):
    values = np.array(values, copy=True)
    for year in np.unique(years):
        m = years == year
        values[..., m] -= values[..., m].mean(axis=-1, keepdims=True)
    return values


def scan(frame, period_lo=2.0, period_hi=30.0):
    """Fit sine, cosine and a separate floating mean for each observing year."""
    t = frame.jd.to_numpy() - JD_REF
    y = frame.delta_mag.to_numpy()
    years = frame.year.to_numpy()
    centered = center_by_year(y, years)
    freqs = np.arange(1 / period_hi, 1 / period_lo, 1 / (np.ptp(t) * 24))
    best = None
    for frequency in np.array_split(freqs, max(1, len(freqs) // 700)):
        a = 2 * np.pi * frequency[:, None] * t
        c = center_by_year(np.cos(a), years)
        s = center_by_year(np.sin(a), years)
        cc, ss, cs = (c*c).sum(1), (s*s).sum(1), (c*s).sum(1)
        yc, ys = c @ centered, s @ centered
        det = cc * ss - cs * cs
        bc = (yc * ss - ys * cs) / det
        bs = (ys * cc - yc * cs) / det
        gains = bc * yc + bs * ys
        j = np.argmax(gains)
        trial = (gains[j], frequency[j], np.arctan2(bs[j], bc[j]) / (2*np.pi))
        if best is None or trial[0] > best[0]:
            best = trial
    return dict(period_d=1/best[1], phase_cycle=best[2], r2=best[0]/(centered @ centered))


def phase_test(frame, fit):
    """Rotate phase independently per calendar year, preserving intra-year patterns."""
    year = frame.year.to_numpy()
    y = center_by_year(frame.delta_mag.to_numpy(), year)
    phase = (frame.jd.to_numpy()-JD_REF)/fit['period_d']-fit['phase_cycle']
    z = np.exp(2j*np.pi*phase)
    sums = np.array([(y[year==v] * z[year==v]).sum() for v in np.unique(year)])
    observed = sums.real.sum()
    rng = np.random.default_rng(19478)
    offsets = rng.uniform(0, 2*np.pi, size=(N_DRAW,len(sums)))
    null = (np.exp(1j*offsets) @ sums).real
    inside = np.abs((phase+.5)%1-.5) < .16
    values = frame.delta_mag.to_numpy()
    return dict(n=len(y), n_inside=int(inside.sum()),
                median_in_minus_out_mag=float(np.median(values[inside])-np.median(values[~inside])),
                seasonal_phase_randomization_p=(1+int((null>=observed).sum()))/(1+N_DRAW))


def run():
    d = pd.read_csv(HERE/'nightly.csv',dtype={'star':str})
    d['year'] = pd.to_datetime(d.mjd_night+2400000.5,origin='julian',unit='D').dt.year
    target = d[d.star=='target']
    r = target[target['filter']=='zr']
    active = d.year.between(2019,2022)
    fit = scan(r[r.year.between(2019,2022)])
    results = dict(fit=fit, active_checks={}, season_r=[], field_controls={}, leave_one_year_out=[])
    for band in ['zg','G']:
        q = target[(target['filter']==band) & target.year.between(2019,2022)]
        results['active_checks'][band] = phase_test(q,fit)
        if band=='G':
            # Midnight crossing can split a Gaia visit; preserve the original alert exclusion.
            results['active_checks']['G_without_alert'] = phase_test(q[abs(q.jd-2458731.2)>2],fit)
    for year, a in r.groupby('year'):
        b=scan(a)
        results['season_r'].append(dict(year=int(year), nights=len(a), deep_nights=int((a.delta_mag>1).sum()),
                                        q05_q95_span_mag=a.delta_mag.quantile(.95)-a.delta_mag.quantile(.05),
                                        strongest_period_d=b['period_d'],yearly_r2=b['r2']))
    for year in range(2019,2023):
        train=r[r.year.between(2019,2022)&(r.year!=year)]
        results['leave_one_year_out'].append(dict(omitted=year, **scan(train)))
    for star,a in d[active & (d['filter']=='zr') & (d.star!='target')].groupby('star'):
        results['field_controls'][star]=dict(phase_test=phase_test(a,fit),
             robust_span_mag=a.delta_mag.quantile(.95)-a.delta_mag.quantile(.05))
    return results


if __name__=='__main__':
    print(json.dumps(run(), indent=2))
