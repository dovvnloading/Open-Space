"""Evaluate the released NEOWISE periods at the actual SPHEREx exposure times."""
import csv
import json
import numpy as np

from replay import ROOT, predictions, tablewrite


def predicted_ratios(visits, period, dates):
    w2 = predictions(visits, period, dates.ravel())[:,1].reshape(dates.shape)
    return np.median(10**(-0.4 * (w2[:,1:]-w2[:,0,None])), axis=0)


def leave_one_out(visits, period):
    # Keep each solution in its original frequency neighborhood. Both bands
    # share the period; W2 alone supplies the subsequent flux prediction.
    frequencies = np.linspace(1/2000, 1/50, 12000)
    window = abs(frequencies - 1/period) <= 0.5/np.ptp(visits[:,0])
    trial = 1/frequencies[window]
    fits = []
    for omitted in range(len(visits)):
        v = np.delete(visits, omitted, axis=0)
        centered = v[:,1:3] - v[:,1:3].mean(axis=0)
        norm = np.sum(centered**2, axis=0)
        def score(p):
            resid = v[:,1:3] - predictions(v, p, v[:,0])
            return np.mean(1-np.sum(resid**2, axis=0)/norm)
        refit = max(trial, key=score)
        fits.append((v, refit))
    return fits


def main():
    visits = np.genfromtxt(ROOT/'visits.csv',delimiter=',',skip_header=1)
    with (ROOT/'six-detector-flux.csv').open() as f:
        exposures = {(r['detector'],r['obs_id']): r for r in csv.DictReader(f)}
    with (ROOT/'six-detector-matches.csv').open() as f:
        groups = list(csv.DictReader(f))
    aliases = json.loads((ROOT/'result.json').read_text())['alias_fits']
    result = []
    for alias in aliases:
        period = alias['period_d']
        fits = leave_one_out(visits,period)
        for d in ('5','6'):
            rows = [[exposures[d,t[key+'_id']] for key in ('first','middle','late')]
                    for t in groups if t['detector']==d]
            dates = np.array([[float(r['mjd']) for r in group] for group in rows])
            flux = np.array([[float(r['flux_jy']) for r in group] for group in rows])
            obs = np.median(flux[:,1:]/flux[:,0,None],axis=0)
            expected = predicted_ratios(visits,period,dates)
            repeats = np.array([predicted_ratios(v,p,dates) for v,p in fits])
            for j,visit in enumerate(('middle','late')):
                result.append([d,f'{period:.8f}',visit,len(rows),f'{obs[j]:.6f}',f'{expected[j]:.6f}',
                               f'{repeats[:,j].min():.6f}', f'{repeats[:,j].max():.6f}'])
    tablewrite(ROOT/'alias-extrapolation.csv',
               ['detector','baseline_period_d','visit_relative_to_first','triplets','observed_ratio',
                'predicted_w2_ratio','leave_one_visit_out_min','leave_one_visit_out_max'],result)


if __name__ == '__main__':
    main()
