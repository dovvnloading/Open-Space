"""Fit M51 source transitions and test the two fixed eclipse-period families."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import xlogy
from scipy.stats import binom

DAY = 86400
EGRESS_2012, INGRESS_2012, EGRESS_2021 = 13813, 13814, 23479
INGRESS_2001 = 1622


def visits(folder):
    data = {}
    for file in folder.glob('*-gti.csv'):
        obs = int(file.name.split('-')[0])
        intervals = pd.read_csv(file).iloc[:, 2:4].to_numpy()
        table = pd.read_csv(folder / f'{obs}.csv.gz')
        data[obs] = (table[table.energy_ev <= 7000], intervals)
    return data


def exposure_at(t, intervals):
    return np.clip(t[:, None] - intervals[:, 0], 0, np.diff(intervals).ravel()).sum(axis=1)


def step_fit(times, intervals):
    times = np.sort(times)
    starts = np.arange(intervals[0, 0], intervals[-1, 1], 10 / DAY)
    left = exposure_at(starts, intervals)
    right = np.diff(intervals).sum() - left
    keep = (left > 0.008) & (right > 0.008)
    if not keep.any():
        return None
    left, right, starts = left[keep], right[keep], starts[keep]
    counts = np.searchsorted(times, starts)
    rest = len(times) - counts
    score = xlogy(counts, counts / left) + xlogy(rest, rest / right)
    best = np.argmax(score)
    near = starts[score >= score[best] - 1.92]
    fraction = left[best] / (left[best] + right[best])
    n, m = counts[best], rest[best]
    llr = xlogy(n, n / max((n + m)*fraction, 1e-12)) + xlogy(m, m / max((n + m)*(1-fraction), 1e-12))
    return dict(edge=starts[best], profile_low=near.min(), profile_high=near.max(),
                before=n, after=m, before_ks=left[best]*86.4, after_ks=right[best]*86.4,
                two_delta_log_likelihood=2*llr, cuts=len(starts))


def counts(table, time, mask):
    early = table.mjd_tdb_geo < time
    return int(sum(mask & early)), int(sum(mask & ~early))


def eclipse_fraction(intervals, epoch, period, width):
    ends = epoch + period*np.arange(np.floor((intervals[0, 0]-epoch)/period),
                                    np.ceil((intervals[-1, 1]-epoch+width)/period)+1)
    live = sum(np.maximum(0, np.minimum(stop, ends)-np.maximum(start, ends-width)).sum()
               for start, stop in intervals)
    return np.clip(live/np.diff(intervals).sum(), 0, 1)


def alias_grid(data, anchors):
    start, ingress, stop = [anchors[k] for k in (EGRESS_2012, INGRESS_2012, EGRESS_2021)]
    rows = []
    for family in (1, 2):
        # The 2012 pair must span one (long family) or two (short family) cycles.
        for cycles in range(230, 540):
            period = (stop-start)/cycles
            duration = family*period-(ingress-start)
            if not 1.04 <= duration <= 1.70:
                continue
            scores = {r: 0. for r in (0.02, 0.05, 0.15)}
            for obs, (tab, intervals) in data.items():
                if obs in anchors:
                    continue
                times = tab.loc[tab.r1_arcsec <= 1.8, 'mjd_tdb_geo'].to_numpy()
                frac = eclipse_fraction(intervals, start, period, duration)
                num = np.count_nonzero((times-start) % period > period-duration)
                for ratio in scores:
                    # Conditional-on-total likelihood fits an independent normalization per visit.
                    scores[ratio] += num*np.log(ratio) - len(times)*np.log1p(frac*(ratio-1))
            err = (ingress-anchors[INGRESS_2001]+period/2) % period-period/2
            rows.append(dict(cycles_2012=family, cycles_2012_2021=cycles, period_days=period,
                             eclipse_days=duration, predicted_minus_observed_2001_ingress_d=err,
                             **{f'other_visits_logLR_{ratio}':s for ratio, s in scores.items()}))
    return rows


def run(folder, out):
    out.mkdir(parents=True, exist_ok=True)
    data = visits(folder)
    edges, summaries, controls = [], [], []
    anchors = {}
    multiplicity = 0
    for obs, (tab, g) in sorted(data.items()):
        summaries.append(dict(obsid=obs, start=g[0,0], stop=g[-1,1], gti_ks=np.diff(g).sum()*86.4,
                              n_ulx1=sum(tab.r1_arcsec<=1.8), n_ulx2=sum(tab.r2_arcsec<=1.8)))
        for source, column in [(1, tab.r1_arcsec), (2, tab.r2_arcsec)]:
            fit = step_fit(tab.loc[column<=1.8, 'mjd_tdb_geo'].to_numpy(), g)
            if fit is None:
                continue
            multiplicity += fit['cuts']
            edges.append(dict(obsid=obs, source=source, **fit))
            if source==1 and obs in (INGRESS_2001, EGRESS_2012, INGRESS_2012, EGRESS_2021):
                anchors[obs] = fit['edge']
                masks = {f'ulx1_{r}': tab.r1_arcsec<=r for r in (1.2, 1.8, 2.5)}
                masks['ulx2_1.8'] = tab.r2_arcsec<=1.8
                masks['annulus_3_6'] = (tab.r1_arcsec>3)&(tab.r1_arcsec<6)&(tab.r2_arcsec>2.5)
                for region, mask in masks.items():
                    n,m = counts(tab, fit['edge'], mask)
                    controls.append(dict(obsid=obs, region=region, before=n, after=m))
    for name, rows in [('visits', summaries), ('step-fits',edges), ('controls',controls),
                       ('period-grid',alias_grid(data,anchors))]:
        pd.DataFrame(rows).to_csv(out / f'{name}.csv', index=False, float_format='%.9f')
    best = next(e for e in edges if e['source']==1 and e['obsid']==EGRESS_2021)
    exposure_fraction = best['before_ks'] / (best['before_ks']+best['after_ks'])
    p = binom.cdf(best['before'],best['before']+best['after'],exposure_fraction)
    print(f"Constant-rate model, two tails times {multiplicity} tested cuts: p <= {2*multiplicity*p:.3g}")
    plot_visit(data[EGRESS_2021], best['edge'], out)


def plot_visit(visit, edge, out):
    tab, g = visit
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 3.9), sharex=True, constrained_layout=True)
    bins = np.r_[np.arange(g[0,0],g[-1,1],2000/DAY),g[-1,1]]
    for ax,r,title,color in zip(axes,[tab.r1_arcsec,tab.r2_arcsec],
                                ['M51 ULX-1','M51 ULX-2 (nearby control)'],['#9e341a','#256583']):
        hist,_=np.histogram(tab.loc[r<=1.8,'mjd_tdb_geo'],bins)
        exposure=np.diff(exposure_at(bins,g))*86.4
        ax.stairs(hist/exposure,(bins-g[0,0])*86.4,color=color,linewidth=1.7)
        ax.axvline((edge-g[0,0])*86.4,ls='--',color='0.45',linewidth=1)
        ax.set_title(title, loc='left', fontsize=10)
        ax.set(ylabel='Counts / ks', ylim=(0,None))
        ax.spines[['top','right']].set_visible(False)
    axes[-1].set_xlabel('Kiloseconds from start of Chandra ObsID 23479')
    fig.savefig(out/'june-egress.png',dpi=140)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('photons', type=Path)
    parser.add_argument('out', type=Path)
    args = parser.parse_args()
    run(args.photons,args.out)
