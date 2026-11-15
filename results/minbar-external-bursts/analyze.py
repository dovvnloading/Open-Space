"""Replay the survey, independently clocked profiles and controls from frozen counts."""
import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from geometry import geometry

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'


def rows(name):
    with (DATA / name).open() as stream:
        return list(csv.DictReader(stream))


def grouped(name):
    groups = defaultdict(list)
    for row in rows(name):
        groups[int(row['seq'])].append(row)
    return groups


def timescreen(events, epochs, offset, width=10):
    """Count all cross-pairs inside the inclusive half-window, allowing duplicates."""
    return int(np.sum(np.searchsorted(epochs, events+offset+width, side='right') -
                      np.searchsorted(epochs, events+offset-width)))


def waveform(pca_rows, kw_rows):
    # Relative time at each spacecraft; clock difference is fitted only in the
    # +/-10 s lag search. Compare all 45 consecutive seconds in each window.
    kw_t = np.array([float(r['start_s']) for r in kw_rows])
    kw_counts = np.array([float(r['counts']) for r in kw_rows])
    second = np.arange(46)
    # The archived Konus bin stamps are start times separated by 64 ms.
    # Average counts per native bin to handle 15 or 16 bins in a second.
    kval = (np.histogram(kw_t,second,weights=kw_counts)[0] /
            np.histogram(kw_t,second)[0])
    grid = np.arange(-704,704,0.25)
    idx = np.rint((np.array([float(r['start_s']) for r in pca_rows])+704)*4).astype(int)
    answers = {}
    for band in ('all','hard','soft'):
        a = np.full(len(grid), np.nan)
        a[idx] = [int(r[band]) for r in pca_rows]
        def packet(shift):
            i = round((704+shift)*4)
            return a[i:i+180].reshape(45,4).sum(axis=1)
        def cor(v,w):
            if not np.isfinite(w).all():
                return float('nan')
            a0, b0 = v - v.mean(), w - w.mean()
            scale = np.linalg.norm(a0) * np.linalg.norm(b0)
            return np.dot(a0,b0)/scale if scale else float('nan')
        def fit(shift):
            ls = np.arange(-10,10.25,0.25)
            rr = np.array([cor(kval,packet(shift+lag)) for lag in ls])
            if np.all(np.isnan(rr)):
                return None
            j = np.nanargmax(rr)
            return float(ls[j]), float(rr[j])
        fitted = fit(0)
        if fitted is None:
            answers[band] = {'status':'no complete 45 s window within lag range'}
            continue
        lag, best = fitted
        control = [fit(s) for s in list(range(-600,-59,60))+list(range(60,601,60))]
        nulls = [x[1] for x in control if x]
        # Remove a common burst envelope as a separate shape check.
        corr_diff = cor(np.diff(kval),np.diff(packet(lag)))
        answers[band] = dict(zero_r=round(float(cor(kval,packet(0))),4),
                            best_lag_s=lag, best_r=round(best,4),
                            diff_r_at_best=round(float(corr_diff),4),
                            control_windows=len(nulls),max_control_r=round(max(nulls),4))
    return answers


def run():
    minbar = rows('minbar.csv')
    meta = rows('events.csv')
    event = grouped('rxte_profiles.csv')
    trigger = grouped('konus_profiles.csv')
    mjd_unix = 40587
    t = np.array([(float(r['Time'])-mjd_unix)*86400 for r in minbar])
    kw = rows('konus_triggers.csv')
    ts = np.array([datetime.fromisoformat(r['utc']).replace(tzinfo=timezone.utc).timestamp() for r in kw])
    grbs = np.array([r['type'] in ('GRB','sGRB') for r in kw])
    offsets = [k*60 for k in range(-720,721) if k]
    null = np.array([timescreen(t,ts[grbs],x) for x in offsets])
    near = []
    for tm,m in zip(t,minbar):
        lo,hi=np.searchsorted(ts,[tm-35,tm+35],side='left')
        for q,k in zip(ts[lo:hi],kw[lo:hi]):
            near.append([int(m['Seq']),m['Name'],m['sflag'],k['utc'],k['type'],round(tm-q,2)])
    trains = defaultdict(list)
    for m in minbar:
        if m['Name']!='IGR J17480-2446':
            trains[m['Name'],m['obsid'],m['Inst']].append(m)
    triples = []
    for key,data in trains.items():
        data.sort(key=lambda r:float(r['Time']))
        for a,b,c in zip(data,data[1:],data[2:]):
            gaps = [(float(y['Time'])-float(x['Time']))*1440 for x,y in [(a,b),(b,c)]]
            if max(gaps)<=45:
                triples.append([round(sum(gaps),4),[int(z['Seq']) for z in (a,b,c)],key[0]])
    triples.sort()
    results={'minbar_rows':len(minbar),'pca_rows':sum(r['Inst'].startswith('XP') for r in minbar),
             'konus_rows':len(kw),'konus_grbs':int(grbs.sum()),'matches_35s':sorted(near),
             'grb_count_10s':timescreen(t,ts[grbs],0),
             'time_shifts':len(null),'shift_count_histogram':dict(zip(
                 map(str,np.unique(null)),map(int,np.unique(null,return_counts=True)[1]))),
             'ranked_three_event_windows':triples,'waveforms':{}}
    for row in meta:
        seq=int(row['seq'])
        results['waveforms'][str(seq)] = waveform(event[seq],trigger[seq])
    results['geometric_checks'] = geometry()
    return results


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
