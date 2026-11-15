"""Analyze two independent Apollo event catalogs. See screen.md for the input URLs.

Usage: python analyze.py path/to/GradeA_thermal_mq_catalog_final.csv path/to/levent.csv
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

FULL = pd.to_datetime(["1976-09-08T12:52", "1976-10-08T04:56", "1976-11-06T23:15",
                       "1976-12-06T18:15", "1977-01-05T12:10", "1977-02-04T03:56",
                       "1977-03-05T17:13", "1977-04-04T04:09"])
CONTACTS = {2: ("1976-11-06T20:48:17", "1976-11-07T01:14:05"),
            7: ("1977-04-04T02:06:52", "1977-04-04T06:29:44")}
GRADESHA = "c3ee72f706b1f8bbbced6ee81f98ed10f3e23c8bb99616242f89a6ff2d3ba82d"
DEEPSHA = "d23bdf8e35e2ff95ac690d1f07c4ba5300bca80af4e405de3853776627daafdf"


def checked(path, sha):
    if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        raise ValueError(f"Input hash mismatch for {path}")
    return pd.read_csv(path)


def eclipses(catalog):
    catalog['time'] = pd.to_datetime(catalog.ft_arrival_time, format='mixed')
    events = catalog.groupby('evid').time.min().sort_values()
    bearing = catalog.groupby('evid').theta_mean.first()
    off_lander = events[bearing.reindex(events.index).notna() & ~bearing.reindex(events.index).between(80,120)]
    def count(lo, hi, sample=events):
        return int(sample.searchsorted(hi) - sample.searchsorted(lo))
    results = {}
    for index, contacts in CONTACTS.items():
        lo, hi = map(pd.Timestamp, contacts)
        width = hi - lo
        observed = count(lo, hi)
        controls = {}
        off_lander_controls = {}
        for j, phase in enumerate(FULL):
            if j not in CONTACTS:
                shift = phase - FULL[index]
                key = phase.date().isoformat()
                controls[key] = count(lo + shift, hi + shift)
                off_lander_controls[key] = count(lo + shift, hi + shift, off_lander)
        adjacent = [count(lo + k*width, hi + k*width) for k in (-2, -1, 1, 2)]
        results[lo.date().isoformat()] = dict(
            utc_contacts=contacts, count=observed, phase_controls=controls,
            upper_rank_fraction=(1+sum(n >= observed for n in controls.values()))/(1+len(controls)),
            adjacent_counts=adjacent,
            outside_lander_sector_count=count(lo,hi,off_lander),
            outside_lander_sector_controls=off_lander_controls,
            event_ids=events[(events >= lo) & (events < hi)].index.tolist())
    results['april_umbral_contacts'] = dict(utc=['03:30:58','05:05:42'],
        count=count(pd.Timestamp('1977-04-04 03:30:58'), pd.Timestamp('1977-04-04 05:05:42')))
    return dict(geophone_rows=len(catalog), events=len(events), tests=results)


def pairs(catalog, trials=4999):
    labels = catalog['2004Class'].str.replace(' ', '')
    select = labels.str.fullmatch('A[1-9][0-9]*', na=False)
    rows = catalog[select]
    dates = pd.to_datetime((1900+rows.Year).astype(str)+'-'+rows.DOY.astype(str), format='%Y-%j')
    minutes = dates.to_numpy().astype('datetime64[m]').astype('int64')
    minutes += (rows.StartTime//100 * 60 + rows.StartTime%100).to_numpy()
    order = np.argsort(minutes, kind='stable')
    minutes = minutes[order]
    nests, labels = np.unique(labels[select].to_numpy()[order], return_inverse=True)
    base = len(nests)
    # Edge endpoints and ordering are fixed, including pairs straddling midnight.
    count = np.searchsorted(minutes, minutes+360, side='right')-np.arange(len(minutes))-1
    left = np.repeat(np.arange(len(minutes)), count)
    right = np.concatenate([np.arange(i+1,i+1+n) for i,n in enumerate(count)])
    delta = minutes[right]-minutes[left]
    sort = np.argsort(delta)
    left, right, delta = left[sort], right[sort], delta[sort]
    lags = [5, 20, 60, 360]
    stops = np.searchsorted(delta, lags, side='right')

    def measure(ids):
        x, y = ids[left], ids[right]
        pair_ids = np.minimum(x,y)*base + np.maximum(x,y)
        maxima, totals, winners = [], [], []
        for stop in stops:
            codes = pair_ids[:stop][x[:stop] != y[:stop]]
            totals.append(len(codes))
            counts = np.bincount(codes, minlength=base*base)
            peak = counts.max()
            maxima.append(peak)
            winners.append(np.flatnonzero(counts==peak) if peak else [])
        return maxima, totals, winners

    observed, totals, winners = measure(labels)
    _, day_starts, day_counts = np.unique(minutes//1440, return_index=True, return_counts=True)
    blocks = [starts[:,None]+np.arange(n) for n in np.unique(day_counts) if n>1
              for starts in [day_starts[day_counts==n]]]
    shuffled = labels.copy()
    rng = np.random.default_rng(65171976)
    null = np.zeros((trials,len(lags)), dtype=int)
    for k in range(trials):
        for block in blocks:
            shuffled[block] = rng.permuted(labels[block], axis=1)
        null[k] = measure(shuffled)[0]
    result = {}
    for i, lag in enumerate(lags):
        result[lag] = dict(distinct_nest_event_pairs=totals[i], max_recurrences=int(observed[i]),
            leading_pairs=[[nests[c//base],nests[c%base]] for c in winners[i]],
            null_median=float(np.median(null[:,i])), null_95_percentile=float(np.quantile(null[:,i],.95)),
            upper_p=float((1+(null[:,i]>=observed[i]).sum())/(trials+1)))
    return dict(rows=len(catalog), classified_nest_events=len(labels), nests=base,
                seed=65171976, trials=trials, lag_minutes=result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('grade_catalog', type=Path)
    parser.add_argument('deep_catalog', type=Path)
    args = parser.parse_args()
    report = dict(eclipses=eclipses(checked(args.grade_catalog, GRADESHA)),
                  deep=pairs(checked(args.deep_catalog, DEEPSHA)))
    print(json.dumps(report, indent=2, allow_nan=False))
