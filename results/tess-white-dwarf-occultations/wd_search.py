"""Native-cadence, event-seeded periodic box and individual dip search for a single visit."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import binned_statistic


DURATIONS_MIN = (2., 4., 8., 16., 32., 64., 120.)
P_MIN, P_MAX = 0.08, 7.0
SEEDS, REFINE = 18, 16


def mad(x):
    return 1.4826 * np.median(abs(x - np.median(x)))


def flatten(time, flux):
    """Divide by an interpolated 6-hour median in each uninterrupted segment."""
    baseline = np.full_like(flux, np.nan, dtype=float)
    for ix in np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > 0.12) + 1):
        if len(ix) < 20 or np.ptp(time[ix]) < 0.1:
            continue
        t = time[ix] - time[ix[0]]
        edges = np.arange(0, t[-1] + .25, .25)
        level = binned_statistic(t, flux[ix], 'median', bins=edges).statistic
        center = binned_statistic(t, t, 'median', bins=edges).statistic
        ok = np.isfinite(level) & (level > 0)
        if ok.any():
            baseline[ix] = np.interp(t, center[ok], level[ok])
    return baseline


def box_scores(t, y, dt, width, noise):
    """Evaluate a central box against two flanks per time, including coverage cuts."""
    idx = [np.searchsorted(t, t + d*width, side=side) for d, side in
           ((-1.8, 'left'), (-.8, 'right'), (-.5, 'left'), (.5, 'right'),
            (.8, 'left'), (1.8, 'right'))]
    sums = np.r_[0., np.cumsum(y)]
    l, m, a, b, n, r = idx
    nin, nleft, nright = b-a, m-l, r-n
    required = max(1, int(np.ceil(.58 * width / dt - 1e-5)))
    ok = (nin >= required) & (nleft >= required) & (nright >= required)
    ns = np.maximum(1, nleft+nright)
    depth = (sums[r]-sums[n]+sums[m]-sums[l])/ns - (sums[b]-sums[a])/np.maximum(1, nin)
    z = depth / (noise * np.sqrt(1 / ns + 1 / np.maximum(nin, 1)))
    if ok.sum() > 10:
        z /= max(1., mad(z[ok]))
    z[~ok] = np.nan
    return depth, z


def event_seeds(t, z, width):
    """Non-overlapping local maxima for candidate periods; isolated hits are uncapped."""
    # find_peaks excludes valid end points next to masked gaps unless missing scores are -inf.
    from scipy.signal import find_peaks
    rank = find_peaks(np.r_[-np.inf, np.nan_to_num(z, nan=-np.inf), -np.inf],
                      height=2.5)[0] - 1
    rank = rank[np.argsort(-z[rank], kind='stable')]
    selected = []
    for idx in rank:
        if len(selected) >= SEEDS and z[idx] < 7.:
            break
        if all(abs(t[idx] - t[j]) >= width for j in selected):
            selected.append(idx)
    return np.array(selected, dtype=int)


def period_seeds(epochs, strength, width, span):
    """Use aligned triples of single-event seeds; then refine only the strongest votes."""
    if len(epochs) < 3:
        return []
    ix, jx = np.triu_indices(len(epochs), 1)
    delta = abs(epochs[jx]-epochs[ix])
    lo = np.maximum(1, np.ceil(delta/P_MAX).astype(int))
    hi = np.floor(delta/max(P_MIN, 4.0*width)).astype(int)
    sizes = np.maximum(0, hi-lo+1)
    if not sizes.any():
        return []
    pair = np.repeat(np.arange(len(lo)), sizes)
    # Each possible cycle count between the anchoring events defines a blind ephemeris.
    cycles = np.arange(len(pair)) - np.repeat(np.cumsum(sizes)-sizes, sizes) + lo[pair]
    periods = delta[pair] / cycles
    origin = epochs[ix[pair]]
    votes = np.zeros(len(pair))
    count = np.zeros(len(pair), dtype=int)
    for time, power in zip(epochs, strength):
        d = (time-origin) / periods
        aligned = abs(d-np.rint(d)) * periods <= .7*width
        votes += aligned * min(power, 12.)
        count += aligned
    eligible = np.flatnonzero(count >= 3)
    # Penalize dense ephemerides that place very many empty transits between seeds.
    merit = votes / np.sqrt(np.maximum(3., span/periods))
    ranked = eligible[np.argsort(-merit[eligible], kind='stable')]
    seeds = []
    for k in ranked:
        p, e = periods[k], origin[k]
        if any(abs(p-q)*span/p <= .4*width and
               abs((e-c+.5*q) % q-.5*q) <= width for q, c in seeds):
            continue
        seeds.append((p, e))
        if len(seeds) >= REFINE:
            break
    return seeds


def nearest(t, centers):
    pos = np.searchsorted(t, centers)
    a = np.clip(pos-1, 0, len(t)-1)
    b = np.clip(pos, 0, len(t)-1)
    return np.where(abs(t[a]-centers) <= abs(t[b]-centers), a, b)


def periodic_candidates(t, depth, z, dt, width, seeds):
    if not seeds:
        return []
    rows = []
    for period, epoch in seeds:
        # A small period change shifts either end of the visit by <= 1.2 box widths.
        ps = period * (1 + np.linspace(-1.2, 1.2, 7) * width / np.ptp(t))
        ps = ps[(ps >= P_MIN) & (ps <= P_MAX)]
        if not len(ps):
            continue
        periods = np.tile(ps, 5)[:, None]
        origins = np.repeat(epoch + np.linspace(-.6, .6, 5)*width, len(ps))[:, None]
        k = np.arange(np.floor((t[0]-epoch)/ps.min())-2, np.ceil((t[-1]-epoch)/ps.min())+3)
        centers = origins + periods*k
        col = nearest(t, centers)
        vals = z[col]
        good = np.isfinite(vals) & (abs(t[col]-centers) <= .53*dt)
        n = good.sum(axis=1)
        sig = ((vals >= 2.5) & good).sum(axis=1)
        values = np.where(good, vals, 0.)
        score = values.sum(axis=1)/np.sqrt(np.maximum(1, n))
        # Demand recurrence after removing the single deepest event.
        other = (values.sum(axis=1)-values.max(axis=1)) / np.sqrt(np.maximum(1, n-1))
        valid = (n >= 3) & (sig >= 3) & (score >= 7.) & (other >= 4.)
        if not valid.any():
            continue
        winner = np.argmax(np.where(valid, score, -np.inf))
        p, e = periods[winner,0], origins[winner,0]
        rows.append(dict(period_days=float(p), epoch_btjd=float(t[0]+(e-t[0])%p),
                         duration_min=float(width*1440), score=float(score[winner]),
                         events_covered=int(n[winner]), events_above_2p5=int(sig[winner]),
                         depth=float(np.mean(depth[col[winner,good[winner]]]))))
    return rows


def search(time_btjd, aperture_flux, cadence_s):
    """Blind search on one quality-filtered visit. Scores are descriptive, not FAPs."""
    t, f = np.asarray(time_btjd, float), np.asarray(aperture_flux, float)
    if t.ndim != 1 or t.shape != f.shape or not np.isfinite(cadence_s) or cadence_s <= 0:
        raise ValueError('Supply matched vectors and a positive cadence in seconds')
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    t, idx = np.unique(t, return_index=True)
    f = f[idx]
    dt = cadence_s / 86400.
    report = {'periodic': [], 'isolated': [], 'usable_cadences': 0}
    if len(t) < 20:
        return report
    base = flatten(t, f)
    ok = np.isfinite(base)
    t, y = t[ok], f[ok] / base[ok] - 1.
    report['usable_cadences'] = len(t)
    if len(t) < 20 or np.ptp(t) < .1:
        return report
    noise = mad(y)
    if not np.isfinite(noise) or noise <= 0:
        return report
    rows, hits = [], []
    for width in np.unique(np.maximum(dt, np.array(DURATIONS_MIN) / 1440.)):
        depth, z = box_scores(t, y, dt, width, noise)
        ix = event_seeds(t, z, width)
        for j in ix[z[ix] >= 7.]:
            hits.append(dict(epoch_btjd=float(t[j]), duration_min=float(width*1440),
                             depth=float(depth[j]), score=float(z[j])))
        use = ix[:SEEDS]
        candidates = period_seeds(t[use], z[use], width, np.ptp(t))
        rows.extend(periodic_candidates(t, depth, z, dt, width, candidates))
    for row in sorted(rows, key=lambda x:-x['score']):
        p, e = row['period_days'], row['epoch_btjd']
        width = row['duration_min']/1440.
        if any(abs(p-r['period_days'])*np.ptp(t)/p < width and
               abs((e-r['epoch_btjd']+.5*p)%p-.5*p) < width for r in report['periodic']):
            continue
        report['periodic'].append(row)
    for row in sorted(hits, key=lambda x:-x['score']):
        if all(abs(row['epoch_btjd']-r['epoch_btjd']) * 1440 >
               .5 * (row['duration_min']+r['duration_min']) for r in report['isolated']):
            report['isolated'].append(row)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('curve', type=Path, help='NPZ with time and flux, optional quality')
    parser.add_argument('--cadence', type=float, required=True, help='seconds')
    args = parser.parse_args()
    with np.load(args.curve) as data:
        ok = data['quality'] == 0 if 'quality' in data else slice(None)
        print(json.dumps(search(data['time'][ok], data['flux'][ok], args.cadence), allow_nan=False))
