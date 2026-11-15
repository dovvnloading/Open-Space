"""Measure recurring asteroid dimmings using Pál et al. (2020) TSSYS-DR1.

Run with Python 3.10+, numpy and scipy. Supply the original release directory
containing release.merge and lightcurves_spectra/. No network access required.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

TARGETS = (29489, 3177, 2037, 1830, 1803)
# Starting event spacings from the blind search, in hours. The optimization is local.
SEEDS = (12.337, 18.618, 13.128, 12.154, 14.272)
WIDTHS = (0.08, 0.12, 0.16, 0.20, 0.24)
PHASE_BINS = 80


def residuals(root, number, spin_seed):
    path = root / 'lightcurves_spectra' / f'{number}.lc'
    a = np.loadtxt((line for line in path.open() if len(line.split()) == 13),
                   usecols=(1, 2, 3, 4, 5, 6, 7, 9, 10, 11))
    ok = np.all(np.isfinite(a), axis=1) & (a[:, 7] == 0) & (a[:, 4] > 0) & (a[:, 4] < .15)
    a = a[ok]
    if len(a) < 40:
        raise ValueError('too few measurements')
    a = a[a[:, 4] < max(.02, 3 * np.median(a[:, 4]))]
    t = a[:, 0] - a[0, 0]
    x = t / np.ptp(t) - .5
    # Remove distance law, a slowly changing phase function, and primary spin.
    mag = a[:, 3] - 5 * np.log10(a[:, 8] * a[:, 9])
    def fit(f, keep):
        angles = 2*np.pi * t[:, None] * f * np.arange(1, 7)
        design = np.column_stack([np.ones(len(t)), x, x**2, np.sin(angles), np.cos(angles)])
        coefs = np.linalg.lstsq(design[keep], mag[keep], rcond=None)[0]
        return mag - design @ coefs
    keep = np.ones(len(t), dtype=bool)
    f0 = 24 / spin_seed
    for _ in range(2):
        f = minimize_scalar(lambda v: np.sum(fit(v, keep)[keep]**2),
                            bounds=(f0-.45/np.ptp(t), f0+.45/np.ptp(t)), method='bounded').x
        r = fit(f, keep)
        scatter = np.median(np.abs(r[keep] - np.median(r[keep]))) / .67449
        # Remove single-epoch excursions. Dimming half-amplitudes are smaller than the cut.
        keep &= np.abs(r - np.median(r[keep])) < max(.1, 6*scatter)
    r = fit(f, keep)
    return a[keep], r[keep], 24 / f


def search(t, residual, seed, span):
    """Positive-magnitude box search across a local frequency window."""
    scatter = np.median(np.abs(residual - np.median(residual))) / .67449
    y = np.clip(residual, np.median(residual)-4*scatter, np.median(residual)+4*scatter)
    y -= np.mean(y)
    best = (-np.inf, 0, 0, 0)
    widths = np.rint(np.asarray(WIDTHS)*PHASE_BINS).astype(int)
    for f in np.linspace(24/seed-.8/span, 24/seed+.8/span, 161):
        ix = ((t*f % 1)*PHASE_BINS).astype(int)
        n = np.bincount(ix, minlength=PHASE_BINS)
        total = np.bincount(ix, weights=y, minlength=PHASE_BINS)
        c_n = np.r_[0, np.cumsum(np.tile(n, 2))]
        c_y = np.r_[0, np.cumsum(np.tile(total, 2))]
        start = np.arange(PHASE_BINS)
        for w in widths:
            ni = c_n[start+w]-c_n[start]
            sums = c_y[start+w]-c_y[start]
            valid = (ni > 4) & (ni < len(y)-4)
            power = np.full(PHASE_BINS, -np.inf)
            power[valid] = sums[valid] / np.sqrt(ni[valid]*(1-ni[valid]/len(y)))
            k = np.argmax(power)
            if power[k] > best[0]:
                best = (power[k], f, ((k+w/2)/PHASE_BINS) % 1, w/PHASE_BINS)
    _, freq, center, width = best
    return 24/freq, center, width


def event_rows(t, r, period, center, width):
    position = t*24/period-center
    cycles = np.rint(position).astype(int)
    phase = np.abs(position-cycles)
    out = []
    for cycle in np.unique(cycles):
        inside = (cycles == cycle) & (phase < width/2)
        away = (cycles == cycle) & (phase > .27)
        if np.sum(inside) >= 2 and np.sum(away) >= 4:
            depth = np.median(r[inside]) - np.median(r[away])
            out.append((cycle, np.sum(inside), np.sum(away), depth))
    return np.asarray(out, dtype=float).reshape(-1, 4)


def interval(d, rng):
    if not len(d):
        return np.full(3, np.nan)
    boot = np.median(rng.choice(d, (5000, len(d)), replace=True), axis=1)
    return [np.median(d), *np.quantile(boot, [.025, .975])]


def write(path, fields, rows):
    with path.open('w', newline='') as handle:
        w = csv.writer(handle, lineterminator='\n')
        w.writerow(fields.split())
        w.writerows(rows)


def run(root, out, comparisons):
    out.mkdir(parents=True, exist_ok=True)
    metadata = {x[0]: x for x in (line.split() for line in (root/'release.merge').read_text().splitlines())}
    rng = np.random.default_rng(317729489)
    rows, visit_rows, bins, checks, control_rows = [], [], [], [], []
    for i, seed in zip(TARGETS, SEEDS):
        rec = metadata[str(i)]
        a, r, spin = residuals(root, i, float(rec[3]))
        origin = a[0, 0]
        t = a[:, 0] - origin
        span = np.ptp(t)
        p, center, width = search(t, r, seed, span)
        events = event_rows(t, r, p, center, width)
        mid, low, high = interval(events[:, 3], rng)
        even = events[events[:, 0] % 2 == 0, 3]
        odd = events[events[:, 0] % 2 == 1, 3]
        rows.append((i,len(r),origin,span,spin,p,origin+center*p/24,width*p,
                     len(events),sum(events[:, 3]>0),mid,low,high,np.median(even),np.median(odd)))
        for cyc, n, nout, dep in events:
            visit_rows.append((i,int(cyc),origin+(cyc+center)*p/24,int(n),int(nout),dep))
        phase = (t*24/p-center+.5)%1-.5
        for lo in np.linspace(-.5,.5,41)[:-1]:
            use = (phase>=lo) & (phase<lo+.025)
            if np.any(use):
                bins.append((i,lo+.0125,np.sum(use),np.median(r[use]),np.mean(r[use])))
        # Fit the ephemeris using the first time half, then score the untouched flux in the second.
        # The original target and local period range were selected on the entire release.
        first = t < span/2
        a_first = search(t[first],r[first],seed,span/2)
        val = event_rows(t[~first],r[~first],*a_first)
        dval = interval(val[:,3],rng)
        checks.append((i,*a_first,len(val),sum(val[:,3]>0),*dval))
        if comparisons and i in TARGETS[:2]:
            for j, ctrl in metadata.items():
                if j == str(i) or not j.isdigit() or ctrl[6:9] != rec[6:9]:
                    continue
                if int(ctrl[5])<300 or abs(float(ctrl[9])-float(rec[9]))>1.5:
                    continue
                try:
                    b, rb, _ = residuals(root,int(j),float(ctrl[3]))
                except ValueError:
                    continue
                tc = b[:,0]-origin
                if min(tc[-1],span)-max(tc[0],0) < .70*span:
                    continue
                er = event_rows(tc,rb,p,center,width)
                if len(er)<12:
                    continue
                dc = er[:,3]
                control_rows.append((i,j,ctrl[9],len(rb),len(dc),sum(dc>0),np.median(dc)))
    write(out/'measurements.csv','mpc n jd_start span_d spin_h recurrence_h first_center_jd box_width_h events positive median_depth_mag bootstrap_low bootstrap_high even_depth_mag odd_depth_mag',rows)
    write(out/'events.csv','mpc cycle jd_center n_in n_out delta_mag',visit_rows)
    write(out/'fold.csv','mpc phase n median_mag mean_mag',bins)
    write(out/'time_split.csv','mpc train_period_h train_phase train_width_fraction test_events test_positive test_depth_mag test_low test_high',checks)
    if comparisons:
        write(out/'detector_controls.csv','target comparison tess_mag n events positive median_depth_mag',control_rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('release',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--detector-controls',action='store_true')
    args = parser.parse_args()
    run(args.release,args.output,args.detector_controls)
