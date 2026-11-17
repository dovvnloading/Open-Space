"""Explore event morphology, recurrence, apertures, background and detector peers.

Descriptive diagnostics only. Peer selection uses the full available sample on
that camera/CCD. Correlations compare local residuals after linear trend removal.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lombscargle

from scan import BIN, Binning, read_curve, triangle, write_csv


def read_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def average(bins, curve, data):
    valid = curve.good & np.isfinite(curve.time) & np.isfinite(curve.flux)
    a = data[valid]
    good = np.isfinite(a)
    count = np.bincount(bins.bin_ids[good], minlength=len(bins.count))
    sums = np.bincount(bins.bin_ids[good], weights=a[good], minlength=len(count))
    out = np.full(len(count), np.nan)
    use = bins.valid & (count > 0)
    out[use] = sums[use] / count[use]
    return out


def residual(x, t):
    design = np.column_stack([np.ones(len(t)), t - np.mean(t)])
    return x - design @ np.linalg.lstsq(design, x, rcond=None)[0]


def profile(t, y, peak, duration, fraction):
    shape = triangle(t, peak, duration, fraction)
    w = np.isfinite(y) & (np.abs(t - peak) < 1.5 * duration)
    t, y, shape = t[w], y[w], shape[w]
    if len(t) < 8 or (shape > 0).sum() < 3:
        return None
    r = residual(shape, t)
    if r @ r < 1e-6:
        return None
    return float(-np.dot(r, y) / np.dot(r, r))


def correlation(t, y, other, loc):
    if other is None:
        return None
    good = loc & np.isfinite(y) & np.isfinite(other)
    if good.sum() < 8:
        return None
    a, b = residual(y[good], t[good]), residual(other[good], t[good])
    if min(np.std(a), np.std(b)) < 1e-10:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def shape_fit(t, flux, peak, duration, fraction):
    t = t - peak
    valid = np.isfinite(flux) & (np.abs(t) < duration * 1.5)
    t, y = t[valid], flux[valid]
    if len(y) < 12:
        return {}
    sigma = max(.0001, np.median(np.abs(np.diff(y))) / .954)
    def err(par, sym):
        base, slope, depth, center, width = par[:5]
        f = .5 if sym else par[5]
        return (base + slope*t - depth*triangle(t, center, width, f) - y) / sigma
    # The winning search template supplies a starting time, not a measured skew.
    args = ([np.median(y), 0, .02, 0, duration],
            [-np.inf,-np.inf,0,-duration*.6,duration*.45],
            [np.inf,np.inf,np.inf,duration*.6,duration*1.8])
    def fit(sym):
        starts = [.5] if sym else [.2, .5, .8]
        opts = []
        for f in starts:
            x, lo, hi = [list(a) for a in args]
            if not sym:
                x.append(f)
                lo.append(.03)
                hi.append(.97)
            opts.append(least_squares(err, x, args=(sym,), bounds=(lo,hi), loss='soft_l1'))
        return min(opts, key=lambda opt: opt.cost)
    sym = fit(True)
    skew = fit(False)
    p = skew.x
    return dict(fit_peak=peak+p[3],fit_duration_h=24*p[4],fit_ingress_fraction=p[5],
                fit_depth=p[2],fit_baseline_at_peak=p[0]+p[1]*p[3],
                symmetric_loss=sym.cost,asymmetric_loss=skew.cost)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--triggers', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    stars = {}
    for row in read_rows(args.manifest):
        curve = read_curve(args.manifest.parent/row['path'], row)
        b = Binning(curve)
        flux = b.bin_flux(b.flux)
        aux = {k: average(b, curve, x) for k,x in curve.diagnostics.items()}
        stars[(row['tic'], row['sector'])] = (curve.meta, b.centers, flux, aux)
    out = vet_rows(stars, read_rows(args.triggers))
    write_csv(args.out, out, sorted(set().union(*(row.keys() for row in out))) if out else ['tic','sector','peak_btjd'])


def vet_rows(stars, triggers):
    out = []
    for event in triggers:
        star = event['tic'], event['sector']
        meta, t, y, aux = stars[star]
        peak, d, f = (float(event[k]) for k in ('peak_btjd', 'duration_h','ingress_fraction'))
        d /= 24
        peak_loc = np.abs(t - peak) < d*1.5
        row = dict(tic=star[0], sector=star[1], peak_btjd=peak)
        row.update(shape_fit(t, y, peak,d,f))
        # Values as stored in the QLP product. Flux columns are normalized.
        for k in ('SAP_FLUX','SAP_BKG','SAP_X','SAP_Y','DET_FLUX','DET_FLUX_SML','DET_FLUX_LAG'):
            if k in aux:
                row[k + '_r'] = correlation(t, y, aux[k], peak_loc)
                row[k + '_amplitude'] = profile(t, aux[k], peak,d,f)
        rlist, deg = [], []
        for key,(m,tx,yx,_) in stars.items():
            if key[1] == star[1] and key != star and meta['camera'] != '' and meta['ccd'] != '' and (
                    meta['camera'],meta['ccd']) == (m['camera'],m['ccd']):
                aligned = np.interp(t, tx, yx, left=np.nan, right=np.nan)
                cc = correlation(t, y, aligned, peak_loc)
                if cc is not None:
                    rlist.append(cc)
                    dra = (meta['ra_obj'] - m['ra_obj'])*np.cos(np.deg2rad(meta['dec_obj']))
                    deg.append(np.hypot(dra, meta['dec_obj']-m['dec_obj']))
        row.update(peer_n=len(rlist),peer_r_gt_07=int(np.sum(np.array(rlist) > .7)),
                   peer_r_median=np.median(rlist) if rlist else None,
                   peer_r_max=max(rlist,default=None),nearest_peer_deg=min(deg,default=None))
        # A sine-wave peak is a variability diagnostic. Window function and
        # multiple testing preclude using it as an automatic periodicity veto.
        reg = aux.get('DET_FLUX',y)
        valid = np.isfinite(reg)
        if valid.sum() > 40:
            x = t[valid]-t[valid][0]
            values = reg[valid] - np.median(reg[valid])
            freq = np.arange(2./np.ptp(x), 12, 1./(6*np.ptp(x)))
            power = lombscargle(x, values, 2*np.pi*freq, normalize=True)
            peak_i = np.argmax(power)
            row.update(ls_period_d=1/freq[peak_i],ls_normalized_power=power[peak_i])
        out.append(row)
    return out


if __name__ == '__main__':
    main()
