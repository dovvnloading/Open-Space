"""CHEOPS visit-level single-box search. Inputs contain one extracted source each."""
import hashlib
from pathlib import Path

import numpy as np
from astropy.io import fits

DURATIONS_HOURS = (0.5, 1., 2., 3., 4., 6., 8.)
DEPTHS = (0.0005, 0.0015, 0.005, 0.02)
INJECT_HOURS = (1., 3., 6.)
TRIALS_PER_CELL = 2
SNR_LIMIT = 8.
AUX = {
    'roll': ('ROLL_ANGLE', 'ROLL', 'ROLLANGLE'),
    'background': ('BACKGROUND', 'BG', 'BKG', 'BKG_LC', 'BACKGROUND_LC'),
    'smear': ('SMEARING_LC', 'SMEARING', 'SMEAR', 'SMEAR_LC'),
    'contamination': ('CONTA_LC', 'CONTAM_LC', 'CONTAMINATION', 'CONTAM', 'CONTAMINATION_LC'),
}
FLUX_COLS = ('FLUX', 'PSF_FLUX', 'SAP_FLUX')
TIME_COLS = ('BJD_TIME', 'BJD_TDB', 'BJD', 'BTJD', 'TIME')
QUALITY_COLS = ('EVENT', 'QUALITY', 'FLAG', 'FLAGS', 'SAP_QUALITY')


def _one(names, candidates):
    return next((names[c] for c in candidates if c in names), None)


def _table(hdus):
    for h in hdus:
        if isinstance(h, (fits.BinTableHDU, fits.TableHDU)):
            names = {n.upper(): n for n in h.columns.names}
            if _one(names, TIME_COLS) and _one(names, FLUX_COLS):
                return h, names
    raise ValueError('No scalar photometric light-curve table with recognized time and flux columns')


def _scalar(data, name):
    arr = np.asarray(data[name], dtype=float)
    if arr.ndim != 1:
        raise ValueError('Unsupported multichannel column ' + name)
    return arr


def _coord(header, keys, ra=False):
    value = next((header.get(k) for k in keys if header.get(k) is not None), None)
    if value is None:
        return None
    try:
        v = float(value)
    except (ValueError, TypeError):
        from astropy.coordinates import Angle
        import astropy.units as u
        try:
            v = Angle(str(value), unit=u.hourangle if ra else u.deg).degree
        except ValueError:
            return None
    return float(v) if np.isfinite(v) else None


def metadata(path):
    """Read source identifiers before asking the caller for catalog ephemerides."""
    with fits.open(path, memmap=False) as hdus:
        table, names = _table(hdus)
        head = hdus[0].header.copy()
        head.update(table.header)
        return dict(star=str(head.get('TARGNAME') or head.get('OBJECT') or head.get('TARGET') or ''),
                    ra_deg=_coord(head, ('RA_TARG', 'RA_OBJ', 'RA'), ra=True),
                    dec_deg=_coord(head, ('DEC_TARG', 'DEC_OBJ', 'DEC')),
                    visit=str(head.get('ARCHID') or head.get('OBSID') or head.get('OBS_ID') or ''),
                    product=Path(path).name, subproduct=str(head.get('EXTNAME', '')),
                    columns=sorted(names))


def _read(path):
    with fits.open(path, memmap=False) as hdus:
        table, names = _table(hdus)
        head = hdus[0].header.copy()
        head.update(table.header)
        time_col = _one(names, TIME_COLS)
        t = _scalar(table.data, time_col)
        finite = t[np.isfinite(t)]
        if not len(finite):
            raise ValueError('No finite times')
        col_index = table.columns.names.index(time_col) + 1
        scale = str(head.get(f'TCTYP{col_index}', head.get('TIMESYS', 'TDB'))).strip().upper()
        if scale not in ('TDB', 'BJD_TDB', 'TT'):
            raise ValueError('Unsupported barycentric time scale: ' + scale)
        tunit = str(table.columns[time_col].unit or 'd').strip().lower()
        if tunit not in ('d', 'day', 'days', 'jd', 'bjd'):
            raise ValueError('Time unit is not days: ' + tunit)
        if time_col.upper() == 'TIME' and str(head.get('TIMEREF', head.get('TREFPOS', ''))).upper() not in (
                'SOLARSYSTEM', 'BARYCENTER', 'BARYCENTRIC'):
            raise ValueError('TIME has no positive barycentric reference')
        if np.median(finite) < 2.4e6:
            ref = float(head.get('BJDREF', 0.) or 0.) or (
                float(head.get('BJDREFI', 0.) or 0.) + float(head.get('BJDREFF', 0.) or 0.))
            if not ref and time_col.upper() == 'BTJD':
                ref = 2457000.
            if not ref:
                raise ValueError('Relative BJD without a documented reference')
            t = t + ref
        if scale == 'TT':
            from astropy.time import Time
            valid = np.isfinite(t)
            t = t.copy()
            t[valid] = Time(t[valid], format='jd', scale='tt').tdb.jd
        if not 2457000 < np.nanmedian(t) < 2500000:
            raise ValueError('Times are outside CHEOPS operational epochs')
        flux = _scalar(table.data, _one(names, FLUX_COLS))
        aux = {k: _scalar(table.data, c) for k, aliases in AUX.items() if (c := _one(names, aliases))}
        roll_col = _one(names, AUX['roll'])
        if roll_col and str(table.columns[roll_col].unit or '').strip().lower() in ('rad', 'radian', 'radians'):
            aux['roll'] = np.rad2deg(aux['roll'])
        quality = np.ones(len(t), dtype=bool)
        for col in QUALITY_COLS:
            if col in names:
                quality &= (_scalar(table.data, names[col]) == 0)
        for key, cols in (('xc', ('CENTROID_X', 'XC', 'X_CENTROID')),
                          ('yc', ('CENTROID_Y', 'YC', 'Y_CENTROID'))):
            c = _one(names, cols)
            if c:
                aux[key] = _scalar(table.data, c)
        return t, flux, aux, quality


def _nearest_aux(t, reference):
    rt, _, cols, rq = _read(reference)
    good = np.isfinite(rt) & rq
    ix = np.flatnonzero(good)[np.argsort(rt[good])]
    rt, keep = np.unique(rt[ix], return_index=True)
    ix = ix[keep]
    if len(ix) < 2:
        return {}
    j = np.searchsorted(rt, t).clip(1, len(ix) - 1)
    near = np.where(abs(rt[j] - t) < abs(rt[j - 1] - t), j, j - 1)
    tol = min(2. / 1440., .6 * np.median(np.diff(rt)))
    return {k: np.where(abs(rt[near] - t) <= tol, v[ix[near]], np.nan) for k, v in cols.items()
            if k in AUX}


def _positive(x):
    try:
        f = float(x)
        return f if np.isfinite(f) and f > 0 else None
    except (ValueError, TypeError):
        return None


def _mask(t, ephemerides):
    mask = np.zeros(len(t), bool)
    unresolved = []
    for row in ephemerides:
        t0 = _positive(row.get('epoch_bjd'))
        period = _positive(row.get('period_days'))
        dur = _positive(row.get('duration_days'))
        if t0 is None or dur is None:
            unresolved.append(str(row.get('label', 'unknown')))
            continue
        n = np.rint((t - t0) / period) if period else np.zeros(len(t))
        # Conservative linear propagation, including a buffer for ingress and timing variations.
        e = (_positive(row.get('epoch_error_days')) or 0.) + abs(n) * (
            _positive(row.get('period_error_days')) or 0.)
        mask |= abs(t - (t0 + n * (period or 0.))) <= .65 * dur + 3 * e
    return mask, unresolved


def _mad(a):
    return 1.4826 * np.median(abs(a - np.median(a))) if len(a) else 0.


def _design(t, aux):
    phi = np.deg2rad(aux['roll'])
    cols = [np.ones(len(t)), (t - np.mean(t)) / np.ptp(t)]
    for harmonic in (1, 2, 3, 4):
        cols.extend([np.sin(harmonic * phi), np.cos(harmonic * phi)])
    for key in ('background', 'smear', 'contamination'):
        v = aux[key] - np.median(aux[key])
        s = _mad(v) or np.std(v)
        if s > 0:
            cols.append(v / s)
    return np.stack(cols, axis=1)


def _residual(f, X):
    good = np.ones(len(f), bool)
    for _ in range(4):
        beta = np.linalg.lstsq(X[good], f[good], rcond=1e-10)[0]
        resid = f - X @ beta
        scale = _mad(resid[good])
        if scale <= np.finfo(float).eps:
            break
        new = abs(resid - np.median(resid[good])) < 3.5 * scale
        if new.sum() < max(len(beta) * 3, len(f) // 3) or np.array_equal(new, good):
            break
        good = new
    return resid


def _windows(t, cadence):
    # Before and after means give equal weight to the two flanks despite spacecraft gaps.
    all_windows = []
    for h in DURATIONS_HOURS:
        d = h / 24
        bounds = [np.searchsorted(t, t + k * d) for k in (-1.5, -.5, .5, 1.5)]
        a, b, c, e = bounds
        counts = (b - a, c - b, e - c)
        good = np.ones(len(t), bool)
        for n in counts:
            good &= n >= max(4, .22 * d / cadence)
        good &= t[np.maximum(c - 1, 0)] - t[np.minimum(b, len(t) - 1)] >= .45 * d
        idx = np.flatnonzero(good)
        if len(idx):
            all_windows.append((d, idx, *(x[idx] for x in bounds)))
    return all_windows


def _scan(res, windows, t):
    prefix = np.concatenate(([0.], np.cumsum(res)))
    white = max(_mad(np.diff(res)) / np.sqrt(2), np.finfo(float).eps)
    peaks = []
    for d, idx, a, b, c, e in windows:
        low = (prefix[c] - prefix[b]) / (c - b)
        base = .5 * ((prefix[b] - prefix[a]) / (b - a) + (prefix[e] - prefix[c]) / (e - c))
        depth = base - low
        sample_noise = white * np.sqrt(1. / (c - b) + .25 / (b - a) + .25 / (e - c))
        sigma = np.maximum(sample_noise, _mad(depth))
        snr = depth / sigma
        hits = np.flatnonzero(snr >= SNR_LIMIT)
        # Reduce adjacent trial centers without dropping peaks separated by gaps in time.
        starts = np.flatnonzero((np.diff(hits) > 1) | (np.diff(t[idx[hits]]) > d / 2)) + 1
        for group in np.split(hits, starts):
            if len(group):
                k = group[np.argmax(snr[group])]
                peaks.append((float(snr[k]), float(t[idx[k]]), d, float(depth[k]),
                              float(base[k]), int(a[k]), int(b[k]), int(c[k]), int(e[k])))
    kept = []
    for p in sorted(peaks, reverse=True):
        if all(abs(p[1] - other[1]) >= .5 * (p[2] + other[2]) for other in kept):
            kept.append(p)
    return kept


def _assess(peak, res, t, aux):
    snr, time, dur, dep, base, a, b, c, e = peak
    inside = res[b:c]
    flanks = np.r_[res[a:b], res[c:e]]
    # Re-estimate the mean depth after removing the most influential downward excursion.
    trimmed = (inside.sum() - inside.min()) / (len(inside) - 1)
    coherent = (base - trimmed >= .7 * dep and
                np.count_nonzero(inside < base - .25 * dep) >= max(4, .55 * len(inside)))
    left, right = np.mean(res[a:b]), np.mean(res[c:e])
    stable = min(left, right) - np.mean(inside) > .45 * dep
    # Compare the same 10-degree roll sectors at independently sampled epochs.
    # A full orbit outside the event avoids comparing immediate ingress/egress neighbors.
    phase = (np.mod(aux['roll'], 360) / 10).astype(int).clip(0, 35)
    far = abs(t - time) > .5 * dur + 100. / 1440.
    in_counts = np.bincount(phase[b:c], minlength=36)
    ref_levels, weights = [], []
    for phase_bin in np.flatnonzero(in_counts):
        pts = res[far & (phase == phase_bin)]
        if len(pts) >= 3:
            ref_levels.append(np.median(pts))
            weights.append(in_counts[phase_bin])
    phase_fraction = sum(weights) / len(inside)
    roll_repeat = None
    if phase_fraction >= .7:
        repeat_depth = base - np.average(ref_levels, weights=weights)
        roll_repeat = bool(repeat_depth < .5 * dep)

    # A telemetry or centroid excursion large relative to the visit is a veto.
    # Astrophysical source identification requires the separate image/catalog checks.
    anomalous = []
    for key in ('background', 'smear', 'contamination', 'xc', 'yc'):
        if key not in aux:
            continue
        v = aux[key]
        good = np.isfinite(v)
        if good.sum() == 0:
            continue
        training = v[good]
        flank = np.r_[v[a:b], v[c:e]]
        local = v[b:c]
        if not np.isfinite(flank).any() or not np.isfinite(local).any():
            continue
        noise = max(_mad(training), _mad(flank[np.isfinite(flank)]))
        delta = abs(np.nanmedian(local) - np.nanmedian(flank))
        if delta > 5 * max(noise, np.finfo(float).eps * max(1, abs(np.median(training)))):
            anomalous.append(key)
    checks = dict(multiple_cadences=bool(coherent), both_flanks=bool(stable),
                  roll_independent=roll_repeat, no_large_auxiliary_shift=not anomalous)
    return dict(epoch_bjd=time, duration_days=dur, depth=dep, depth_ppm=dep * 1e6,
                snr=snr, points=int(c - b), roll_reference_fraction=float(phase_fraction),
                measured_left_baseline=float(left), measured_right_baseline=float(right),
                telemetry_excursions=anomalous, checks=checks,
                photometry_pass=all(v is True for v in checks.values()))


def _inject(t, norm_flux, design, windows, aux, exposure_span, base_events, seed):
    # Blindly search every duration/center in the injected raw series and refit trends.
    # The sampled interval includes occultation gaps and points removed by phase masks.
    rng = np.random.default_rng(seed)
    bins = []
    for depth in DEPTHS:
        for h in INJECT_HOURS:
            dur = h / 24.
            recovered = confused = sampled = 0
            for _ in range(TRIALS_PER_CELL):
                center = rng.uniform(*exposure_span)
                offset = abs(t - center)
                # Trapezoid, ingress and egress each 10% of full duration.
                shape = np.clip((.5 - offset / dur) / .1, 0., 1.)
                if (shape > 0).sum() < 4:
                    continue
                sampled += 1
                if any(abs(v['epoch_bjd'] - center) < .5 * (v['duration_days'] + dur)
                       for v in base_events):
                    confused += 1
                    continue
                model = (1. + norm_flux) * (1. - depth * shape) - 1.
                res = _residual(model, design)
                for peak in _scan(res, windows, t):
                    if abs(peak[1] - center) > .35 * dur or not .5 <= peak[2] / dur <= 2:
                        continue
                    if _assess(peak, res, t, aux)['photometry_pass']:
                        recovered += 1
                        break
            bins.append(dict(depth_ppm=int(round(depth * 1e6)), duration_hours=h,
                             trials=TRIALS_PER_CELL, sampled=sampled, confused=confused,
                             recovered=recovered))
    return bins


def analyze(path, ephemerides, seed=0, reference_path=None):
    """Return serializable visit diagnostics, all threshold crossings after merging, and recovery bins."""
    result = dict(status='unresolved_input', events=[], injections=[], limits=[])
    try:
        result.update(metadata(path))
        t, flux, aux, quality = _read(path)
        if reference_path and any(k not in aux for k in AUX):
            aligned = _nearest_aux(t, reference_path)
            result['reference_telemetry'] = [k for k in AUX if k not in aux and k in aligned]
            aux.update({k: aligned[k] for k in result['reference_telemetry']})
        missing = sorted(set(AUX) - aux.keys())
        if missing:
            result['limits'].append('Missing telemetry: ' + ', '.join(missing))
            return result
        good = np.isfinite(t) & np.isfinite(flux) & quality
        for k in AUX:
            good &= np.isfinite(aux[k])
        result['input_rows'] = len(t)
        result['good_rows'] = int(good.sum())
        if good.sum() < 48:
            result['status'] = 'insufficient_sampling'
            return result
        idx = np.flatnonzero(good)[np.argsort(t[good])]
        t, keep = np.unique(t[idx], return_index=True)
        idx = idx[keep]
        flux = flux[idx]
        aux = {k: v[idx] for k, v in aux.items()}
        if len(t) < 48:
            result['status'] = 'insufficient_sampling'
            return result
        span = (float(t[0]), float(t[-1]))
        cadence = float(np.median(np.diff(t)))
        cut, unresolved = _mask(t, ephemerides)
        result['epoch_masks_unresolved'] = unresolved
        result['masked_rows'] = int(cut.sum())
        result['ephemerides_supplied'] = len(ephemerides)
        t, flux = t[~cut], flux[~cut]
        aux = {k: v[~cut] for k, v in aux.items()}
        result.update(rows_searched=len(t), span_bjd=list(span), cadence_days=cadence)
        if len(t) < 48 or np.median(flux) <= 0:
            result['status'] = 'insufficient_sampling'
            return result
        design = _design(t, aux)
        norm_flux = flux / np.median(flux) - 1.
        windows = _windows(t, cadence)
        if not windows:
            result['status'] = 'insufficient_sampling'
            return result
        residuals = _residual(norm_flux, design)
        result['durations_searched_hours'] = [round(24 * w[0], 4) for w in windows]
        result['eligible_boxes'] = sum(len(w[1]) for w in windows)
        result['scatter_ppm'] = float(_mad(residuals) * 1e6)
        result['events'] = [_assess(p, residuals, t, aux) for p in _scan(residuals, windows, t)]
        # Deterministic independent visit/product seeds even when the caller uses seed=0.
        digest = hashlib.sha256((str(seed) + ':' + str(span) + ':' + result['product']).encode())
        derived_seed = int.from_bytes(digest.digest()[:8], 'big')
        result['injection_seed'] = derived_seed
        result['injections'] = _inject(t, norm_flux, design, windows, aux, span,
                                       result['events'], derived_seed)
        result['status'] = 'searched'
    except (OSError, ValueError, TypeError, KeyError, IndexError, np.linalg.LinAlgError) as exc:
        result['limits'].append(type(exc).__name__ + ': ' + str(exc))
    return result



def forced_depth(path, event, ephemerides, reference_path=None):
    """Measure the preselected event in an alternate aperture at fixed epoch and duration."""
    t, flux, aux, quality = _read(path)
    if reference_path and set(AUX)-aux.keys():
        aligned = _nearest_aux(t, reference_path)
        aux.update({k:v for k,v in aligned.items() if k not in aux})
    if set(AUX)-aux.keys():
        return dict(status='incomplete telemetry')
    good = quality & np.isfinite(t) & np.isfinite(flux)
    for k in AUX:
        good &= np.isfinite(aux[k])
    mask,_ = _mask(t, ephemerides)
    good &= ~mask
    ix = np.flatnonzero(good)
    ix = ix[np.argsort(t[ix])]
    if len(ix)<48 or np.median(flux[ix])<=0:
        return dict(status='insufficient coverage')
    _, unique = np.unique(t[ix],return_index=True)
    ix = ix[unique]
    if len(ix)<48:
        return dict(status='insufficient coverage')
    t,flux = t[ix],flux[ix]
    aux = {k:v[ix] for k,v in aux.items()}
    res = _residual(flux/np.median(flux)-1., _design(t, aux))
    d, t0 = event['duration_days'],event['epoch_bjd']
    a,b,c,e = np.searchsorted(t, t0+np.array([-1.5,-.5,.5,1.5])*d)
    if min(b-a,c-b,e-c)<4:
        return dict(status='insufficient fixed-window samples')
    dep = .5*(np.mean(res[a:b])+np.mean(res[c:e]))-np.mean(res[b:c])
    ratio=float(dep/event['depth'])
    return dict(status='agreement' if .5<=ratio<=2. else 'depth discrepancy',
                depth=float(dep),points=int(c-b),depth_ratio=ratio)
