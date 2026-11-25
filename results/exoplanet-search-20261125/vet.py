#!/usr/bin/env python3
"""Fixed-ephemeris checks on unfiltered native TESS samples.

This is a follow-up measurement, not a period search, physical transit fit,
false-alarm probability or planet-probability calculation.  Each predicted event
is measured against robust, local linear flanks; no in-transit value is clipped.
Times are BTJD (BJD_TDB - 2457000), durations and periods are days. A positive
photometric depth means dimming, in ppm of the local fitted flux. Instrumental
changes use the same deficit convention in their original units.

Nominal eligibility: QUALITY=0, finite samples, >=3 in-event and each flank,
>=60% expected in-event coverage and >=40% in each flank. Expected counts use
that sector's median native cadence. Flanks lie 0.8--2.5 full durations from
center, capped at 0.45 period to avoid adjacent events. Quoted formal errors
include the local baseline fit. A descriptive beta factor uses duration/2
residual bins in the flanks (floor 1). Neither beta-adjusted errors nor off-phase
ranks are search-calibrated significances. SAP and PDCSAP share photons.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

FIELDS = ('time', 'pdcsap', 'pdcsap_err', 'sap', 'sap_err', 'quality', 'bkg', 'cx', 'cy')


def mad(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826022185 * np.median(np.abs(x - np.median(x)))) if x.size else float('nan')


def cadence(time):
    ds = np.diff(np.unique(np.asarray(time)[np.isfinite(time)]))
    ds = ds[ds > 0]
    if not len(ds):
        return float('nan')
    mid = np.median(ds)
    return float(np.median(ds[ds < 3 * mid]))


def _fit_baseline(x, y, e):
    """IRLS baseline, downweighting flanks only; errors include excess scatter."""
    design = np.column_stack((np.ones(len(x)), x))
    initial = np.linalg.lstsq(design, y, rcond=None)[0]
    sigma = mad(y - design @ initial)
    good = np.isfinite(e) & (e > 0)
    fallback = float(np.median(e[good])) if np.any(good) else sigma
    if not np.isfinite(fallback) or fallback <= 0:
        fallback = max(float(np.max(np.abs(y))) * 1e-12, 1e-12)
    ee = np.where(good, e, fallback)
    weight = 1 / ee**2
    robust = np.ones(len(x))
    for _ in range(4):
        normal = design.T @ ((weight * robust)[:, None] * design)
        cov = np.linalg.pinv(normal)
        coef = cov @ (design.T @ (weight * robust * y))
        residual = y - design @ coef
        sigma = max(mad(residual), float(np.median(ee)))
        robust = np.minimum(1.0, 3 * sigma / np.maximum(np.abs(residual), sigma * 1e-12))
    inflation = max(1.0, sigma / float(np.median(ee)))
    cov = cov * inflation**2
    return coef, cov, residual, sigma, inflation, ee


def _beta(x, residual, sigma, width, dt):
    # Distinct time bins preserve the actual gaps and time ordering.
    bins = np.floor(x / width).astype(np.int64)
    means, ns = [], []
    for b in np.unique(bins):
        values = residual[bins == b]
        if len(values) >= max(3, 0.4 * width / dt):
            means.append(float(np.mean(values)))
            ns.append(len(values))
    if len(means) < 4 or sigma <= 0:
        return 1.0, len(means)
    expected = sigma * math.sqrt(float(np.mean(1 / np.asarray(ns, float))))
    return max(1.0, float(np.std(means, ddof=1)) / expected), len(means)


def measure_event(data, center, duration, period, field='pdcsap',
                  inner=0.8, outer=2.5, paired=False, normalize=True):
    """Return all coverage outcomes; statistical numbers only when eligible."""
    t = np.asarray(data['time'], float)
    y = np.asarray(data[field], float)
    q = np.asarray(data['quality'])
    x = t - center
    dt = cadence(t)
    extent = min(outer * duration, 0.45 * period)
    result = dict(center_btjd=float(center), duration_days=float(duration),
                  field=field, eligible=False, native_cadence_days=dt,
                  flank_inner_durations=inner, flank_outer_days=extent)
    if not np.isfinite(dt) or dt <= 0 or extent <= inner * duration:
        result['reason'] = 'invalid_cadence_or_overlapping_flanks'
        return result
    ins = np.isfinite(t) & (np.abs(x) <= 0.5 * duration)
    left = np.isfinite(t) & (x >= -extent) & (x <= -inner * duration)
    right = np.isfinite(t) & (x >= inner * duration) & (x <= extent)
    valid = np.isfinite(t) & np.isfinite(y) & (q == 0)
    if paired:
        valid &= np.isfinite(data['sap']) & np.isfinite(data['pdcsap'])
    n_raw = int(np.sum(ins))
    result['raw_in_event_samples'] = n_raw
    result['flagged_in_event_samples'] = int(np.sum(ins & (q != 0)))
    i, l, r = ins & valid, left & valid, right & valid
    ns = [int(np.sum(z)) for z in (i, l, r)]
    expected = [duration / dt, (extent - inner * duration) / dt, (extent - inner * duration) / dt]
    result['counts'] = dict(zip(('in', 'left', 'right'), ns))
    result['coverage'] = dict(zip(('in', 'left', 'right'), (min(1.0, n / e) for n, e in zip(ns, expected))))
    if any(n < 3 for n in ns) or any(n < frac * e for n, frac, e in zip(ns, (.6, .4, .4), expected)):
        result['reason'] = 'insufficient_event_or_two_sided_baseline_coverage'
        return result
    e = np.asarray(data.get(field + '_err', np.full(len(t), np.nan)), float)
    sides = l | r
    coef, cov, residual, sigma, inflation, ee = _fit_baseline(x[sides], y[sides], e[sides])
    level = float(coef[0])
    scale = level if normalize else 1.0
    if not np.isfinite(level) or (normalize and level <= 0):
        result['reason'] = 'nonpositive_photometric_baseline'
        return result
    mean_x = float(np.mean(x[i]))
    vector = np.array([1.0, mean_x])
    difference = float(vector @ coef - np.mean(y[i]))
    ei = e[i]
    baseline_e = float(np.median(ee))
    ei = np.where(np.isfinite(ei) & (ei > 0), ei, baseline_e) * inflation
    variance = float(np.sum(ei**2) / ns[0]**2 + vector @ cov @ vector)
    formal = math.sqrt(max(variance, 0.0)) / abs(scale)
    beta, nblocks = _beta(x[sides], residual, sigma, duration / 2, dt)
    depth = difference / scale
    multiplier = 1e6 if normalize else 1.0
    result.update(eligible=True, reason='measured', depth=depth * multiplier,
                  formal_error=formal * multiplier, red_error=formal * beta * multiplier,
                  snr_formal=depth / formal if formal > 0 else None,
                  snr_red=depth / (formal * beta) if formal > 0 else None,
                  beta=beta, flank_residual_blocks=nblocks, baseline=level,
                  baseline_slope_per_day=float(coef[1]),
                  flank_scatter=sigma / abs(scale) * multiplier,
                  uncertainty_source='pipeline_errors_with_flank_scatter_floor' if np.any(np.isfinite(e[sides]) & (e[sides] > 0)) else 'empirical_flank_scatter',
                  units='ppm' if normalize else ('pixel' if field in ('cx', 'cy') else 'native_flux'))
    return result


def aggregate(events):
    valid = [e for e in events if e.get('eligible') and e.get('red_error', 0) > 0]
    result = {'opportunities': len(events), 'covered_events': len(valid)}
    if not valid:
        return result
    d = np.asarray([e['depth'] for e in valid])
    er = np.asarray([e['red_error'] for e in valid])
    n = len(d)
    eq = float(np.mean(d))
    propagated = float(np.sqrt(np.sum(er**2)) / n)
    between = float(np.std(d, ddof=1) / np.sqrt(n)) if n > 1 else None
    eqerror = max(propagated, between or 0)
    w = er**-2
    weighted = float(np.average(d, weights=w))
    chi2 = float(np.sum(((d - weighted) / er)**2))
    we = float(np.sqrt(max(1, chi2 / max(1, n-1)) / np.sum(w)))
    result.update(equal_event_depth=eq, equal_event_error=eqerror,
                  equal_event_snr=eq / eqerror if eqerror else None,
                  propagated_equal_event_error=propagated,
                  between_event_standard_error=between,
                  weighted_depth=weighted, weighted_error=we,
                  weighted_snr=weighted / we if we else None,
                  positive_events=int(np.sum(d > 0)), positive_events_over_2_red_sigma=int(np.sum(d / er > 2)),
                  depth_consistency_chi2=chi2, depth_consistency_dof=n-1,
                  median_event_beta=float(np.median([e['beta'] for e in valid])),
                  units=valid[0]['units'])
    if n > 1:
        dropped = []
        for j in range(n):
            dd, ee = np.delete(d, j), np.delete(er, j)
            unc = max(float(np.sqrt(np.sum(ee**2)) / len(dd)), float(np.std(dd, ddof=1) / np.sqrt(len(dd))) if len(dd) > 1 else 0)
            dropped.append(float(np.mean(dd) / unc))
        result['drop_one_min_equal_event_snr'] = min(dropped)
        result['maximum_absolute_event_fraction_of_sum'] = float(np.max(np.abs(d)) / abs(np.sum(d))) if np.sum(d) else None
    return result


def all_events(data, period, epoch, duration, field='pdcsap', inner=.8, outer=2.5, paired=False):
    t = np.asarray(data['time'], float)
    good = t[np.isfinite(t)]
    if not len(good):
        return []
    first = math.ceil((float(np.min(good)) - duration / 2 - epoch) / period)
    last = math.floor((float(np.max(good)) + duration / 2 - epoch) / period)
    rows = []
    for k in range(first, last + 1):
        event = measure_event(data, epoch + k * period, duration, period, field, inner, outer, paired)
        event['cycle'] = k
        rows.append(event)
    return rows


def odd_even(events):
    a = aggregate([e for e in events if e['cycle'] % 2 == 0])
    b = aggregate([e for e in events if e['cycle'] % 2])
    out = {'even': a, 'odd': b}
    if a['covered_events'] and b['covered_events']:
        error = math.hypot(a['equal_event_error'], b['equal_event_error'])
        out['even_minus_odd_depth'] = a['equal_event_depth'] - b['equal_event_depth']
        out['difference_error'] = error
        out['difference_snr'] = out['even_minus_odd_depth'] / error if error else None
    return out


def phase_controls(data, period, epoch, duration, actual, field, nphase=96):
    # Phases are diagnostic alternatives, not independent realizations.
    grid = np.arange(1, nphase) / nphase
    guard = min(2.5 * duration, .45 * period) / period + .5 * duration / period
    grid = grid[np.minimum(grid, 1-grid) > guard]
    controls = []
    for phase in grid:
        summary = aggregate(all_events(data, period, epoch + phase * period, duration, field))
        comparable = summary['covered_events'] >= max(2, math.ceil(.8 * actual['covered_events']))
        controls.append({'phase': float(phase), 'comparable_coverage': comparable, **summary})
    comparable = [c for c in controls if c['comparable_coverage'] and 'equal_event_depth' in c]
    out = {'tested': len(controls), 'comparable': len(comparable), 'controls': controls,
           'interpretation': 'conditional off-phase rank; correlated and not a full-search or survey false-alarm probability'}
    if comparable and 'equal_event_depth' in actual:
        depths = np.asarray([c['equal_event_depth'] for c in comparable])
        scale = mad(depths)
        out.update(as_deep_or_deeper=int(np.sum(depths >= actual['equal_event_depth'])),
                   median_control_depth=float(np.median(depths)), control_depth_mad=scale,
                   depth_minus_control_median_in_mad=(actual['equal_event_depth'] - float(np.median(depths))) / scale if scale > 0 else None)
    return out


def _clean(value):
    if isinstance(value, dict):
        return {k: _clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_clean(x) for x in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def vet(cache, tic, period, epoch, duration, nphase=96):
    if not (np.isfinite(period) and np.isfinite(epoch) and np.isfinite(duration) and 0 < duration < .25 * period):
        raise ValueError('Require finite ephemeris with 0 < duration < period/4.')
    files = sorted(Path(cache).glob(f'tic{tic}_s*.npz'), key=lambda p: int(re.search(r'_s(\d+)', p.name).group(1)))
    if not files:
        raise FileNotFoundError(f'No tic{tic}_s*.npz files in {cache}')
    result = {'tic': int(tic), 'period_days': period, 'epoch_btjd': epoch, 'duration_days': duration,
              'units': 'BTJD; days; photometric depths in ppm, positive for dimming',
              'status': 'diagnostics_only', 'sectors': [],
              'limitations': ['Fixed-ephemeris follow-up chosen after a search; no calibrated detection FAP.',
                             'SAP/PDCSAP and aperture variants share photons.',
                             'No physical transit fit, source localization, contamination model, or planetary FPP.',
                             'Local beta estimates have few bins; equal-event scatter and off-phase behavior must also be inspected.',
                             'No ephemeris optimization in other sectors; phase drift can cause nondetections.']}
    total = {f: [] for f in ('pdcsap', 'sap')}
    total_secondary = {f: [] for f in total}
    for path in files:
        sector = int(re.search(r'_s(\d+)', path.name).group(1))
        with np.load(path, allow_pickle=False) as z:
            missing = [key for key in ('time', 'quality', 'pdcsap', 'sap') if key not in z]
            if missing:
                raise ValueError(f'{path.name}: missing {missing}')
            data = {k: np.array(z[k], copy=True) if k in z else np.full(len(z['time']), np.nan) for k in FIELDS}
        if any(np.ndim(x) != 1 or len(x) != len(data['time']) for x in data.values()):
            raise ValueError(f'{path.name}: mismatched or nonscalar columns')
        native_t = data['time'][np.isfinite(data['time'])]
        midpoint = float((np.min(native_t) + np.max(native_t)) / 2) if len(native_t) else float('nan')
        record = {'sector': sector, 'npz_file': path.name, 'npz_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'products': {}}
        metadata = path.with_suffix('.json')
        if metadata.exists():
            record['metadata'] = json.loads(metadata.read_text())
            record['metadata_sha256'] = hashlib.sha256(metadata.read_bytes()).hexdigest()
        for field in total:
            events = all_events(data, period, epoch, duration, field)
            secondary = all_events(data, period, epoch + .5 * period, duration, field)
            for e in events + secondary:
                e['sector'] = sector
            summary = aggregate(events)
            rec = {'summary': summary, 'events': events, 'odd_even': odd_even(events),
                   'first_half': aggregate([e for e in events if e['center_btjd'] < midpoint]),
                   'second_half': aggregate([e for e in events if e['center_btjd'] >= midpoint]),
                   'secondary_half_phase': {'summary': aggregate(secondary), 'events': secondary},
                   'paired_cadences': aggregate(all_events(data, period, epoch, duration, field, paired=True)),
                   'sensitivity': {}}
            for label, factor, outer in [('narrow_flanks', 1, 1.8), ('wide_flanks', 1, 3.5), ('short_duration', .75, 2.5), ('long_duration', 1.25, 2.5)]:
                rec['sensitivity'][label] = aggregate(all_events(data, period, epoch, duration * factor, field, outer=outer))
            if nphase > 0:
                rec['phase_controls'] = phase_controls(data, period, epoch, duration, summary, field, nphase)
            record['products'][field] = rec
            total[field].extend(events)
            total_secondary[field].extend(secondary)
        instrument = []
        for e in record['products']['sap']['events']:
            row = {'center_btjd': e['center_btjd'], 'cycle': e['cycle']}
            for field in ('bkg', 'cx', 'cy'):
                row[field] = measure_event(data, e['center_btjd'], duration, period, field=field, normalize=False)
            instrument.append(row)
        record['instrument_diagnostics'] = instrument
        result['sectors'].append(record)
    result['combined'] = {f: {'summary': aggregate(total[f]), 'odd_even': odd_even(total[f]), 'secondary_half_phase': aggregate(total_secondary[f])} for f in total}
    return _clean(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--tic', type=int, required=True)
    parser.add_argument('--period', type=float, required=True)
    parser.add_argument('--epoch', type=float, required=True)
    parser.add_argument('--duration', type=float, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--phase-grid', type=int, default=96, help='Diagnostic phase grid; 0 disables phase controls.')
    args = parser.parse_args()
    if args.phase_grid not in (0,) and args.phase_grid < 8:
        parser.error('--phase-grid must be 0 or at least 8')
    report = vet(args.cache, args.tic, args.period, args.epoch, args.duration, args.phase_grid)
    args.out.mkdir(parents=True, exist_ok=True)
    output = args.out / f'tic{args.tic}_vet.json'
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'output': str(output), 'combined': report['combined']}, allow_nan=False))


if __name__ == '__main__':
    main()
