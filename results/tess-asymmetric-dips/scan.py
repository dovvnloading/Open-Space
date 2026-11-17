"""Search FFI SAP photometry and measure event-level selection on the same curves."""
import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import chi2


BIN = 1 / 48  # days
DURATIONS = (6, 12, 24)  # total duration in hours
DEPTHS = (.01, .02, .05)
FRACTIONS = (.2, .5, .8)  # time to minimum as fraction of total duration
MIN_DEPTH = .01
MIN_SCORE = 7


def mad(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def triangle(t, peak, duration, fraction):
    dt = (t - peak) / duration
    return np.maximum(0, np.minimum(1 + dt / fraction, 1 - dt / (1 - fraction)))


@lru_cache(maxsize=None)
def templates():
    out = []
    for hours in DURATIONS:
        n = round(hours / (24 * BIN))
        for fraction in FRACTIONS:
            # One half-duration of shoulder on either side of the entire event.
            idx = np.arange(-math.ceil((fraction + .5) * n),
                            math.ceil((1.5 - fraction) * n) + 1)
            signal = triangle(idx, 0, n, fraction)
            design = np.column_stack((np.ones(len(idx)), idx))
            residual = signal - design @ np.linalg.lstsq(design, signal, rcond=None)[0]
            kernel = -residual / np.dot(residual, residual)
            out.append((hours, fraction, idx, kernel))
    return out


@dataclass
class Curve:
    time: np.ndarray
    flux: np.ndarray
    good: np.ndarray
    diagnostics: dict
    meta: dict


AUX_COLUMNS = ('SAP_BKG', 'BKG', 'DET_FLUX', 'DET_FLUX_SML', 'DET_FLUX_LAG',
               'KSPSAP_FLUX', 'PDCSAP_FLUX', 'MOM_CENTR1', 'MOM_CENTR2',
               'X_CENTROID', 'Y_CENTROID', 'SAP_X', 'SAP_Y', 'SYS_RM_FLUX')


def read_curve(path, row):
    if path.suffix.lower() == '.npz':
        with np.load(path, allow_pickle=False) as arr:
            for field, key in [('TICID', 'tic'), ('SECTOR', 'sector')]:
                if int(arr[field].item()) != int(row[key]):
                    raise ValueError(f'{path}: {field} does not match manifest')
            meta = {k.lower(): arr[k].item() if k in arr else '' for k in
                    ('CAMERA', 'CCD', 'TESSMAG', 'RA_OBJ', 'DEC_OBJ')}
            aux = {k: arr[k].astype(float) for k in AUX_COLUMNS if k in arr}
            return Curve(arr['TIME'].astype(float), arr['SAP_FLUX'].astype(float),
                         arr['QUALITY'] == 0, aux, meta)
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdus:
        primary = hdus[0].header
        table = hdus[1].data
        hdr = hdus[1].header
        for col, key in [('TICID', 'tic'), ('SECTOR', 'sector')]:
            value = primary.get(col, hdr.get(col))
            if value is None or int(value) != int(row[key]):
                raise ValueError(f'{path}: {col}={value}, manifest={row[key]}')
        names = {n.upper(): n for n in table.names}
        def data(key):
            return np.asarray(table[names[key]], dtype=float)
        t, f = data('TIME'), data('SAP_FLUX')
        good = data('QUALITY') == 0
        aux = {}
        for key in AUX_COLUMNS:
            if key in names:
                aux[key] = data(key)
        meta = {key.lower(): primary.get(key, hdr.get(key, '')) for key in
                ('CAMERA', 'CCD', 'TESSMAG', 'RA_OBJ', 'DEC_OBJ')}
    return Curve(t, f, good, aux, meta)


class Binning:
    def __init__(self, curve):
        finite_times = curve.time[np.isfinite(curve.time)]
        if len(finite_times) < 2 or np.any(np.diff(finite_times) <= 0):
            raise ValueError('Timestamps must be finite at least twice and strictly ordered')
        self.lo, self.hi = finite_times[[0, -1]]
        self.span = self.hi - self.lo
        use = curve.good & np.isfinite(curve.time) & np.isfinite(curve.flux)
        self.time, self.flux = curve.time[use], curve.flux[use]
        self.bin_ids = np.floor((self.time - self.lo) / BIN).astype(int)
        n = int(np.floor(self.span / BIN)) + 1
        self.centers = self.lo + (np.arange(n) + .5) * BIN
        self.count = np.bincount(self.bin_ids, minlength=n)
        cadence = np.median(np.diff(finite_times))
        self.valid = self.count >= max(1, math.ceil(.5 * BIN / cadence))
        self.norm = np.median(self.flux) if len(self.flux) else np.nan
        if not np.isfinite(self.norm) or self.norm <= 0:
            self.valid[:] = False

    def bin_flux(self, flux):
        y = np.full(len(self.count), np.nan)
        if np.any(self.valid):
            sums = np.bincount(self.bin_ids, weights=flux, minlength=len(y))
            y[self.valid] = sums[self.valid] / (self.count[self.valid] * self.norm)
        return y


def detect(y, centers):
    events = []
    for hours, fraction, idx, kernel in templates():
        if len(y) < len(kernel):
            continue
        depth = np.correlate(y, kernel, mode='valid')
        ok = np.isfinite(depth)
        if ok.sum() < 30:
            continue
        # Empirical noise on this template's timescale in this star-sector.
        noise = max(mad(depth[ok]), 1e-8)
        score = np.where(ok, depth / noise, -np.inf)
        peaks = find_peaks(np.r_[-np.inf, score, -np.inf], height=MIN_SCORE)[0] - 1
        for i in peaks:
            if depth[i] < MIN_DEPTH:
                continue
            peak = centers[i - idx[0]]
            duration = hours / 24
            events.append(dict(peak_btjd=float(peak), duration_h=hours,
                               ingress_fraction=fraction, depth=float(depth[i]),
                               score=float(score[i]),
                               start=float(peak - duration * fraction),
                               stop=float(peak + duration * (1 - fraction))))
    kept = []
    for event in sorted(events, key=lambda e: -e['score']):
        if all(event['stop'] <= e['start'] or event['start'] >= e['stop'] for e in kept):
            kept.append(event)
    return sorted(kept, key=lambda e: e['peak_btjd'])


def diagnostics(curve, event):
    start, stop = event['start'], event['stop']
    duration = stop - start
    t = curve.time
    inside = curve.good & (t >= start) & (t <= stop)
    outside = curve.good & (((t >= start - duration/2) & (t < start)) |
                            ((t > stop) & (t <= stop + duration/2)))
    values = {}
    for key, raw in curve.diagnostics.items():
        a, b = raw[inside], raw[outside]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) >= 3 and len(b) >= 3:
            shift = np.median(a) - np.median(b)
            values[key] = {'inside_minus_shoulder': float(shift),
                           'shoulder_mad': float(mad(b))}
    return json.dumps(values, sort_keys=True)


def matched(events, peak, duration):
    return any(abs(e['peak_btjd'] - peak) <= max(BIN, .25 * duration)
               for e in events)


def scan_one(curve, tic, sector, seed, trials):
    bins = Binning(curve)
    baseline = detect(bins.bin_flux(bins.flux), bins.centers)
    for ev in baseline:
        ev['diagnostics'] = diagnostics(curve, ev)
        ev.update(tic=tic, sector=sector, **curve.meta)
    # Stable if the manifest is reordered or shards of it are run separately.
    rng = np.random.default_rng(np.random.SeedSequence([seed, int(tic), int(sector)]))
    injections = []
    for depth in DEPTHS:
        for hours in DURATIONS:
            duration = hours / 24
            for _ in range(trials):
                peak = rng.uniform(bins.lo, bins.hi)
                fraction = float(rng.choice([.2, .8]))
                nearby = matched(baseline, peak, duration)
                model = triangle(bins.time, peak, duration, fraction)
                recovered = False
                if not nearby and np.any(model > 0) and bins.valid.any():
                    y = bins.bin_flux(bins.flux * (1 - depth * model))
                    recovered = matched(detect(y, bins.centers), peak, duration)
                injections.append(dict(tic=tic, sector=sector, depth=depth, duration_h=hours,
                                       peak_btjd=peak, ingress_fraction=fraction,
                                       baseline_coincidence=int(nearby), recovered=int(recovered)))
    audit = dict(span_d=bins.span, start_btjd=bins.lo, stop_btjd=bins.hi,
                 rows=len(curve.time), good_flux_rows=len(bins.time),
                 usable_halfhour_bins=int(bins.valid.sum()), n_triggers=len(baseline))
    return audit, baseline, injections


def rates(audit, controls, n_events):
    # Independent uniform injection times in each curve. Weighted Hoeffding
    # lower bound allows heterogeneous probabilities and unequal spans.
    span = {(a['tic'], a['sector']): a['span_d'] for a in audit}
    cells = {}
    for i in controls:
        key = (i['depth'], i['duration_h'])
        counts = cells.setdefault(key, {}).setdefault((i['tic'], i['sector']), [0, 0])
        counts[0] += i['recovered']
        counts[1] += 1
    results = []
    for (depth, hours), counts in cells.items():
        exposure = sum(span[k] * hits / n for k, (hits, n) in counts.items())
        margin = math.sqrt(.5 * math.log(1/.025) * sum(span[k]**2/n for k, (_, n) in counts.items()))
        lower = max(0, exposure - margin)
        numerator = .5 * chi2.ppf(.975, 2 * (1 + n_events))
        results.append(dict(depth=depth, duration_h=hours, trials=sum(x[1] for x in counts.values()),
                            recovered=sum(x[0] for x in counts.values()),
                            effective_star_days=exposure, exposure_lower_97p5=lower,
                            all_trigger_count=n_events, rate_upper_95_per_star_day=(
                                float(numerator / lower) if lower > 0 else None)))
    return results


def write_csv(path, records, fields):
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--injections', type=int, default=8, help='independent trials per depth-duration cell per star-sector')
    parser.add_argument('--seed', type=int, default=8462852)
    args = parser.parse_args()
    if args.injections < 1:
        parser.error('--injections must be positive')
    with args.manifest.open(newline='') as fh:
        rows = list(csv.DictReader(fh))
    keys = [(int(r['tic']), int(r['sector'])) for r in rows]
    if not rows or len(keys) != len(set(keys)):
        parser.error('Manifest must contain unique star-sectors and at least one row')
    args.out.mkdir(parents=True, exist_ok=True)
    audit, events, controls = [], [], []
    for row, (tic, sector) in zip(rows, keys):
        path = args.manifest.parent / row['path']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if row.get('sha256') and row['sha256'].lower() != digest:
            raise ValueError(f'Checksum mismatch: {path}')
        curve = read_curve(path, row)
        a, e, c = scan_one(curve, tic, sector, args.seed, args.injections)
        audit.append(dict(tic=tic, sector=sector, sha256=digest, **curve.meta, **a))
        events.extend(e)
        controls.extend(c)
    for event in events:
        event['same_detector_other_stars'] = sum(
            e['tic'] != event['tic'] and e['sector'] == event['sector']
            and event['camera'] != '' and event['ccd'] != ''
            and e['camera'] == event['camera'] and e['ccd'] == event['ccd']
            and abs(e['peak_btjd'] - event['peak_btjd']) <= 2 * BIN for e in events)
    write_csv(args.out/'audit.csv', audit, list(audit[0]))
    write_csv(args.out/'triggers.csv', events, ['tic', 'sector','camera', 'ccd', 'tessmag',
              'ra_obj', 'dec_obj', 'peak_btjd', 'start', 'stop', 'duration_h', 'ingress_fraction',
              'depth', 'score', 'same_detector_other_stars', 'diagnostics'])
    write_csv(args.out/'injections.csv', controls, list(controls[0]))
    report = dict(n_star_sectors=len(rows), n_distinct_tics=len(set(k[0] for k in keys)),
                  all_trigger_count=len(events), seed=args.seed, trials_per_cell=args.injections,
                  manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                  limits=rates(audit, controls, len(events)))
    (args.out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
