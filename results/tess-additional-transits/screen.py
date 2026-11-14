"""Replay the 90-target search from SHA-verified SPOC light curves."""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import groupby
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

HERE = Path(__file__).resolve().parent
CATALOG = pd.read_csv(HERE / 'ephemerides.csv')
SELECTION = pd.read_csv(HERE / 'selection.csv')
HARMONICS = np.array([.25, 1/3, .5, 2/3, 1, 1.5, 2, 3, 4])


def light_curve(path, tic, column):
    with np.load(path) as data:
        valid = (data['QUALITY'] == 0) & np.isfinite(data['TIME']) & np.isfinite(data[column])
        time = data['TIME'][valid]
        flux = data[column][valid].astype(float)
    flux /= np.median(flux)
    bins = np.floor((time - time.min()) * 144 + 1e-5).astype(int)
    counts = np.bincount(bins)
    good = counts >= 3
    time = np.bincount(bins, weights=time)[good] / counts[good]
    flux = np.bincount(bins, weights=flux)[good] / counts[good]
    for segment in np.split(np.arange(len(time)), np.flatnonzero(np.diff(time) > .12) + 1):
        if len(segment) <= 60:
            flux[segment] = np.nan
            continue
        width = min(145, len(segment) // 2 * 2 - 1)
        flux[segment] /= median_filter(flux[segment], size=width, mode='nearest')
    flux -= 1
    scatter = 1.4826 * np.nanmedian(abs(flux - np.nanmedian(flux)))
    good = np.isfinite(flux) & (flux < 5 * scatter) & (flux > -15 * scatter)
    for row in CATALOG[CATALOG.tid == tic].itertuples():
        period = row.pl_orbper
        epoch = row.pl_tranmid - 2457000
        if period > 0 and np.isfinite(epoch):
            duration = row.pl_trandurh / 24
            half_width = min(.14 * period, max(.035, duration if np.isfinite(duration) else .09))
            good &= abs((time - epoch + .5 * period) % period - .5 * period) > half_width
    return time[good], flux[good]


def depth(time, flux, period, epoch, duration):
    phase = abs((time - epoch + period / 2) % period - period / 2)
    inside = phase < duration / 2
    outside = (phase >= .8 * duration) & (phase < 2.5 * duration)
    if inside.sum() < 3 or outside.sum() < 5:
        return np.nan, np.nan, 0
    measured = flux[outside].mean() - flux[inside].mean()
    scatter = 1.4826 * np.median(abs(flux - np.median(flux)))
    error = scatter * np.sqrt(1 / inside.sum() + 1 / outside.sum())
    events = np.unique(np.rint((time[inside] - epoch) / period)).size
    return measured, measured / error, events


def screen_host(record, cache):
    tic, selected = record
    sectors = sorted(int(s) for s in selected.split())
    curves = {s: light_curve(cache / f'{tic}-{s}.npz', tic, 'PDCSAP_FLUX') for s in sectors}
    consecutive = [list(n for _, n in group) for _, group in
                   groupby(enumerate(sectors), key=lambda item: item[1] - item[0])]
    run = max(consecutive, key=lambda group: (len(group), group[-1]))
    if len(run) >= 3 and len(run) == len(sectors):
        run = run[:-1]
    others = [s for s in sectors if s not in run]
    time = np.concatenate([curves[s][0] for s in run])
    flux = np.concatenate([curves[s][1] for s in run])
    start = time[0]
    time -= start
    baseline = np.ptp(time)
    scatter = 1.4826 * np.median(abs(flux - np.median(flux)))
    model = BoxLeastSquares(time, flux, dy=scatter)
    periods = 1 / np.arange(1/13, 1/.55, .012/baseline)[::-1]
    coarse = model.power(periods, [.042, .063, .083, .111, .139], oversample=4)
    power = coarse.power.copy()
    catalog_periods = CATALOG[CATALOG.tid == tic].pl_orbper.dropna().to_numpy()
    rows, diagnostics = [], []
    for rank in range(1, 8):
        initial_period = periods[np.argmax(power)]
        grid = initial_period * np.linspace(1 - .03*initial_period/baseline,
                                            1 + .03*initial_period/baseline, 31)
        fine = model.power(grid, [.0417, .0556, .0694, .0833, .1042, .125, .1458], oversample=10)
        best = np.argmax(fine.power)
        period = grid[best]
        duration = fine.duration[best]
        epoch = fine.transit_time[best] + start
        matches = any(np.min(abs(period / p - HARMONICS)) / max(1, period / p) < .008
                      for p in catalog_periods if p > 0)
        signal = {s: depth(*curves[s], period, epoch, duration) for s in sectors}
        passes = (not matches and fine.depth_snr[best] >= 8 and
                  sum(signal[s][1] >= 3 for s in run) >= 2 and
                  any(signal[s][1] >= 4 for s in others))
        rows.append([tic, rank, period, epoch, duration * 24, fine.depth[best] * 1e6,
                     fine.depth_snr[best], ' '.join(map(str, run)),
                     ' '.join(map(str, others)), int(matches), int(passes)])
        for s in sectors:
            dep, snr, n = signal[s]
            diagnostics.append([tic, rank, s, 'PDCSAP_FLUX', n, dep * 1e6, snr])
            if passes:
                raw = light_curve(cache / f'{tic}-{s}.npz', tic, 'SAP_FLUX')
                raw_depth, raw_snr, raw_n = depth(*raw, period, epoch, duration)
                diagnostics.append([tic, rank, s, 'SAP_FLUX', raw_n, raw_depth * 1e6, raw_snr])
        power[abs(1/periods - 1/period) < .8/baseline] = -np.inf
    return rows, diagnostics


def save(path, labels, rows):
    with path.open('w', newline='') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(labels.split())
        for row in rows:
            writer.writerow([format(n, '.8f') if isinstance(n, (float, np.floating)) else n for n in row])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    jobs = sorted((int(row.tic), row.used_sectors) for row in SELECTION.dropna(subset=['used_sectors']).itertuples())
    with ProcessPoolExecutor(3) as pool:
        results = list(pool.map(partial(screen_host, cache=args.cache), jobs))
    save(args.output / 'peaks.csv', 'tic rank period_days epoch_btjd duration_hours depth_ppm bls_depth_snr training_sectors heldout_sectors harmonic_match passed_triage',
         [row for peaks, _ in results for row in peaks])
    save(args.output / 'sector_checks.csv', 'tic rank sector flux observed_epochs depth_ppm white_snr',
         [row for _, diagnostics in results for row in diagnostics])
