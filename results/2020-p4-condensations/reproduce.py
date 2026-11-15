"""Deblend the western C/2020 P4 condensation in native 2020 August 7 LASCO C2 pixels."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates, maximum_filter
from scipy.optimize import least_squares
from scipy.special import erf

HERE = Path(__file__).resolve().parent


def materialize(row, cache):
    path = cache / row['name']
    if not path.exists():
        with urlopen(row['url'], timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ValueError(f"Archive bytes changed for {path.name}")
        path.write_bytes(data)
    if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
        raise ValueError(f"Digest mismatch for {path.name}")
    return path


def hours(stamp):
    h, m, s = map(float, stamp.split(':'))
    return h + m / 60 + s / 3600


def image(path):
    # These fixed, uncompressed Level 0.5 products use a signed big-endian primary HDU.
    data = path.read_bytes()
    header = {}
    for offset in range(0, len(data), 80):
        card = data[offset:offset + 80].decode('ascii')
        if card[:8].strip() == 'END':
            break
        if card[8:10] == '= ':
            header[card[:8].strip()] = card[10:].split('/')[0].strip().strip("'").strip()
    start = (offset + 80 + 2879) // 2880 * 2880
    if [int(header[k]) for k in ('BITPIX', 'NAXIS1', 'NAXIS2')] != [16, 1024, 1024]:
        raise ValueError(path)
    if abs(float(header['CROTA2'])) < 170 or header['FILTER'] != 'Orange' or header['POLAR'] != 'Clear':
        raise ValueError(path)
    raw = np.frombuffer(data, '>i2', 1024**2, offset=start).reshape(1024, 1024)
    # With this roll the public upper-left report coordinates are (1023 - column, row).
    cut = raw[:220, :484][:, ::-1].astype(float) / float(header['EXPTIME'])
    return cut, hours(header['TIME-OBS']), float(header['CDELT1'])


def guide(path, times):
    rows = [line.split() for line in path.read_text().splitlines() if ' C2 ' in line]
    t = np.array([hours(r[1]) for r in rows]) - 15
    xy = np.array([(float(r[3]), 1023 - float(r[4])) for r in rows])
    return np.stack([np.polyval(coef, times - 15) for coef in np.polyfit(t, xy, 2).T], axis=1)


def gaussian(x, y, x0, y0, sx, sy, angle):
    c, s = np.cos(angle), np.sin(angle)
    u, v = c * (x - x0) + s * (y - y0), -s * (x - x0) + c * (y - y0)
    return np.exp(-0.5 * ((u / sx)**2 + (v / sy)**2))


def evaluate(p, x, y, double):
    bg, bx, by, amplitude, mx, my, sx, sy, angle, *extra = p
    core = gaussian(x, y, mx, my, sx, sy, angle)
    if double:
        height, dx, dy, width, length, rotation = extra
        flux = amplitude * core + height * gaussian(x, y, mx + dx, my + dy, width, length, rotation)
    else:
        flux = amplitude * core * (1 + erf((extra[0] * (x - mx) + extra[1] * (y - my)) / 2))
    return bg + bx * x + by * y + flux


def fit(x, y, z, double):
    amplitude = max(2, 0.8 * np.quantile(z, .97))
    initial = [0, 0, 0, amplitude, -2, -1, 1.5, 1, 0]
    lower = [-7, -1, -1, .1, -5, -3, .5, .5, -np.pi]
    upper = [7, 1, 1, 100, 1.5, 2.5, 3.5, 3, np.pi]
    if double:
        initial += [amplitude * .65, 2.5, 1, .9, .9, 0]
        lower += [.01, .7, -1, .55, .55, -np.pi]
        upper += [80, 5, 4, 2.5, 2.5, np.pi]
    else:
        initial += [.4, .3]
        lower += [-8, -8]
        upper += [8, 8]
    return least_squares(lambda p: evaluate(p, x, y, double) - z, initial,
                         bounds=(lower, upper), loss='soft_l1', f_scale=1.5, max_nfev=500).x


def held_error(p, x, y, z, double):
    return np.mean(np.minimum((evaluate(p, x, y, double) - z)**2, 16))


def measure(frame, center):
    cx, cy = center - [540, 0]
    xx = np.arange(round(cx) - 8, round(cx) + 5)
    yy = np.arange(round(cy) - 6, round(cy) + 7)
    x, y = np.meshgrid(xx - cx, yy - cy)
    mask = (x <= 3.4) & (abs(y) <= 5.5)
    parity = (np.add.outer(yy, xx)[mask] % 2).astype(bool)
    x, y, z = x[mask], y[mask], frame[np.ix_(yy, xx)][mask]
    double = fit(x, y, z, True)
    scores = []
    for multi in (False, True):
        total = 0
        for train in (parity, ~parity):
            p = fit(x[train], y[train], z[train], multi)
            total += (~train).sum() * held_error(p, x[~train], y[~train], z[~train], multi)
        scores.append(np.sqrt(total / len(z)))
    west = center + double[4:6]
    east = west + double[10:12]
    return [*west, *east, np.linalg.norm(east - west), *scores]


def stack(frames, centers, indices, x, y, order):
    return np.median([map_coordinates(frames[i], [y + centers[i, 1], x + centers[i, 0] - 540],
                                      order=order) for i in indices], axis=0)


def summarize_image(im, x, y):
    # Search for two local maxima in the western complex only. A plateau counts once per grid peak here.
    zone = (x >= -5) & (x <= 2) & (abs(y) <= 3)
    peak = (im == maximum_filter(im, size=5)) & zone & (im > 3)
    return [[float(x[i, j]), float(y[i, j]), round(float(im[i, j]), 2)] for i, j in zip(*np.where(peak))]


def field_control(frames, times, indices, anchor):
    # Bright features in a preliminary stack aligned with the common stellar drift.
    x, y = np.meshgrid(np.arange(-5, 6.), np.arange(-5, 6.))
    fits = []
    for i in indices:
        cx, cy = np.array(anchor) + np.array([12, .9]) * (times[i] - 15)
        px, py = round(cx), round(cy)
        patch = frames[i, py-5:py+6, px-5:px+6]
        initial = [0, max(1, patch[4:7, 4:7].max()), cx-px, cy-py, 1.2, 1.2, 0]
        params = least_squares(
            lambda q: (q[0] + q[1] * gaussian(x, y, *q[2:]) - patch).ravel(), initial,
            bounds=([-4, 0, -2, -2, .5, .5, -3.2], [4, 100, 2, 2, 3, 3, 3.2]),
            loss='soft_l1').x
        fits.append([times[i] - 15, px + params[2], py + params[3]])
    fitted = np.asarray(fits)
    motion = np.polyfit(fitted[:, 0], fitted[:, 1:], 1)
    residual = np.linalg.norm(fitted[:, 1:] - (fitted[:, :1] * motion[0] + motion[1]), axis=1)
    use = residual < np.quantile(residual, .75)
    motion = np.polyfit(fitted[use, 0], fitted[use, 1:], 1)
    ys, xs = np.mgrid[-6:6.01:.25, -6:6.01:.25]
    centers = (times[:, None] - 15) * motion[0] + motion[1] + [540, 0]
    aligned = stack(frames, centers, indices, xs, ys, 3)
    peaks = (aligned == maximum_filter(aligned, 7)) & (aligned > .45 * aligned.max())
    return {'xy_at_15_utc': (motion[1] + [540, 0]).tolist(),
            'drift_px_per_hour': motion[0].tolist(), 'peaks_above_45pct': int(peaks.sum())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    with (HERE / 'inputs.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    with ThreadPoolExecutor(5) as pool:
        paths = list(pool.map(lambda r: materialize(r, args.cache), rows))
    fits = [p for p in paths if p.suffix == '.fts']
    images, times, scales = zip(*(image(p) for p in fits))
    times, frames = np.array(times), np.stack(images)
    frames -= np.median(frames, axis=0)
    frames -= gaussian_filter(frames, (0, 9, 9))
    centers = guide(args.cache / 'soho4047_xy.txt', times)
    reference_a = guide(args.cache / 'soho4046_xy.txt', times)
    reference_c = guide(args.cache / 'soho4049_xy.txt', times)
    indices = np.flatnonzero((times > 12.55) & (times < 17.7))
    native = np.array([measure(frames[i], centers[i]) for i in indices])
    with (args.output / 'astrometry.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['fits', 'ut_hours', 'west_x', 'west_y', 'east_x', 'east_y', 'sep_px',
                         'single_held_rms', 'double_held_rms'])
        for i, row in zip(indices, native):
            writer.writerow([fits[i].name, f'{times[i]:.6f}', *[f'{n:.4f}' for n in row]])
    y, x = np.mgrid[-9:10:.25, -9:13:.25]
    splits = [('All exposures', indices), ('Odd exposures', indices[::2]), ('Even exposures', indices[1::2])]
    splits += [(label, group) for label, group in zip(['First third', 'Middle third', 'Last third'], np.array_split(indices, 3))]
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.3), constrained_layout=True, sharex=True, sharey=True)
    peaks = {}
    for ax, (label, group) in zip(axes.flat, splits):
        im = stack(frames, centers, group, x, y, 3)
        ax.pcolormesh(x, y, im, shading='auto', cmap='magma', vmin=0, vmax=18)
        ax.contour(x, y, im, levels=[4, 8, 10, 12, 14, 18, 24], colors='white', linewidths=.45)
        ax.set_title(f'{label} ({len(group)})', fontsize=11)
        ax.set_aspect('equal')
        peaks[label] = {'cubic': summarize_image(im, x, y),
                        'linear': summarize_image(stack(frames, centers, group, x, y, 1), x, y)}
    fig.suptitle('C/2020 P4  |  LASCO C2  |  7 August 2020', fontsize=16)
    fig.supxlabel('Column offset from component B guide (C2 pixels)')
    fig.supylabel('Row offset, positive toward image bottom (C2 pixels)')
    fig.savefig(args.output / 'condensations.png', dpi=145)
    plt.close(fig)
    dt = np.median(np.diff(times[indices]))
    summary = {'fit_exposures': len(indices), 'background_exposures_including_fits': len(fits),
               'plate_scale_arcsec': float(np.median(scales)), 'median_separation_px': np.median(native[:, 4]),
               'separation_quartiles_px': np.quantile(native[:, 4], [.25, .75]).tolist(),
               'median_separation_arcsec': np.median(native[:, 4]) * np.median(scales),
               'double_better_held_pixel_frames': int(np.sum(native[:, 6] < native[:, 5])),
               'median_held_rms_single_double': np.median(native[:, 5:], axis=0).tolist(),
               'median_guide_separation_B_A_C_px': [float(np.median(np.linalg.norm((p - centers)[indices], axis=1))) for p in (reference_a, reference_c)],
               'index_guide_moves_px': float(np.linalg.norm(centers[indices[-1]] - centers[indices[0]])),
               'field_controls': [field_control(frames, times, indices, anchor) for anchor in
                                  ([170, 189], [353, 39], [198, 133], [82, 109])],
               'typical_exposure_spacing_hr': float(dt), 'stack_maxima_x_y_dn_per_s': peaks}
    (args.output / 'morphology.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
