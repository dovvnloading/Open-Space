"""Reproduce the HST image comparisons and the full-field jet-sector counts."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter, shift
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
# Position tested on the 2024 chart and a blind prediction from the 2016 table
# using the typical frame correction. Difference-image centroids refine the latter.
POSITIONS = {
    23: ('12h30m47.88s +12d23m54.01s', '12h30m47.88s +12d22m53.95s'),
    32: ('12h30m46.92s +12d22m24.20s', '12h30m46.92s +12d22m34.26s'),
    38: ('12h30m45.47s +12d23m48.08s', '12h30m45.46s +12d24m48.12s'),
}
SCALE = .05  # arcsec per pixel


def load_images():
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    images = {}
    y, x = np.mgrid[-24:25, -24:25]
    with np.load(ROOT / 'cutouts.npz', allow_pickle=False) as cutouts:
        for key, info in manifest.items():
            raw = cutouts[key]
            digest = hashlib.sha256(raw.astype('>f4').tobytes()).hexdigest()
            if digest != info['sha256_pixels_be_float32']:
                raise ValueError(f'Pixel digest mismatch for {key}')
            number, site, visit, band = key.split('-')
            coord = SkyCoord(POSITIONS[int(number)][int(site)])
            wcs = WCS(fits.Header.fromstring(info['cutout_header'], sep='\n'), naxis=2)
            dec = coord.dec.deg + y * SCALE / 3600
            ra = coord.ra.deg + x * SCALE / (3600 * np.cos(coord.dec.radian))
            xx, yy = wcs.all_world2pix(ra, dec, 0)
            if min(xx.min(), yy.min()) < 2 or max(xx.max(), yy.max()) > 74:
                raise ValueError(f'Cutout edge for {key}')
            images[int(number), int(site), band, visit] = map_coordinates(
                raw.astype(float), [yy, xx], order=3)
    return images, x, y


def register_images(frames, reference, comparison):
    """Translate images by optimizing robust background structure outside the source."""
    highpass = {key: frame - median_filter(frame, size=11)
                for key, frame in frames.items()}

    def normalized(a):
        scale = np.median(np.abs(a - np.median(a))) * 1.4826
        return np.clip(a / scale, -3, 8)

    ref = normalized(highpass[reference])
    registered = {}
    for key, frame in highpass.items():
        if key == reference:
            registered[key] = frame
            continue
        norm = normalized(frame)

        def loss(offset):
            test = shift(norm, offset, order=1)
            return -np.corrcoef(ref[comparison], test[comparison])[0, 1]

        fit = minimize(loss, (0, 0), method='Powell', bounds=[(-2, 2)] * 2)
        registered[key] = shift(frame, fit.x, order=3)
    return registered


def compare(images, x, y):
    tests = {
        23: [('f606w', ['12'], ['31']), ('f606w', ['11'], ['31']),
             ('f814w', ['12'], ['31']), ('f814w', ['11'], ['31'])],
        32: [('f606w', ['01'], ['a1']), ('f814w', ['01'], ['a1']),
             ('f814w', ['02'], ['a1'])],
        38: [('f814w', ['27', '31', '44', 'a1'], ['rjcl403']),
             ('f814w', ['27', '31', '44', 'a1'], ['01', '02']),
             ('f814w', ['31', '44'], ['rjcl403']),
             ('f814w', ['18', '27'], ['01', '02'])]
    }
    r = np.hypot(x, y)
    field = (r > 5) & (r < 22)
    measured, panels = [], {}
    for number, comparisons in tests.items():
        for site, name in enumerate(('2024', '2016-adjusted')):
            for band, bright, faint in comparisons:
                frames = {v: images[number, site, band, v] for v in bright + faint}
                if number == 38:
                    frames['31'] = images[number, site, band, '31']
                    frames = register_images(frames, '31', field)
                diff = (np.median([frames[v] for v in bright], axis=0)
                        - np.median([frames[v] for v in faint], axis=0))
                if number != 38:
                    diff -= median_filter(diff, size=9)
                sm = gaussian_filter(diff, .9)
                noise = 1.4826 * np.median(np.abs(sm[field] - np.median(sm[field])))
                score = (sm - np.median(sm[field])) / noise
                peak = np.unravel_index(np.argmax(score[21:28, 21:28]), (7, 7))
                iy, ix = np.array(peak) + 21
                region = ((x - x[iy, ix])**2 + (y - y[iy, ix])**2 < 5)
                weights = np.maximum(score, 0) * region
                east = SCALE * np.sum(x * weights) / weights.sum()
                north = SCALE * np.sum(y * weights) / weights.sum()
                measured.append(dict(nova=number, site=name, filter=band,
                                     bright='+'.join(bright), reference='+'.join(faint),
                                     east_arcsec=round(east, 3), north_arcsec=round(north, 3),
                                     peak_local_contrast=round(score[iy, ix], 2),
                                     control_field_mad=round(noise, 6)))
                panels.setdefault((number, name), score)
    return measured, panels


def wedge_counts(corrections):
    rows = list(csv.DictReader((ROOT / 'positions.csv').open()))
    centers = SkyCoord('12h30m49.423s +12d23m28.04s')
    ra, dec = np.array([(float(r['ra_deg']), float(r['dec_deg'])) for r in rows]).T
    result = {}
    for stage in ('published', 'remeasured'):
        if stage == 'remeasured':
            for i, row in enumerate(rows):
                n = int(row['id'])
                if row['survey'] == 'ACS' and n in corrections:
                    ra[i], dec[i] = corrections[n]
        east = (ra - centers.ra.deg) * np.cos(centers.dec.radian)
        north = dec - centers.dec.deg
        angle = np.degrees(np.arctan2(east, north)) % 360
        radius = 3600 * np.hypot(east, north)
        bins = np.floor(((289.713 - angle + 18) % 360) / 36).astype(int)
        result[stage] = {
            'full_field': np.bincount(bins, minlength=10).tolist(),
            'inner_31_5_arcsec': np.bincount(bins[radius < 31.5], minlength=10).tolist(),
            'novae': {int(row['id']): {'PA_deg': round(a, 2), 'radius_arcsec': round(r, 2)}
                      for row, a, r in zip(rows, angle, radius)
                      if row['survey'] == 'ACS' and int(row['id']) in corrections}
        }
    return result


def main():
    images, x, y = load_images()
    measurements, panels = compare(images, x, y)
    coords = {}
    for n, sites in POSITIONS.items():
        records = [r for r in measurements if r['nova'] == n and r['site'] == '2016-adjusted']
        east, north = np.median([(r['east_arcsec'], r['north_arcsec']) for r in records], axis=0)
        center = SkyCoord(sites[1])
        coords[n] = (center.ra.deg + east / (3600 * np.cos(center.dec.radian)),
                     center.dec.deg + north / 3600)
    with (ROOT / 'photometry.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, measurements[0].keys(), lineterminator='\n')
        w.writeheader()
        w.writerows(measurements)
    result = {'centroid_degrees': coords, 'sectors': wedge_counts(coords)}
    (ROOT / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    if '--plot' in sys.argv:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2, figsize=(6.5, 6.8), layout='constrained')
        for ax, (number, site) in zip(axes.flat, panels):
            ax.imshow(panels[number, site][11:38, 11:38], origin='lower', cmap='RdBu_r',
                      vmin=-8, vmax=8, extent=[-.675, .675, -.675, .675])
            ax.plot(0, 0, marker='o', ms=12, markerfacecolor='none', markeredgecolor='k', lw=.7)
            ax.set_title(f'Nova {number} · {site}', fontsize=10)
            ax.set_xticks([-.5, 0, .5])
            ax.set_yticks([-.5, 0, .5])
        fig.supxlabel('East offset (arcsec)')
        fig.supylabel('North offset (arcsec)')
        fig.savefig(ROOT / 'difference_images.png', dpi=150)
    print(result['centroid_degrees'], result['sectors']['remeasured']['full_field'])


if __name__ == '__main__':
    main()
