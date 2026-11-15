"""Recalculate closest approach and the optional ZTF image diagnostics."""
import argparse
import csv
import json
import math
from datetime import datetime, timedelta
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
AU_KM = 149597870.7
DAY = 86400.0


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def subtract(a, b):
    return [x - y for x, y in zip(a, b)]


def local_minimum(first, second):
    r = subtract(first['xyz_au'], second['xyz_au'])
    v = subtract(first['vxyz_au_per_day'], second['vxyz_au_per_day'])
    dt = -dot(r, v) / dot(v, v)
    closest = [x + dt * y for x, y in zip(r, v)]
    return dt, math.sqrt(dot(closest, closest)) * AU_KM, math.sqrt(dot(v, v)) * AU_KM / DAY


def geometry():
    evidence = json.loads((ROOT / 'ephemerides.json').read_text())
    for service, epoch in [('jpl', 'jd_tdb'), ('miriade', 'jd_tt')]:
        left, right = (evidence[service][key] for key in evidence['objects'])
        if service == 'jpl':
            epochs = [row[epoch] for row in left]
            if epochs != [row[epoch] for row in right]:
                raise ValueError('Unmatched vector epochs.')
            distances = [math.dist(a['xyz_au'], b['xyz_au']) * AU_KM for a, b in zip(left, right)]
            i = distances.index(min(distances))
            if not 0 < i < len(distances) - 1:
                raise ValueError('Closest approach is outside the sampled bracket.')
            print('JPL separation samples in km', [round(x) for x in distances])
            left, right = left[i], right[i]
        dt, distance, speed = local_minimum(left, right)
        jd = left[epoch] + dt
        when = datetime(2000, 1, 1, 12) + timedelta(days=jd - 2451545)
        print(service, when.isoformat(timespec='seconds'), epoch[3:].upper(),
              f'{distance:.1f} km, {speed:.5f} km/s')


def measure_frames(folder):
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    from scipy.ndimage import map_coordinates, median_filter
    from scipy.signal import fftconvolve

    rows = list(csv.DictReader((ROOT / 'ztf_frames.csv').open()))
    folder.mkdir(exist_ok=True, parents=True)
    print('exposure,band,comet_snr,asteroid_snr,asteroid_peak_offset_px,peak_on_search_boundary')
    for r in rows:
        path = folder / (r['exposure'] + '.fits')
        if not path.exists():
            with urlopen(r['cutout_url'], timeout=120) as response:
                blob = response.read()
            if not blob.startswith(b'SIMPLE'):
                raise ValueError('IRSA did not return a FITS image.')
            path.write_bytes(blob)
        with fits.open(path) as hdus:
            pixels = hdus[0].data.astype(float)
            header = hdus[0].header
        wcs = WCS(header)
        pixels -= median_filter(pixels, size=15)
        y, x = np.mgrid[-8:9, -8:9]
        sigma = header['SEEING'] / (2.355 * 1.012)
        kernel = np.exp(-(x*x + y*y) / (2*sigma*sigma))
        amplitudes = fftconvolve(pixels, kernel / np.sum(kernel**2), mode='same')
        med = np.nanmedian(amplitudes)
        noise = 1.4826 * np.nanmedian(np.abs(amplitudes - med))
        snrs = []
        offset = ['', '']
        for key in ['comet', 'asteroid']:
            pos = [float(r[key + '_' + axis + '_deg']) for axis in ['ra', 'dec']]
            px, py = wcs.all_world2pix(*pos, 0)
            if 8 <= px < pixels.shape[1]-8 and 8 <= py < pixels.shape[0]-8:
                value = map_coordinates(amplitudes, [[py], [px]], order=1)[0]
                snrs.append(f'{(value-med)/noise:.2f}')
                if key == 'asteroid':
                    dy, dx = np.mgrid[-2:2.01:0.2, -2:2.01:0.2]
                    window = map_coordinates(amplitudes, [py + dy, px + dx], order=3)
                    peak = np.unravel_index(np.argmax(window), window.shape)
                    shift = math.hypot(dx[peak], dy[peak])
                    edge = max(abs(dx[peak]), abs(dy[peak])) > 1.95
                    offset = [f'{shift:.2f}', str(edge).lower()]
            else:
                snrs.append('')
        print(','.join([r['exposure'], r['band'], *snrs, *offset]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frames', type=Path, help='Use/download FITS cutouts here and measure')
    args = parser.parse_args()
    if args.frames:
        measure_frames(args.frames)
    else:
        geometry()
