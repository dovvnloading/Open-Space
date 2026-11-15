"""Rebuild the TESS difference images and eclipse-by-eclipse controls.

Run with Python, numpy, scipy and astropy. Input FITS files from TESScut
are downloaded into the cache directory passed on the command line.
"""
import argparse
import hashlib
import io
import json
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent


def fit_image(image, target):
    """Fit an elliptical Gaussian and a sloped background near the target."""
    y_grid, x_grid = np.indices(image.shape)
    x0, y0 = target
    mask = (abs(x_grid - x0) < 2.8) & (abs(y_grid - y0) < 2.8)
    x, y, z = x_grid[mask], y_grid[mask], image[mask]
    amplitude = np.ptp(z)
    initial = [amplitude, x0, y0, 1, 1, 0, np.median(z), 0, 0]
    lower = [0, x0-2, y0-2, .4, .4, -.7, -np.inf, -np.inf, -np.inf]
    upper = [4*amplitude, x0+2, y0+2, 1.9, 1.9, .7, np.inf, np.inf, np.inf]

    def residual(params):
        amp, cx, cy, sx, sy, rho, offset, gx, gy = params
        u, v = (x-cx)/sx, (y-cy)/sy
        exponent = -.5 * (u*u + v*v - 2*rho*u*v) / (1-rho*rho)
        return amp*np.exp(exponent) + offset + (x-x0)*gx + (y-y0)*gy - z

    return least_squares(residual, initial, bounds=(lower, upper)).x[1:3]


def difference_events(times, flux, period, width):
    """Bracket each eclipse with two same-orbit local out-of-eclipse windows."""
    orbit = np.rint(times / period)
    phase_time = times - orbit*period
    differences = []
    cycles = []
    for cycle in np.unique(orbit):
        same = orbit == cycle
        inside = same & (abs(phase_time) < .37*width)
        left = same & (phase_time > -1.4*width) & (phase_time < -.7*width)
        right = same & (phase_time > .7*width) & (phase_time < 1.4*width)
        if min(inside.sum(), left.sum(), right.sum()) < 2:
            continue
        baseline = (np.median(flux[left], axis=0) + np.median(flux[right], axis=0))/2
        differences.append(baseline - np.median(flux[inside], axis=0))
        cycles.append(int(cycle))
    return np.array(differences), cycles


def aperture_flux(image, position):
    """Small source aperture, corrected by the median of the outer pixels."""
    y, x = np.indices(image.shape)
    r2 = (x-position[0])**2 + (y-position[1])**2
    return float((image[r2<2.0] - np.median(image[r2>12])).sum())


def measure(payload, target):
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        member, = archive.namelist()
        fits_bytes = archive.read(member)
    with fits.open(io.BytesIO(fits_bytes)) as hdus:
        data = hdus[1].data
        time = data['TIME'] + hdus[1].header['BJDREFI'] - target['t0_bjd']
        flux = data['FLUX'].astype(float)
        good = (data['QUALITY'] == 0) & np.isfinite(time) & np.isfinite(flux).all(axis=(1,2))
        time, flux = time[good], flux[good]
        wcs = WCS(hdus[2].header, naxis=2)
        nominal = np.array(wcs.world_to_pixel_values(*target['toi_position']))
        eb = np.array(wcs.world_to_pixel_values(*target['binary_position']))
    period = target['period_days']
    duration = target['duration_hours'] / 24
    events, cycles = difference_events(time, flux, period, duration)
    image = np.median(events, axis=0)
    centroid = fit_image(image, nominal)
    axis = (eb - nominal) / np.linalg.norm(eb - nominal)
    rng = np.random.default_rng(514)
    projected = []
    for _ in range(200):
        sample = events[rng.integers(len(events), size=len(events))]
        center = fit_image(np.median(sample, axis=0), nominal)
        projected.append(float((center-nominal) @ axis))
    signal = {str(shift): None for shift in [0, -.25, .25]}
    for shift in [0, -.25, .25]:
        ds, _ = difference_events(time-shift*period, flux, period, duration)
        if len(ds):
            signal[str(shift)] = aperture_flux(np.median(ds, axis=0), eb)
    return {
        'product': member, 'fits_sha256': hashlib.sha256(fits_bytes).hexdigest(),
        'cycles': cycles, 'good_cadences': int(good.sum()),
        'toi_xy': nominal.tolist(), 'binary_xy': eb.tolist(),
        'deficit_xy': centroid.tolist(), 'toi_to_binary_px': np.linalg.norm(eb-nominal),
        'distance_to_toi_px': np.linalg.norm(centroid-nominal),
        'distance_to_binary_px': np.linalg.norm(centroid-eb),
        'projected_offset_px': (centroid-nominal) @ axis,
        'projected_bootstrap_5_95_px': np.percentile(projected, [5,95]).tolist(),
        'phase_aperture_electrons_per_second': signal,
        'difference_image_electrons_per_second': image.round(5).tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('--output', type=Path, default=HERE/'measurements.json')
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    target = json.loads((HERE/'target.json').read_text())
    results = {}
    for sector in target['all_available_sectors']:
        params = dict(zip(('ra','dec'), target['toi_position']))
        params.update(x=9,y=9,units='px',sector=sector)
        url = 'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?' + urllib.parse.urlencode(params)
        archive = args.cache / f'3143-s{sector}.zip'
        if not archive.exists():
            with urllib.request.urlopen(url, timeout=180) as response:
                contents = response.read()
            with zipfile.ZipFile(io.BytesIO(contents)) as check:
                if check.testzip():
                    raise ValueError(f'Archive failed CRC checks for sector {sector}')
            archive.write_bytes(contents)
        entry = measure(archive.read_bytes(), target)
        entry['url'] = url
        results[sector] = entry
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
