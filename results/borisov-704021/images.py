"""Reproduce forced aperture measurements, source positions and inter-frame motion."""
from pathlib import Path
import csv
import json
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve

DATA = Path(__file__).resolve().parent / 'data'


def read_image(row):
    with fits.open(DATA / (row['image'] + '.fits')) as hdus:
        frame = hdus[0].data.astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)
            wcs = WCS(hdus[0].header)
    coords = np.array([float(row['ra_deg']), float(row['dec_deg'])])
    x, y = wcs.all_world2pix(coords[None, :], 0)[0]
    # Use the same compensated circular aperture and control offsets in every image.
    yy, xx = np.mgrid[-14:15, -14:15]
    rad = np.hypot(xx, yy)
    inner = rad <= 6
    outer = (rad >= 9) & (rad <= 14)
    kernel = inner.astype(float) - outer * (inner.sum() / outer.sum())
    flux = fftconvolve(frame - np.median(frame), kernel, mode='same')
    offsets = np.arange(-40, 41, 4)
    cx, cy = np.meshgrid(offsets, offsets)
    use = np.hypot(cx, cy) >= 16
    controls = map_coordinates(flux, [y + cy[use], x + cx[use]], order=1)
    forced = map_coordinates(flux, [[y], [x]], order=1)[0]
    med = np.median(controls)
    noise = 1.4826 * np.median(abs(controls - med))
    # Centroid is a second measurement, independent from the fixed aperture score.
    frame_y, frame_x = np.indices(frame.shape)
    separation = np.hypot(frame_x - x, frame_y - y)
    source = separation <= 9
    sky = np.median(frame[(separation >= 11) & (separation <= 17)])
    weights = (frame - sky)[source]
    measured = np.array([np.sum(weights * axis[source]) / weights.sum() for axis in (frame_x, frame_y)])
    skypos = wcs.all_pix2world(measured[None, :], 0)[0]
    difference = (skypos - coords) * [np.cos(np.deg2rad(coords[1])), 1] * 3600
    return {'image': row['image'], 'mjd_utc': float(row['mjd_utc']),
            'forced_flux_adu': round(forced, 2), 'local_sigma_adu': round(noise, 2),
            'forced_score': round((forced - med) / noise, 2),
            'centroid_ra_deg': round(skypos[0], 8), 'centroid_dec_deg': round(skypos[1], 8),
            'residual_ra_cosdec_arcsec': round(difference[0], 3), 'residual_dec_arcsec': round(difference[1], 3)}, (controls - med) / noise


def main():
    with open(DATA / 'imaging-ephemeris.csv') as stream:
        rows = list(csv.DictReader(stream))
    measurements = []
    for r in rows:
        result, controls = read_image(r)
        measurements.append((result, controls))
    nights = []
    for day in sorted(set(int(i[0]['mjd_utc']) for i in measurements)):
        sub = [(i, c) for i, c in measurements if int(i['mjd_utc']) == day]
        scores = np.array([i['forced_score'] for i, c in sub])
        # A shared offset shifts the predicted moving track in every constituent exposure.
        controls = np.mean([c for i, c in sub], axis=0)
        nights.append({'mjd': day, 'number_images': len(sub), 'offset_tracks': len(controls),
                       'controls_at_least_target': int(np.sum(controls >= scores.mean()))})
    print(json.dumps({'images': [i for i, c in measurements], 'control_tracks': nights}, indent=2))


if __name__ == '__main__':
    main()
