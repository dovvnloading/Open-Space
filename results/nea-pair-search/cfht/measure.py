"""Measure the 2017 RB15 sky track in three CFHT MegaCam cutouts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

HERE = Path(__file__).resolve().parent
PHASES = np.linspace(-0.5, 0.5, 41)


def trail(x, y, dx, dy, sigma):
    return np.mean(np.exp(-((x[..., None] - dx * PHASES)**2
                           + (y[..., None] - dy * PHASES)**2) / (2*sigma**2)), axis=-1)


def measure(frame, meta, path):
    with fits.open(path) as hdus:
        image_hdu = next(h for h in hdus if h.data is not None)
        image = image_hdu.data.astype(float)
        header = image_hdu.header.copy()
    wcs = WCS(header)
    sky = np.array([meta['direction_at_epoch_deg'], meta['direction_after_one_min_deg']])
    exposure = header['EXPTIME']
    midpoint_mjd = header['MJD-OBS'] + exposure / (2*86400)
    offset_minutes = (midpoint_mjd + 2400000.5 - meta['ssos_epoch_jd_utc']) * 1440
    prediction = sky[0] + (sky[1] - sky[0]) * offset_minutes
    pix = wcs.all_world2pix([prediction, prediction + (sky[1] - sky[0])], 0)
    center = np.rint(pix[0]).astype(int)
    dx, dy = (pix[1] - pix[0]) * exposure / 60
    x, y = np.meshgrid(np.arange(-30, 31), np.arange(-17, 18))
    patch = image[center[1]+y, center[0]+x]

    def residual(parameters):
        east_pix, north_pix, height, width, vx, vy, background = parameters
        return (background + height*trail(x-east_pix, y-north_pix, vx, vy, width) - patch).ravel()

    fit = least_squares(residual, [0,0,600,1.8,dx,dy,np.median(patch)],
                        bounds=([-12,-12,0,.6,-40,-20,0], [12,12,20000,8,0,20,10000]),
                        loss='soft_l1', f_scale=18)
    fitted = center + fit.x[:2]
    ra, dec = wcs.all_pix2world([fitted], 0)[0]
    cosdec = np.cos(np.radians(dec))
    off = ([ra, dec] - prediction) * [3600*cosdec, 3600]
    # Measure AT the frozen orbit position; no maximization or translation is used for this score.
    gx, gy = np.meshgrid(np.arange(-30,31), np.arange(-30,31))
    kernel = trail(gx, gy, dx, dy, 1.9)
    kernel /= kernel.sum()
    filtered = fftconvolve(image - gaussian_filter(image, 15), kernel, mode='same')
    yy, xx = np.mgrid[-210:211:16, -210:211:16]
    mask = np.hypot(xx, yy) >= 35
    ctrl = map_coordinates(filtered, [pix[0,1]+yy[mask], pix[0,0]+xx[mask]], order=1)
    at_prediction = map_coordinates(filtered, [[pix[0,1]], [pix[0,0]]], order=1)[0]
    med = np.median(ctrl)
    sd = 1.4826 * np.median(np.abs(ctrl-med))
    score = (at_prediction-med) / sd
    return dict(image=frame, site='568', utc_mid=Time(midpoint_mjd,format='mjd').utc.isot,
                ra_deg=round(ra,8), dec_deg=round(dec,8), predicted_ra_deg=round(prediction[0],8),
                predicted_dec_deg=round(prediction[1],8), east_arcsec=round(off[0],4),
                north_arcsec=round(off[1],4), separation_arcsec=round(np.linalg.norm(off),4),
                trail_length_px=round(float(np.linalg.norm(fit.x[4:6])),3),
                forced_contrast_mad=round(float(score),2), controls=len(ctrl)), (patch, pix[0]-center, fit.x[:2]), (ctrl-med)/sd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    rows, controls = [], []
    fig, axes = plt.subplots(1,3,figsize=(8.1,2.9),layout='constrained')
    for axis, (frame, meta) in zip(axes, json.loads((HERE/'inputs.json').read_text()).items()):
        path = args.cache / (frame + '.fits')
        if not path.exists():
            with urlopen(meta['cutout_url'], timeout=100) as response:
                path.write_bytes(response.read())
        with fits.open(path) as hdus:
            hdu = next(h for h in hdus if h.data is not None)
            props = {k: hdu.header[k] for k in ('MJD-OBS', 'EXPTIME', 'FILTER', 'OBJECT')}
            wcs = dict(WCS(hdu.header).to_header(relax=True))
            metadata = json.dumps([wcs, props, list(hdu.data.shape)], sort_keys=True,
                                  separators=(',', ':')).encode()
            pixels = np.asarray(hdu.data, dtype='<f8').tobytes(order='C')
            if hashlib.sha256(metadata + pixels).hexdigest() != meta['measurement_input_sha256']:
                raise ValueError('CADC image pixels, astrometry or exposure metadata changed: ' + frame)
        row, (patch,predicted,measured), ctrl = measure(frame, meta, path)
        rows.append(row)
        controls.append(ctrl)
        v = np.percentile(patch,[30,99.4])
        axis.imshow(patch,origin='lower',extent=(-30.5,30.5,-17.5,17.5),cmap='gray',vmin=v[0],vmax=v[1])
        axis.plot(*predicted,marker='+',color='#16a3b6',markersize=12,markeredgewidth=1.3)
        axis.set_title(f"{frame}   {row['utc_mid'][11:19]} UTC",fontsize=10)
        axis.set_xticks([]);axis.set_yticks([])
    fig.savefig(args.output/'cfht-trails.png',dpi=160)
    with (args.output/'measurements.csv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    score = np.mean([row['forced_contrast_mad'] for row in rows])
    averaged = np.mean(controls,axis=0)
    print(f'Combined contrast {score:.2f}; matched-offset controls >= target: '
          f'{np.sum(averaged >= score)} / {len(averaged)}')


if __name__ == '__main__':
    main()
