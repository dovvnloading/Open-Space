"""Fetch and measure CFHT MegaPipe g/r cutouts of 5XMM J085618.5+380711."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import center_of_mass, maximum_filter
from scipy.signal import fftconvolve

RA, DEC = 134.077254, 38.119831
SHIFT_EAST, SHIFT_NORTH, R99 = 0.3474233, 0.9488474, 3.204
# Legacy Survey DR10 primary PSF sources, r < 21 within 30 arcsec (AB magnitudes).
STARS = [
    ('9907734376220281', 134.07515680058188, 38.125915802969764, 19.81570, 19.24263),
    ('9907734301708473', 134.07284376420003, 38.11411876200624, 21.58512, 20.55341),
]


def get(directory, band, weight):
    stem = f'MegaPipe.211.256.{band}'
    if weight:
        # The companion compressed plane lacks a SODA sky descriptor. Use the
        # exact FITS section returned by SODA for these co-registered cutouts.
        base = f'https://cadc-west-01.canfar.net/raven/files/cadc:CFHTSG/{stem}.weight.fits.fz'
        url = base + '?' + urlencode({'SUB': '[1][1742:2286,7053:7597]'})
    else:
        base = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync'
        url = base + '?' + urlencode({'ID': f'cadc:CFHTSG/{stem}.fits', 'CIRCLE': f'{RA} {DEC} 0.014'})
    path = directory / ('cfht-' + band + ('.weight.fits' if weight else '.fits'))
    if not path.exists():
        with urlopen(url, timeout=120) as stream:
            raw = stream.read()
        if not raw.startswith(b'SIMPLE'):
            raise ValueError('Expected a FITS cutout')
        path.write_bytes(raw)
    data, header = fits.getdata(path, header=True)
    pixels = np.ascontiguousarray(data, dtype='>f4').tobytes()
    return data, header, {'url': url, 'file_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                          'float32_big_endian_pixel_sha256': hashlib.sha256(pixels).hexdigest()}


def measure(directory, band):
    data, header, image = get(directory, band, False)
    weight, wh, weights = get(directory, band, True)
    wcs = WCS(header)
    for key in ('CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD2_2'):
        assert np.isclose(header[key], wh[key], rtol=1.e-9)
    assert weight.shape == data.shape == (545, 545)
    xp, yp = wcs.world_to_pixel_values(RA, DEC)
    xr, yr = wcs.world_to_pixel_values(RA + SHIFT_EAST / (3600*np.cos(np.radians(DEC))),
                                       DEC + SHIFT_NORTH / 3600)
    y, x = np.indices(data.shape)
    scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600)
    d0 = np.hypot(x - xp, y - yp) * scale
    d1 = np.hypot(x - xr, y - yr) * scale
    center_masks = [d0 <= R99, d1 <= R99]
    area = center_masks[0] | center_masks[1]
    # Include the whole photometry aperture for all possible positions in either circle.
    padded = np.minimum(d0, d1) <= R99 + 4.5 * scale
    valid = np.isfinite(data) & np.isfinite(weight) & (weight > 0) & (data < header['SATURATE'])
    if not np.all(valid[padded]):
        raise ValueError('Masked or saturated science pixels in the localization region')
    # Circular four-pixel (0.74 arcsec) aperture and estimate of total flux.
    dy, dx = np.mgrid[-4:5, -4:5]
    kernel = (dx**2 + dy**2 <= 16).astype(float)
    summed = fftconvolve(np.nan_to_num(data), kernel, mode='same')
    calibrators = []
    zp = float(header['PHOTZP'])
    for sid, ra, dec, gmag, rmag in STARS:
        sx, sy = wcs.world_to_pixel_values(ra, dec)
        dist = np.hypot(x - sx, y - sy)
        sky = np.median(data[(dist > 17) & (dist <= 22)])
        cy, cx = center_of_mass(np.maximum(data - sky, 0) * (dist <= 6))
        rad = np.hypot(x - round(cx), y - round(cy))
        small, total = [(data[rad <= k] - sky).sum(dtype=float) for k in (4, 15)]
        reference_mag = gmag if band.startswith('G') else rmag
        calibrators.append(dict(ls_id=sid, correction_fraction=small / total,
                                large_aperture_ab=float(zp - 2.5*np.log10(total)),
                                legacy_ab=reference_mag))
    fraction = min(s['correction_fraction'] for s in calibrators)
    # Empirical nonoverlapping aperture grid within a 10-40 arcsec annulus. Positive
    # field sources remain in the robust MAD. Check separate quadrants and use the
    # largest scatter to guard against spatial depth variation.
    grid = (x % 10 == 2) & (y % 10 == 2) & (d0 >= 10) & (d0 <= 40)
    usable = fftconvolve(valid.astype(float), kernel, mode='same') >= kernel.sum() - 0.1
    grid &= usable
    background = float(np.median(summed[grid]))
    subsets = [grid & ((x > xp) == east) & ((y > yp) == north)
               for east in (False, True) for north in (False, True)]
    scatters = [float(1.4826 * np.median(np.abs(summed[s] - np.median(summed[s])))) for s in subsets]
    sigma = max(scatters)
    signal = summed - background
    at_centers = [signal[round(float(yy)),round(float(xx))] / sigma for xx,yy in ((xp,yp),(xr,yr))]
    max_sn = [float(np.max(signal[m]) / sigma) for m in center_masks]
    # Report the nearest significant optical peak; inspect radius extends to 15 arcsec.
    peaks = (signal == maximum_filter(signal, size=7)) & (signal >= 5 * sigma) & (d0 <= 15) & usable
    near = []
    for yy, xx in zip(*np.where(peaks)):
        a, d = wcs.pixel_to_world_values(xx, yy)
        near.append(dict(ra=float(a), dec=float(d), distance_native_arcsec=float(d0[yy,xx]),
                         distance_shifted_arcsec=float(d1[yy,xx]), sn=float(signal[yy,xx] / sigma)))
    # Nominal 5-sigma sensitivity and a conservative upper envelope: each possible
    # center has flux <= its measured signal + 5 sigma, conditional on this noise estimate.
    limit = lambda flux: float(zp - 2.5*np.log10(flux / fraction))
    return dict(image=image, weight=weights, filter=band, exptime_header_s=header['EXPTIME'],
                mjd_start_header=header['MJD-OBS'], photzp_ab=zp, plate_scale_arcsec=scale,
                aperture_radius_pixels=4, total_aperture_radius_pixels=15,
                grid_apertures=int(grid.sum()), quadrant_mad_sigmas=scatters,
                weight_range_padded=[float(weight[padded].min()),float(weight[padded].max())],
                calibrators=calibrators, adopted_encircled_fraction=fraction,
                forced_snr_native_shifted=[float(v) for v in at_centers],
                max_snr_native_shifted=max_sn, nearest_5sigma_peaks=sorted(near,key=lambda v:v['distance_native_arcsec'])[:3],
                five_sigma_point_sensitivity_ab=limit(5*sigma),
                conservative_point_limit_over_region_ab=limit(max(0,float(np.max(signal[area])))+5*sigma))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache_directory', type=Path)
    parser.add_argument('output_json', type=Path)
    args = parser.parse_args()
    args.cache_directory.mkdir(parents=True, exist_ok=True)
    out = dict(ra=RA, dec=DEC, comparison_shift_arcsec_east_north=[SHIFT_EAST, SHIFT_NORTH],
               radius_arcsec=R99, measurements=[measure(args.cache_directory, b) for b in ('G.MP9401', 'R.MP9601')])
    args.output_json.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
