"""Measure narrow-band residuals in the included WCS-preserving HST cutouts."""
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import convolve, gaussian_filter, map_coordinates, maximum_filter
from scipy.optimize import least_squares

HERE = Path(__file__).parent


def stamp(filename):
    with np.load(HERE / (filename + 'cut.npz')) as data:
        wcs = WCS(fits.Header.fromstring(str(data['wcs']), sep=''))
        return data['arr'].astype(float), wcs, data['np0'], data['xy'] - data['np0']


def measurement(galaxy, pn_coord, single_filter=False):
    narrow, wcs, origin, center = stamp(galaxy + 'narrow')
    yy, xx = np.indices(narrow.shape)
    world = wcs.all_pix2world(xx + origin[0], yy + origin[1], 0)
    pn = np.array(wcs.all_world2pix(*pn_coord, 0)) - origin
    radius = np.hypot(xx - center[0], yy - center[1])
    controls = (radius > 9) & (radius < 32)
    layers = [np.ones(narrow.shape)]
    for band in (['f555w'] if single_filter else ['f555w', 'f814w']):
        continuum, cwcs, corigin, _ = stamp(('ngc628e' if single_filter else galaxy) + band)
        cx, cy = cwcs.all_world2pix(*world, 0)
        layers.append(map_coordinates(continuum, [cy-corigin[1], cx-corigin[0]], order=1, mode='nearest'))
    design = np.stack(layers, axis=-1)
    noise = 1.4826*np.median(abs(narrow[controls]-np.median(narrow[controls])))
    n = len(layers)-1
    coefficients = least_squares(
        lambda beta: (design[controls] @ beta - narrow[controls])/noise,
        [0]+[.03]*n, loss='soft_l1', bounds=([-np.inf]+[0]*n, [np.inf]*(n+1))).x
    residual = narrow - design @ coefficients
    peaks = gaussian_filter(residual, 1)
    peaks[np.hypot(xx-pn[0], yy-pn[1])>5] = -np.inf
    iy, ix = np.unravel_index(peaks.argmax(), peaks.shape)
    ky,kx = np.mgrid[-3:4, -3:4]
    aperture_map = convolve(residual, (np.hypot(kx,ky)<2.5).astype(float))
    obs = aperture_map[iy,ix]
    bg = aperture_map[controls]
    local_best = maximum_filter(aperture_map, footprint=np.hypot(*np.mgrid[-5:6,-5:6])<=5)[controls]
    return dict(residual_sum=float(obs), robust_aperture_z=float((obs-np.median(bg))/(1.4826*np.median(abs(bg-np.median(bg))))),
                local_fraction_with_larger_searched_aperture=float(np.mean(local_best>=obs)),
                peak_to_cluster_arcsec=float(.04*np.hypot(ix-center[0],iy-center[1])))


if __name__ == '__main__':
    targets = json.loads((HERE/'targets.json').read_text())
    print(json.dumps({g:measurement(g,x['pn_radec_deg'],g=='ngc628') for g,x in targets.items()},indent=2))
