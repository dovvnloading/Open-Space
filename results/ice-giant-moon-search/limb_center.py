"""Determine native detector planet center independently from the illuminated limb."""
import argparse
import csv
from pathlib import Path
import warnings

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
import numpy as np
from scipy import ndimage, optimize


def measure(path, radius):
    with fits.open(path, memmap=False) as hdus:
        sci = hdus['SCI']
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)
            w = WCS(sci.header, naxis=2)
        center = w.all_world2pix([[sci.header['MT_RA'], sci.header['MT_DEC']]], 0)[0]
        data = sci.data.astype(float)
        good = np.isfinite(data) & ((hdus['DQ'].data & 3) == 0)
    if not np.all((center > radius+7) & (center < np.array(data.shape[::-1]) - radius-7)):
        return
    # DQ holes receive interpolated values for derivatives only.
    smoothed = ndimage.gaussian_filter(np.where(good, data, 0), 1.2)
    coverage = ndimage.gaussian_filter(good.astype(float), 1.2)
    smoothed /= np.maximum(coverage, 1e-6)
    angle = np.linspace(0, 2*np.pi, 360, endpoint=False)
    radial = np.arange(radius-8, radius+8.1, 0.1)
    x = center[0]+np.cos(angle)[:, None]*radial
    y = center[1]+np.sin(angle)[:, None]*radial
    profile = ndimage.map_coordinates(smoothed, [y, x], order=1)
    gradient = np.gradient(profile, radial, axis=1)
    i = np.argmin(gradient, axis=1)
    strength = -gradient[np.arange(len(angle)), i]
    weight = ndimage.map_coordinates(coverage, [y[np.arange(len(i)), i], x[np.arange(len(i)), i]], order=1)
    valid = ((i > 4) & (i < len(radial)-5) & (strength > 0.05*np.median(strength))
             & (weight > .8))
    if valid.sum() < 180:
        raise ValueError(f'insufficient usable limb angles: {path} ({valid.sum()})')
    rx = (radial[i]*np.cos(angle))[valid]
    ry = (radial[i]*np.sin(angle))[valid]
    # R(theta) includes second-harmonic apparent oblateness; robust radial loss.
    # The harmonic approximation assumes the projected disk is near circular.
    def residual(par):
        dx, dy, r, c, s = par
        theta = np.arctan2(ry-dy, rx-dx)
        return np.hypot(rx-dx, ry-dy) - r - c*np.cos(2*theta) - s*np.sin(2*theta)
    fit = optimize.least_squares(residual, [0, 0, radius, 0, 0], loss='soft_l1', f_scale=.25)
    dx, dy, r, c, s = fit.x
    rms = np.sqrt(np.mean(np.minimum(residual(fit.x)**2, 4)))
    return dict(exposure_id=path.stem, x=center[0]+dx, y=center[1]+dy, wcs_x=center[0],
                wcs_y=center[1], dx=dx, dy=dy, radius_px=r, usable_angles=int(valid.sum()),
                clipped_limb_rms_px=rms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('files', nargs='+', type=Path)
    parser.add_argument('--radius-px', type=float, required=True)
    args = parser.parse_args()
    rows = [row for f in sorted(args.files) if (row := measure(f, args.radius_px)) is not None]
    if not rows:
        raise ValueError('no fully enclosed planet limbs')
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
