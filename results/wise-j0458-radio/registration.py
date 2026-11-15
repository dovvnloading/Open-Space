"""Measure relative field astrometry from bright sources in first- and fourth-epoch FITS."""
import argparse
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter
from replay import Image, csv_out


def mad(image):
    return 1.4826 * np.median(abs(image - np.median(image)))


def register(directory):
    endpoint = []
    for epoch in (1, 4):
        path, = directory.glob(f"*VLASS{epoch}*.fits")
        endpoint.append(Image(path))
    first, last = endpoint
    # Find isolated pixel maxima above six local-noise units in both visits.
    peaks = np.argwhere((maximum_filter(first.image, 9) == first.image)
                        & (first.image > 6*mad(first.image)))
    ny, nx = first.image.shape
    results = []
    for y, x in peaks:
        if min(x,y,nx-1-x,ny-1-y)<8 or np.hypot(x-(nx-1)/2,y-(ny-1)/2)<15:
            continue
        # Use a fixed local fitting window and constant background, as in the frozen field table.
        ra1, dec1, flux1 = center(first, x, y)
        xp, yp = np.round(last.wcs.world_to_pixel_values(ra1, dec1)).astype(int)
        if min(xp,yp,nx-1-xp,ny-1-yp)<9:
            continue
        if np.max(last.image[yp-2:yp+3,xp-2:xp+3]) < 6*mad(last.image):
            continue
        ra4, dec4, flux4 = center(last,xp,yp)
        results.append(dict(ra1=ra1,dec1=dec1,amp1=flux1,ra4=ra4,dec4=dec4,amp4=flux4,
            east_shift_as=(ra4-ra1)*3600*np.cos(np.deg2rad(dec1)),north_shift_as=(dec4-dec1)*3600))
    return results


def center(im, x, y):
    from scipy.optimize import least_squares
    use = np.hypot(im.xx-x, im.yy-y) < 6
    dx, dy, pixels = im.xx[use]-x, im.yy[use]-y, im.image[use]
    def residual(p):
        amplitude, xo, yo, baseline = p
        return amplitude*im.kernel(dx-xo,dy-yo)+baseline-pixels
    res = least_squares(residual,[max(pixels),0,0,0],bounds=([0,-2.5,-2.5,-2],[np.inf,2.5,2.5,2]))
    amplitude, xo, yo = res.x[:3]
    ra, dec = im.wcs.pixel_to_world_values(x+xo,y+yo)
    return float(ra), float(dec), float(amplitude)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    csv_out(args.output, register(args.directory))
