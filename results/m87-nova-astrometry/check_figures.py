"""Reproduce the printed cross-marker and jet-angle measurements from the source PDF."""
import csv
import hashlib
import io
from pathlib import Path
import sys
from urllib.request import urlopen

import fitz
import numpy as np
from PIL import Image
from scipy.ndimage import maximum_filter
from scipy.optimize import least_squares
from scipy.signal import convolve2d
from scipy.spatial import cKDTree

from analyze import POSITIONS, ROOT, SkyCoord

URL = 'https://stsci-opo.org/STScI-01J8T4398VD1G2K90XFNS917HW.pdf'
SHA = 'ea343ed770797f021eaf66df22a833c0fe8be3c2d8943a9b719d1e7c14a37bc2'


def embedded_rgb(doc, page_number):
    images = doc[page_number - 1].get_images()
    data = doc.extract_image(images[0][0])['image']
    return np.asarray(Image.open(io.BytesIO(data)).convert('RGB')).astype(float)


def crosses(channel):
    kernel = np.zeros((11, 11))
    for i in range(-4, 5):
        for j in (-1, 0, 1):
            kernel[i + 5, i + 5 + j] = kernel[i + 5, 5 - i + j] = 1
    score = convolve2d(np.maximum(channel - 12, 0), kernel, mode='same')
    return np.argwhere((score == maximum_filter(score, 9)) & (score > 600))[:, ::-1]


def tangent(ra, dec):
    center = SkyCoord('12h30m49.423s +12d23m28.04s')
    return 3600 * np.column_stack(((ra - center.ra.deg) * np.cos(center.dec.radian),
                                    dec - center.dec.deg))


def warp(xy, params):
    return xy @ params[:4].reshape((2, 2)) + params[4:]


def main():
    payload = Path(sys.argv[1]).read_bytes() if len(sys.argv) > 1 else urlopen(URL, timeout=60).read()
    if hashlib.sha256(payload).hexdigest() != SHA:
        raise ValueError('Unexpected PDF version')
    with fitz.open(stream=payload, filetype='pdf') as doc:
        fig1 = embedded_rgb(doc, 2)[:838]  # only the wide upper panel
        fig2 = embedded_rgb(doc, 3)
    red, green, blue = fig1.transpose(2, 0, 1)
    cyan = cKDTree(crosses((green + blue) / 2 - red))
    pink = cKDTree(crosses((red + blue) / 2 - green))
    with (ROOT / 'positions.csv').open() as f:
        catalog = [r for r in csv.DictReader(f) if r['survey'] == 'WFC3']
    xy = tangent(np.array([float(r['ra_deg']) for r in catalog]),
                 np.array([float(r['dec_deg']) for r in catalog]))
    # Initial estimate from north-up/east-left convention and the printed scale bar.
    params = np.array([-3.51, 0, 0, -3.51, 448., 414.])
    for radius in (25, 12, 9, 6, 4, 4):
        distance, idx = cyan.query(warp(xy, params))
        mask = distance < radius
        params = least_squares(lambda q: (warp(xy[mask], q) - cyan.data[idx[mask]]).ravel(),
                               params, loss='soft_l1').x
    for number, positions in POSITIONS.items():
        for label, site in zip(('2024', '2016-adjusted'), positions):
            c = SkyCoord(site)
            point = warp(tangent(c.ra.deg, c.dec.deg), params)[0]
            print(number, label, 'nearest pink cross (pixels)', round(pink.query(point)[0], 2))
    red, green, blue = fig2.transpose(2, 0, 1)
    y, x = np.nonzero((blue - red > 55) & (green - red > 55) & (abs(green - blue) < 40))
    slope = np.polyfit(x, y, 1)[0]
    print('Jet PA', round(270 + np.degrees(np.arctan(-slope)), 3))


if __name__ == '__main__':
    main()
