from pathlib import Path
import numpy as np, json, csv
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.signal import fftconvolve
from astropy.stats import sigma_clipped_stats
import sys, warnings

warnings.filterwarnings("ignore")
for arg in sys.argv[1:]:
    p = Path(arg)
    m = json.loads(p.with_suffix(".json").read_text())
    coef = 2.6
    h = fits.open(p)
    for hh in h:
        if hh.data is None:
            continue
        a = hh.data.astype(float)
        ylen, xlen = a.shape
        w = WCS(hh.header)
        x, y = w.all_world2pix([[float(m["Object_RA"]), float(m["Object_Dec"])]], 0)[0]
        if 12 < x < xlen - 12 and 12 < y < ylen - 12:
            break
    else:
        continue
    bg = np.nanmedian(a)
    noise = 1.4826 * np.nanmedian(abs(a - bg))
    a = np.where(np.isfinite(a), a, bg)
    a -= bg
    sm = gaussian_filter(a, 1.2)
    sy, sx = np.where((sm == maximum_filter(sm, 17)) & (sm > 6 * noise) & (a < 45000))
    y0, x0 = round(y), round(x)
    idx = np.indices((21, 21)) - 10
    Y, X = idx
    wcs = w.all_pix2world([[x, y], [x + 1, y], [x, y + 1]], 0)
    mat = (wcs[1:] - wcs[0]).T * 3600
    mat[0] *= np.cos(np.deg2rad(float(m["Object_Dec"])))

    def moments(s):
        if s.shape != (21, 21):
            return
        cx = cy = 0
        for i in range(7):
            weight = np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * coef**2))
            I = s * weight
            if I.sum() <= 0:
                return np.full((2, 2), np.nan)
            cx = np.sum(I * X) / I.sum()
            cy = np.sum(I * Y) / I.sum()
        I = s * weight
        a = X - cx
        b = Y - cy
        return (
            mat
            @ np.array(
                [
                    [np.sum(a * a * I), np.sum(a * b * I)],
                    [np.sum(a * b * I), np.sum(b * b * I)],
                ]
            )
            / I.sum()
            @ mat.T
        )

    target = moments(a[y0 - 10 : y0 + 11, x0 - 10 : x0 + 11])
    field = []
    for yy, xx in zip(sy, sx):
        if min(xx, yy, xlen - xx, ylen - yy) < 11 or np.hypot(xx - x, yy - y) < 12:
            continue
        if a[yy, xx] / noise < 12:
            continue
        u = moments(a[yy - 10 : yy + 11, xx - 10 : xx + 11])
        field.append(u)
    field = np.array(field)
    if not len(field):
        continue
    F = np.nanmedian(field, axis=0)
    v, M = np.linalg.eigh(target - F)
    ve = M[:, -1]
    print(
        p.stem,
        round(np.rad2deg(np.arctan2(*ve)) % 180, 2),
        round(float(v[-1] - v[0]), 4),
        "n",
        len(field),
    )
