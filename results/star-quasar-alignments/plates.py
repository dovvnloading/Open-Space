"""Compare historical DSS centroids with Gaia predictions on locally registered plates."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent


def centroid(image, guess, radius):
    cx, cy = np.round(guess).astype(int)
    patch = image[cy-radius:cy+radius+1, cx-radius:cx+radius+1]
    yy, xx = np.mgrid[-radius:radius+1, -radius:radius+1]
    floor = np.median(patch)
    amp = patch.max()-floor
    if amp < 100:
        return np.full(2, np.nan)

    def residual(v):
        dx, dy, height, sigma, background = v
        model = background + height*np.exp(-((xx-dx)**2+(yy-dy)**2)/(2*sigma**2))
        return (model - patch).ravel()

    fit = least_squares(residual, [*(guess-[cx,cy]),amp,1.4,floor],
                        bounds=([-2.5,-2.5,0,.5,0],[2.5,2.5,2*amp,3.8,np.inf]),
                        loss="soft_l1", f_scale=150)
    return np.array([cx, cy]) + fit.x[:2]


def sky_to_plate(wcs, stars, dt):
    east = dt*stars.pmra / (3.6e6*np.cos(np.radians(stars.dec)))
    north = dt*stars.pmdec / 3.6e6
    return np.column_stack(wcs.all_world2pix(stars.ra+east, stars.dec+north, 0))


def measure(image, header, stars, targets):
    wcs = WCS(header)
    year = Time(header["DATE-OBS"][:10]).jyear
    dt = year - 2016
    pos = sky_to_plate(wcs, stars, dt)
    nearest = cKDTree(pos).query(pos, k=2)[0][:,1]
    x, y = pos.T
    h, w = image.shape
    use = ((x>7) & (y>7) & (x<w-8) & (y<h-8) & (nearest>8)
           & (stars.phot_g_mean_mag>13) & (stars.phot_g_mean_mag<18)
           & ~stars.source_id.isin(targets.source_id))
    pred = pos[use]
    measured = np.array([centroid(image, point, 4) for point in pred])
    valid = np.isfinite(measured[:,0])
    measured, pred = measured[valid], pred[valid]
    design = np.column_stack((pred, np.ones(len(pred))))
    keep = np.ones(len(pred), bool)
    for _ in range(4):
        affine = np.linalg.lstsq(design[keep], (measured-pred)[keep], rcond=None)[0]
        diff = np.linalg.norm(measured-pred-design@affine, axis=1)
        keep = diff < max(.35, 3*np.median(diff))
    scale = 3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
    results = []
    for row, pix in zip(targets.itertuples(), sky_to_plate(wcs, targets, dt)):
        center = pix + np.r_[pix, 1] @ affine
        for radius in (2, 3, 4):
            position = centroid(image, center, radius)
            coords = wcs.all_pix2world([position, center], 0)
            delta = (coords[0]-coords[1])*[np.cos(np.radians(row.dec))*3600, 3600]
            results.append([header["PLATEID"], year, row.source_id, radius, int(keep.sum()),
                            np.median(diff[keep])*scale, *delta])
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    cache = parser.parse_args().cache
    cache.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(HERE / "gaia-2109.csv")
    catalog = pd.read_csv(HERE / "plate-reference.csv").dropna(subset=["pmra","pmdec"])
    rows = []
    for item in json.loads((HERE/"sources.json").read_text()):
        if not item["file"].startswith("poss"):
            continue
        path = cache/item["file"]
        if not path.exists():
            with urlopen(item["url"], timeout=120) as stream:
                path.write_bytes(stream.read())
        with fits.open(path) as h:
            im = h[0].data
            if hashlib.sha256(im.tobytes()).hexdigest() != item["pixels_sha256"]:
                raise ValueError(f"Pixel data changed for {path}")
            rows.extend(measure(im, h[0].header, catalog, targets))
    cols = ["plate","year","id","radius_pix","controls","control_median_arcsec",
            "east_arcsec","north_arcsec"]
    result = pd.DataFrame(rows, columns=cols)
    expected = pd.read_csv(HERE/"plate-offsets.csv").sort_values(["plate","id","radius_pix"])
    result = result.sort_values(["plate","id","radius_pix"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(result, expected.reset_index(drop=True), atol=.0001, rtol=0)
    print(result[result.radius_pix==3].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
