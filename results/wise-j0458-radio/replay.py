"""Fit the VLASS restoring beam at moving, stationary and matched comparison positions."""
import csv
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
RA, DEC, MJD, PMRA, PMDEC = 74.7268, 64.5825, 60859., 218., 300.


def predicted(date):
    years = (Time(date).mjd - MJD) / 365.25
    return RA + PMRA * years / 3.6e6 / np.cos(np.deg2rad(DEC)), DEC + PMDEC * years / 3.6e6


class Image:
    def __init__(self, path):
        with fits.open(path) as f:
            self.header = f[0].header.copy()
            self.image = np.asarray(f[0].data.squeeze(), dtype=float) * 1e3
        self.wcs = WCS(self.header).celestial
        self.date = self.header["DATE-OBS"]
        self.yy, self.xx = np.indices(self.image.shape)
        self.jac = self.wcs.pixel_scale_matrix * 3600
        self.pa = np.deg2rad(self.header["BPA"])
        self.axes = np.array([self.header["BMAJ"], self.header["BMIN"]]) * 3600

    def kernel(self, dx, dy):
        east, north = self.jac @ np.array([dx, dy])
        u = (east*np.sin(self.pa) + north*np.cos(self.pa)) / self.axes[0]
        v = (east*np.cos(self.pa) - north*np.sin(self.pa)) / self.axes[1]
        return np.exp(-4*np.log(2)*(u*u + v*v))

    def fit(self, sky, offset=(0., 0.)):
        x, y = self.wcs.world_to_pixel_values(*sky)
        dx, dy = self.xx-x-offset[0], self.yy-y-offset[1]
        use = (abs(dx)<9) & (abs(dy)<9) & np.isfinite(self.image)
        dx, dy = dx[use], dy[use]
        basis = np.column_stack([self.kernel(dx, dy), np.ones_like(dx), dx, dy])
        return np.linalg.lstsq(basis, self.image[use], rcond=None)[0][0]

    def centroid(self, sky):
        x, y = self.wcs.world_to_pixel_values(*sky)
        dx, dy = self.xx-x, self.yy-y
        take = (np.hypot(dx, dy)<8) & np.isfinite(self.image)
        dx, dy, pixels = dx[take], dy[take], self.image[take]
        def residual(p):
            a, xo, yo, bg, gx, gy = p
            return a*self.kernel(dx-xo,dy-yo) + bg + gx*dx + gy*dy - pixels
        sol = least_squares(residual, [.6,0,0,0,0,0],bounds=([0,-3,-3,-2,-1,-1],[10,3,3,2,1,1]))
        a, xo, yo = sol.x[:3]
        e, n = self.jac @ [xo,yo]
        rc, dc = self.wcs.pixel_to_world_values(x+xo,y+yo)
        return float(rc), float(dc), a, e, n


def csv_out(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    images = [Image(p) for p in sorted((ROOT/"inputs").glob("VLASS*.fits"))]
    positions = [predicted(im.date) for im in images]
    measurements, controls = [], []
    for epoch, (im, pos) in enumerate(zip(images, positions), 1):
        offsets = [(a,b) for b in range(-32,33,6) for a in range(-32,33,6) if 7<np.hypot(a,b)<35]
        amplitudes = np.array([im.fit(pos, xy) for xy in offsets])
        med = np.median(amplitudes)
        noise = 1.4826*np.median(np.abs(amplitudes-med))
        a = im.fit(pos)
        r, d, free, east, north = im.centroid(pos)
        measurements.append(dict(epoch=epoch, date=im.date, ra=pos[0], dec=pos[1], forced_mjy=a,
            control_median_mjy=med, control_mad_mjy=noise, z=(a-med)/noise, rank=1+sum(amplitudes>=a),
            at_2017_coord_mjy=im.fit(positions[0]), at_2025_coord_mjy=im.fit(positions[-1]),
            centroid_ra=r if free>.1 else None, centroid_dec=d if free>.1 else None, free_peak_mjy=free,
            east_offset_as=east if free>.1 else None, north_offset_as=north if free>.1 else None))
        controls.extend(dict(epoch=epoch, dx=xy[0], dy=xy[1], forced_mjy=v) for xy,v in zip(offsets,amplitudes))
    csv_out(ROOT/"epochs.csv", measurements)
    csv_out(ROOT/"controls.csv", controls)
    for r in measurements:
        print(r["epoch"],*(round(r[k],3) for k in ("forced_mjy","control_mad_mjy","z",
                                                   "at_2017_coord_mjy","at_2025_coord_mjy")))
    amps = np.array([r["forced_mjy"] for r in measurements[:3]])
    weights = np.array([r["control_mad_mjy"]**-2 for r in measurements[:3]])
    comparison = np.array([r["forced_mjy"] for r in controls]).reshape(4,-1)[:3]
    comparison = np.average(comparison,axis=0,weights=weights)
    print("Prior-only inverse-variance stack, control SD",np.average(amps,weights=weights),comparison.std())


if __name__ == "__main__":
    main()
