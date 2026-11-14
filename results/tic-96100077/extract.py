"""Extract the target clock and difference images from MAST TESSCut FFIs."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from measure import DATA, DEC, OUT, PERIOD, RA, REF, clock
from scipy.interpolate import CubicHermiteSpline

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
AU_LIGHT_SECONDS = 499.00478384


def corrected_times(table, sector):
    # Revert the CCD-centre correction before applying the target-direction correction.
    spacecraft = table["TIME"].astype(float) - table["TIMECORR"].astype(float)
    vectors = pd.read_csv(DATA / f"s{sector:03}-vectors.csv")
    xyz = CubicHermiteSpline(
        vectors.jd - 2457000, vectors[["x", "y", "z"]], vectors[["vx", "vy", "vz"]]
    )(spacecraft)
    ra, dec = np.deg2rad([RA, DEC])
    direction = [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
    if spacecraft.min() < (vectors.jd.min() - 2457000) or spacecraft.max() > (
        vectors.jd.max() - 2457000
    ):
        raise ValueError("Spacecraft ephemeris does not bracket cadences")
    return spacecraft + (xyz @ direction) * AU_LIGHT_SECONDS / 86400


def sector_photometry(sector):
    with fits.open(RAW / f"s{sector:03}.fits") as hdus:
        data = hdus[1].data
        keep = (
            (data["QUALITY"] == 0)
            & np.isfinite(data["TIME"])
            & np.isfinite(data["TIMECORR"])
        )
        keep &= np.isfinite(data["FLUX"]).all(axis=(1, 2))
        data = data[keep]
        t = corrected_times(data, sector)
        cube = data["FLUX"].astype(float)
        wcs = WCS(hdus[2].header)
        tx, ty = np.array(wcs.all_world2pix(RA, DEC, 0), dtype=float)
        y, x = np.indices(cube.shape[1:])
        rim = (x == 0) | (y == 0) | (x == x.max()) | (y == y.max())
        cube -= np.median(cube[:, rim], axis=1)[:, None, None]
        central = (abs(x - round(tx)) <= 1) & (abs(y - round(ty)) <= 1)
        flux = cube[:, central].sum(axis=1)
        receipt = {
            "sector": sector,
            "quality_zero_cadences": len(t),
            "median_net_flux": float(np.median(flux)),
        }
        if (
            receipt["median_net_flux"] < 100
        ):  # e-/s; bright (G=10.88) star must be visible.
            receipt["status"] = "no measurable target in these pixels"
            return receipt, None
        # Nearby times serve as in/out pairs; local detrending protects against scattered light.
        flux = pd.Series(flux, index=t)
        baseline = (
            pd.Series(flux.values, index=pd.to_timedelta(t - t.min(), unit="D"))
            .rolling("1D", center=True, min_periods=10)
            .median()
            .values
        )
        yflux = flux.values / baseline
        valid = np.isfinite(yflux) & (baseline > 100) & (yflux > 0.65) & (yflux < 1.35)
        t, yflux, cube = t[valid], yflux[valid], cube[valid]
        ph = ((t - REF) / PERIOD) % 1
        shift = clock(t, yflux)[0] / (86400 * PERIOD)
        # Per-orbit images limit drift. Difference = quadrature flux minus eclipse flux.
        centered = ((ph - shift + 0.25) % 0.5) - 0.25
        orbit = np.floor((t - REF) / PERIOD - shift + 0.125).astype(int)
        images = []
        for n in np.unique(orbit):
            inside = (orbit == n) & (abs(centered) < 0.065)
            outside = (orbit == n) & (abs(centered) > 0.185)
            if inside.sum() >= 2 and outside.sum() >= 2:
                images.append(
                    np.median(cube[outside], axis=0) - np.median(cube[inside], axis=0)
                )
        im = np.median(images, axis=0)
        # Use a 5x5 central window; faint outer pixels and background gradients bias larger apertures.
        central5 = (abs(x - round(tx)) <= 2) & (abs(y - round(ty)) <= 2)
        weights = np.maximum(im - np.median(im[rim]), 0) * central5
        cx, cy = (
            np.sum(x * weights) / np.sum(weights),
            np.sum(y * weights) / np.sum(weights),
        )
        cra, cdec = wcs.all_pix2world(cx, cy, 0)
        offset = np.hypot((cra - RA) * np.cos(np.deg2rad(DEC)), cdec - DEC) * 3600
        depth = 1 - np.median(yflux[abs(centered) < 0.04]) / np.median(
            yflux[abs(centered) > 0.21]
        )
        receipt.update(
            status="measured",
            usable=len(t),
            n_image_orbits=len(images),
            offset_arcsec=float(offset),
            depth=float(depth),
            target_x=tx,
            target_y=ty,
            deficit_x=cx,
            deficit_y=cy,
        )
        pd.DataFrame(im).to_csv(
            OUT / f"difference-s{sector:03}.csv",
            index=False,
            header=False,
            float_format="%.7f",
        )
        lc = pd.DataFrame({"sector": sector, "time": t, "flux": yflux})
        return receipt, lc


if __name__ == "__main__":
    inventory = json.loads((DATA / "sectors.json").read_text())["results"]
    had = []
    curves = []
    for row in sorted(inventory, key=lambda r: int(r["sector"])):
        receipt, curve = sector_photometry(int(row["sector"]))
        had.append(receipt)
        if curve is not None:
            curves.append(curve)
    pd.DataFrame(had).to_csv(OUT / "sectors.csv", index=False, float_format="%.9f")
    pd.concat(curves).to_csv(RAW / "curves.csv", index=False, float_format="%.12f")
