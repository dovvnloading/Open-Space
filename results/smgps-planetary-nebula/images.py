"""Register the radio image and two optical surveys on the published MASH position."""

import csv
import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
CENTER = SkyCoord(145.505 * u.deg, -52.345 * u.deg)
PIXEL = 1.5
GRID = np.arange(-120, 120.01, PIXEL)
EAST, NORTH = np.meshgrid(GRID, GRID)
R = np.hypot(EAST, NORTH)
SKY = CENTER.spherical_offsets_by(EAST * u.arcsec, NORTH * u.arcsec)


def sample(path):
    layers = []
    with fits.open(path) as hdus:
        for h in hdus:
            if not h.is_image or h.data is None or h.data.ndim != 2:
                continue
            header = h.header.copy()
            if path.name in ("ha.fits", "sr.fits"):
                if h is not hdus[0]:
                    continue
                # The digitized plates carry an additional legacy DSS solution.
                # Select the explicit J2000 TAN/CD calibration of the cutout.
                header = fits.Header(
                    {
                        k: header[k]
                        for k in [
                            "NAXIS",
                            "NAXIS1",
                            "NAXIS2",
                            "CTYPE1",
                            "CTYPE2",
                            "CRPIX1",
                            "CRPIX2",
                            "CRVAL1",
                            "CRVAL2",
                            "CD1_1",
                            "CD1_2",
                            "CD2_1",
                            "CD2_2",
                            "EQUINOX",
                        ]
                    }
                )
            w = WCS(header).celestial
            x, y = w.world_to_pixel(SKY)
            image = h.data.astype(float)
            layers.append(
                map_coordinates(
                    image,
                    [y, x],
                    order=1,
                    mode="constant",
                    cval=np.nan,
                    prefilter=False,
                )
            )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(layers, axis=0)


def smooth(img, sigma):
    good = np.isfinite(img)
    num = gaussian_filter(np.where(good, img, 0), sigma)
    den = gaussian_filter(good.astype(float), sigma)
    return np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0.75)


def suppress_stars(img):
    # Common angular filter on every band before matching the restoring beam.
    val = np.nanmedian(img)
    med = median_filter(np.nan_to_num(img, nan=val), size=5)
    med[~np.isfinite(img)] = np.nan
    return med


def residual(ha, red):
    good = (R > 95) & np.isfinite(ha + red)
    h0 = np.median(ha[good])
    r0 = np.median(red[good])
    h = ha - h0
    r = red - r0
    # Match compact continuum structure on the outskirts using unsaturated stars.
    lim = np.percentile(r[good], [65, 98])
    star = good & (r > lim[0]) & (r < lim[1])
    scale = np.sum(h[star] * r[star]) / np.sum(r[star] ** 2)
    return suppress_stars(h - scale * r), float(scale)


def stack(paths):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian([sample(p) for p in paths], axis=0)


def display(image):
    lo, hi = np.nanpercentile(image[(R < 109) & np.isfinite(image)], [12, 96])
    return np.arcsinh(np.clip((image - lo) / (hi - lo), 0, 1.4) * 2)


def main():
    radio = sample(CACHE / "radio.fits") * 1.0e6
    blue = sample(CACHE / "ha.fits")
    red = sample(CACHE / "sr.fits")
    old, scale1 = residual(blue, red)
    paths = sorted(CACHE.glob("ADP*.fits"))
    ha = [p for p in paths if fits.getheader(p)["FILTER"] == "NB_659"]
    rd = [p for p in paths if fits.getheader(p)["FILTER"] == "r_SDSS"]
    new, scale2 = residual(stack(ha), stack(rd))
    # Blur optical details to roughly the 8 arcsec radio resolution. Cropping and
    # median replacement suppress point-like stars and the bright radio interloper.
    radio = smooth(suppress_stars(radio), 1)
    shs = smooth(old, 2.1)
    vphas = smooth(new, 2.1)
    good = (R < 96) & np.isfinite(radio + shs + vphas)

    def rho(a, b, mask=good):
        ok = mask & np.isfinite(a + b)
        return float(spearmanr(a[ok], b[ok]).statistic)

    # Compare sectors within concentric 12 arcsec annuli (removes generic round morphology).
    def beyond_radius(a):
        out = a.copy()
        for lo in np.arange(0, 108, 12):
            ann = (R >= lo) & (R < lo + 12)
            out[ann] -= np.nanmedian(a[ann])
        return out

    # Rotations on the original square grid about the PN position. No fitting of the center/angle.
    angle = np.linspace(0, 360, 72, endpoint=False)
    controls = []
    for deg in angle:
        co, si = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
        yy = (+si * EAST + co * NORTH - GRID[0]) / PIXEL
        xx = (co * EAST - si * NORTH - GRID[0]) / PIXEL
        for label, arr in [("SHS", shs), ("VPHAS", vphas)]:
            shifted = map_coordinates(
                beyond_radius(arr), [yy, xx], order=1, cval=np.nan
            )
            controls.append(
                dict(survey=label, deg=int(deg), rho=rho(beyond_radius(radio), shifted))
            )
    with (HERE / "rotation-controls.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=controls[0])
        w.writeheader()
        w.writerows(controls)
    results = dict(
        grid_arcsec=PIXEL,
        radio_shs=rho(radio, shs),
        radio_vphas=rho(radio, vphas),
        shs_vphas=rho(shs, vphas),
        area_pixels=int(good.sum()),
        scale_shs=scale1,
        scale_vphas=scale2,
        ha_exposures=len(ha),
        r_exposures=len(rd),
    )
    # Quadrants summarized in absolute sky orientation around the catalog position.
    quad = []
    for name, a in [("MeerKAT", radio), ("SHS", shs), ("VPHAS", vphas)]:
        for ns, ys in [("N", 1), ("S", -1)]:
            for ew, xs in [("E", 1), ("W", -1)]:
                sel = (NORTH * ys > 0) & (EAST * xs > 0) & (R >= 42) & (R <= 88)
                quad.append(
                    dict(
                        survey=name,
                        quadrant=ns + ew,
                        median=float(np.nanmedian(a[sel])),
                    )
                )
    (HERE / "morphology.json").write_text(
        json.dumps({**results, "rim_quadrants": quad}, indent=2) + "\n"
    )
    fig, axs = plt.subplots(1, 3, figsize=(10.6, 3.9), layout="constrained")
    for ax, img, title in zip(
        axs,
        [radio, shs, vphas],
        ["MeerKAT · 2018", "UKST Hα − red · 1999", "VPHAS+ Hα − red · 2014 / 2016"],
    ):
        ax.imshow(
            display(img),
            origin="lower",
            extent=[GRID[0] - 0.75, GRID[-1] + 0.75, GRID[0] - 0.75, GRID[-1] + 0.75],
            cmap="gray_r",
        )
        ax.set(
            xlim=(113, -113),
            ylim=(-113, 113),
            title=title,
            xlabel="East offset (arcsec)",
            xticks=[90, 45, 0, -45, -90],
            yticks=[-90, -45, 0, 45, 90],
        )
        ax.plot(0, 0, marker="+", ms=6, color="#bd4422", mew=0.8)
    axs[0].set_ylabel("North offset (arcsec)")
    fig.savefig(HERE / "registered-shells.png", dpi=155)
    plt.close(fig)
    print(results)


if __name__ == "__main__":
    main()
