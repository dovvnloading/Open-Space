"""Compare registered L 788-34 cutouts, aperture flux, and empirical point sources."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.colors import AsinhNorm
from nircam_check_sources import geometry, read_short
from nircam_reduce import FIELDS, cube, median, reduce_cube, write_table
from scipy.ndimage import map_coordinates


def growth(im, pos, scale, rad, sky=(1.16, 1.5)):
    x, y = pos
    n = int(np.ceil(sky[1] / scale)) + 1
    ix, iy = round(x), round(y)
    if min(ix, iy) < n or ix + n >= im.shape[1] or iy + n >= im.shape[0]:
        return np.nan
    yy, xx = np.mgrid[iy - n : iy + n + 1, ix - n : ix + n + 1]
    z = im[yy, xx]
    rho = np.hypot(xx - x, yy - y) * scale
    core = rho < rad
    background = (rho >= sky[0]) & (rho < sky[1])
    if np.isfinite(z[core]).mean() < 0.8:
        return np.nan
    return core.sum() * (np.nanmean(z[core]) - median(z[background]))


def standard(path):
    with fits.open(path) as hd:
        dat = hd["SCI"].data.astype(float)
        dat[(hd["DQ"].data & 3) != 0] = np.nan
        im = median(dat, axis=0)
        y, x = np.unravel_index(np.nanargmax(im), im.shape)
        assert min(y, x, im.shape[0] - y, im.shape[1] - x) > 50
        _, _, scale, _ = geometry(hd["SCI"].header, im.shape)
        yy, xx = np.mgrid[y - 1 : y + 2, x - 1 : x + 2]
        values = im[yy, xx] - np.nanmin(im[yy, xx])
        pos = np.array([np.nansum(xx * values), np.nansum(yy * values)]) / np.nansum(
            values
        )
    return im, pos, scale


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("."))
    args = ap.parse_args()
    _, f, _ = FIELDS["l788_34"]
    swfiles = sorted(args.inputs.glob("jw08826*nrca2*fits"))
    a, h, _ = cube(args.inputs / f)
    cd, cen, scale, _ = geometry(h, a.shape[1:])
    m = np.load(args.out / "l788_34_rdi_maps.npz")["residual_MJy_sr"].astype(float)
    with (args.out / "nircam_psf_fit.json").open() as ffit:
        fit = json.load(ffit)[0]
    pos = np.array([fit["fit_x0"], fit["fit_y0"]])
    sw, f2 = read_short(args.inputs / f.replace("nrcalong", "nrca2"))
    cd2, c2, scale2, r2 = geometry(f2, sw.shape[1:])
    pos2 = c2 + np.linalg.solve(cd2, cd @ (pos - cen))
    sr = [read_short(p)[0] for p in swfiles if p.name != f.replace("nrcalong", "nrca2")]
    m2 = median(reduce_cube(sw, sr, r2, sky_bounds=(3.9, 4.5)), axis=0)
    s1 = sorted(args.inputs.glob("jw01537026*_nrca2_calints.fits"))
    s2 = sorted(args.inputs.glob("jw01537028*_nrcalong_calints.fits"))
    assert len(s1) == 3 and len(s2) == 5
    bands = [
        ("F200W", m2, pos2, c2, scale2, f2, s1),
        ("F444W", m, pos, cen, scale, h, s2),
    ]
    rows = []
    fig, axs = plt.subplots(2, 3, figsize=(9.6, 5.5), layout="constrained")
    for index, (filt, image, position, center, scl, header, standards) in enumerate(
        bands
    ):
        psfs = [standard(p) for p in standards]
        rap_arcsec = np.array([0.105, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 1.1])
        sep = np.hypot(*(position - center)) * scl
        angles = np.linspace(0, 2 * np.pi, 36, endpoint=False)
        locs = [
            center + sep / scl * np.array([np.cos(th), np.sin(th)]) for th in angles
        ]
        growth_object = []
        growth_point = []
        for rap in rap_arcsec:
            scale_flux = 1e12 * header["PIXAR_SR"]
            raw = growth(image, position, scl, rap) * scale_flux
            samples = np.array(
                [
                    growth(image, l, scl, rap) * scale_flux
                    for l in locs
                    if np.hypot(*(l - position)) * scl >= 3.0
                ]
            )
            samples = samples[np.isfinite(samples)]
            bg = median(samples)
            n = 1.4826 * median(np.abs(samples - bg))
            ee = median(
                [
                    growth(z, loc, sc, rap) / growth(z, loc, sc, 1.1)
                    for z, loc, sc in psfs
                ]
            )
            rows.append(
                dict(
                    band=filt,
                    r_arcsec=rap,
                    flux_uJy=raw - bg,
                    local_blank_median_uJy=bg,
                    blank_robust_sigma_uJy=n,
                    n_overlapping_blank_apertures=len(samples),
                    psf_flux_fraction_relative_1p1=ee,
                )
            )
            growth_object.append(raw - bg)
            growth_point.append(ee)
        for ax, im, loc, pix in [
            (axs[index, 0], image, position, scl),
            (axs[index, 1], *psfs[0]),
        ]:
            y, x = np.mgrid[-1.45:1.46:0.025, -1.45:1.46:0.025]
            from scipy.ndimage import distance_transform_edt

            indices = distance_transform_edt(
                ~np.isfinite(im), return_distances=False, return_indices=True
            )
            out = map_coordinates(
                im[tuple(indices)], [y / pix + loc[1], x / pix + loc[0]], order=1
            )
            out -= median(out[np.hypot(x, y) > 1.25])
            peak = np.nanpercentile(out, 99.98)
            out /= peak
            ax.imshow(
                out,
                norm=AsinhNorm(0.05, vmin=-0.02, vmax=1),
                cmap="magma",
                origin="lower",
                extent=[-1.45, 1.45, -1.45, 1.45],
            )
            ax.set_xticks([-1, 0, 1])
            ax.set_yticks([-1, 0, 1])
        ax = axs[index, 2]
        ax.plot(
            rap_arcsec,
            np.array(growth_object) / growth_object[-1],
            "o-",
            color="#075989",
            label="Source",
            ms=4,
        )
        ax.plot(
            rap_arcsec, growth_point, "s--", color="#9c4815", label="Calibrator", ms=3
        )
        ax.set_ylim(0, 1.07)
        ax.grid(alpha=0.22)
        ax.set_ylabel("Fraction inside 1.1″")
        axs[index, 0].set_ylabel(filt + " offset (arcsec)")
        ax.legend(frameon=False, fontsize=8)
    axs[0, 0].set_title("L 788-34 field source")
    axs[0, 1].set_title("Off-axis standard")
    axs[0, 2].set_title("Aperture growth")
    axs[1, 0].set_xlabel("Detector offset (arcsec)")
    axs[1, 1].set_xlabel("Detector offset (arcsec)")
    axs[1, 2].set_xlabel("Aperture radius (arcsec)")
    fig.savefig(args.out / "l788_34_morphology.png", dpi=150)
    write_table(args.out / "l788_34_growth.csv", list(rows[0]), rows)
    for r in rows:
        if r["r_arcsec"] in (0.2, 0.9, 1.1):
            print(r)


if __name__ == "__main__":
    main()
