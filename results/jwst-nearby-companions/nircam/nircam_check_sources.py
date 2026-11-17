"""Inspect all F444W trigger locations with independent references and paired F200W."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from nircam_reduce import (
    FIELDS,
    cube,
    matched,
    median,
    offaxis_psf,
    reduce_cube,
    thresholds,
    write_table,
)
from scipy import ndimage as ndi
from scipy.optimize import least_squares


def geometry(head, shape):
    cd = (
        np.array([[head["CD1_1"], head["CD1_2"]], [head["CD2_1"], head["CD2_2"]]])
        * 3600
    )
    center = np.array([head["CRPIX1"] - 1, head["CRPIX2"] - 1])
    scale = np.sqrt(abs(np.linalg.det(cd)))
    yy, xx = np.indices(shape)
    return cd, center, scale, np.hypot(xx - center[0], yy - center[1]) * scale


def read_short(p):
    with fits.open(p) as h:
        assert h[0].header["FILTER"] == "F200W"
        assert h[0].header["APERNAME"] == "NRCA2_MASK335R"
        a = h["SCI"].data.astype(float)
        a[(h["DQ"].data & 3) != 0] = np.nan
        return a, h["SCI"].header.copy()


def aperture(image, x, y):
    half = 8
    ix, iy = int(round(x)), int(round(y))
    zy, zx = np.mgrid[iy - half : iy + half + 1, ix - half : ix + half + 1]
    if (
        zx.min() < 0
        or zy.min() < 0
        or zx.max() >= image.shape[1]
        or zy.max() >= image.shape[0]
    ):
        return np.nan
    sample = image[zy, zx]
    r = np.hypot(zx - x, zy - y)
    core = r < 3.5
    sky = (r >= 5) & (r < 8)
    if np.isfinite(sample[core]).mean() < 0.82:
        return np.nan
    return (np.nanmean(sample[core]) - median(sample[sky])) * core.sum()


def fit_primary(image, psf, x0, y0):
    yy, xx = np.indices(image.shape)
    sel = np.isfinite(image) & (np.hypot(xx - x0, yy - y0) <= 16)
    image_vals = image[sel]

    def profile(params, x, y):
        return ndi.map_coordinates(
            psf, [18 + y - params[1], 18 + x - params[0]], order=3, mode="constant"
        )

    initial = [x0, y0, max(0.1, image[y0, x0]) / psf.max(), 0.0]

    def diff(p):
        return image_vals - (p[2] * profile(p, xx[sel], yy[sel]) + p[3])

    best = least_squares(
        diff,
        initial,
        bounds=([x0 - 2, y0 - 2, 0, -np.inf], [x0 + 2, y0 + 2, np.inf, np.inf]),
    ).x
    model = best[2] * profile(best, xx, yy)
    good = sel & (np.hypot(xx - best[0], yy - best[1]) < 4)
    return best, model, float(np.corrcoef(image[good], model[good])[0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--archival", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("."))
    args = ap.parse_args()
    refs = [cube(p)[0] for p in sorted(args.archival.glob("jw05835*long_calints.fits"))]
    assert len(refs) == 9
    standards = sorted(args.inputs.glob("jw01537028*long_calints.fits"))
    psf, _ = offaxis_psf(standards[2:])
    detection_psf, _ = offaxis_psf(standards[:2])
    kernel = detection_psf[14:23, 14:23]
    rows = []
    models = []
    for field in ["l788_34", "ross619"]:
        _, src, _ = FIELDS[field]
        a, h, _ = cube(args.inputs / src)
        cd, center, scale, radius = geometry(h, a.shape[1:])
        alt = median(reduce_cube(a, refs, radius), axis=0)
        alt_amp = matched(alt, kernel)
        alt_cut, alt_score, _ = thresholds(alt_amp, radius, center, scale)
        med = np.load(args.out / (field + "_rdi_maps.npz"))["residual_MJy_sr"].astype(
            float
        )
        after = med.copy()
        if field == "l788_34":
            pars, model, cc = fit_primary(med, psf, 131, 221)
            after -= model
            models.append(
                dict(
                    field=field,
                    fit_x0=float(pars[0]),
                    fit_y0=float(pars[1]),
                    model_flux_uJy=float(pars[2] * h["PIXAR_SR"] * 1e12),
                    core_pearson_r=cc,
                )
            )
        residual_amp = matched(after, kernel)
        residual_cut, residual_snr, _ = thresholds(residual_amp, radius, center, scale)
        short, h2 = read_short(args.inputs / src.replace("nrcalong", "nrca2"))
        cd2, center2, scale2, r2 = geometry(h2, short.shape[1:])
        shortrefs = [
            read_short(p)[0]
            for p in args.inputs.glob("jw08826*nrca2*fits")
            if p.name != src.replace("nrcalong", "nrca2")
        ]
        shortmed = median(
            reduce_cube(short, shortrefs, r2, sky_bounds=(3.9, 4.5)), axis=0
        )
        with (args.out / (field + "_peaks.csv")).open() as fh:
            hits = list(csv.DictReader(fh))
        for peak in hits:
            x, y = int(peak["x0"]), int(peak["y0"])
            east, north = cd @ ([x - center[0], y - center[1]])
            dx, dy = np.linalg.solve(cd2, [east, north])
            sx, sy = center2 + [dx, dy]
            rho = np.hypot(dx, dy)
            angles = np.linspace(0, 2 * np.pi, int(2 * np.pi * rho / 8), endpoint=False)
            positions = [
                center2 + rho * np.array([np.cos(p), np.sin(p)]) for p in angles
            ]
            blank = np.array(
                [
                    aperture(shortmed, tx, ty)
                    for tx, ty in positions
                    if np.hypot(tx - sx, ty - sy) > 12
                ]
            )
            sky = np.nanmedian(blank)
            noise = np.nanstd(blank, ddof=1)
            measured = aperture(shortmed, sx, sy) - sky
            scale2jy = h2["PIXAR_SR"] * 1e12
            # Shared fixed detector locations, independent of the reference fit.
            aalt = alt_amp[y, x]
            lalt = alt_score[y, x]
            gen = peak.copy()
            gen.update(
                field=field,
                archival_rdi_snr=float(lalt),
                archival_rdi_flux_uJy=float(aalt * h["PIXAR_SR"] * 1e12),
                archival_above_cut=int(aalt > alt_cut[y, x]),
                one_source_sub_snr=float(residual_snr[y, x]),
                one_source_sub_above_cut=int(residual_amp[y, x] > residual_cut[y, x]),
                f200w_x0=float(sx),
                f200w_y0=float(sy),
                f200w_ap_uJy=float(measured * scale2jy),
                f200w_ap_err_uJy=float(noise * scale2jy),
                f200w_ap_snr=float(measured / noise),
                f200w_n_control=int(np.isfinite(blank).sum()),
            )
            rows.append(gen)
    write_table(args.out / "nircam_trigger_checks.csv", list(rows[0]), rows)
    (args.out / "nircam_psf_fit.json").write_text(json.dumps(models, indent=2) + "\n")
    print(models)
    for r in rows:
        print(
            r["field"],
            r["x0"],
            r["y0"],
            *[
                round(float(r[k]), 2)
                for k in [
                    "archival_rdi_snr",
                    "one_source_sub_snr",
                    "f200w_ap_snr",
                    "f200w_ap_uJy",
                    "f200w_ap_err_uJy",
                ]
            ],
        )


if __name__ == "__main__":
    main()
