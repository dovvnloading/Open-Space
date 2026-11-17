"""Calibrated MIRI F1065C RDI and injection recovery on HIP 25878 (GO 9056)."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi
from scipy.optimize import lsq_linear
from scipy.stats import norm, t

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nircam"))
from nircam_reduce import fill, fit_at, matched, median, write_table
from standard_photometry import band_jy

INNER = 3.2
OUTER = 6.2
TRANSMISSION_FLOOR = 0.90
RADII = [3.45, 3.95, 4.45, 4.95, 5.45, 5.95]


def load(path):
    with fits.open(path) as hd:
        assert (
            hd[0].header["FILTER"] in ("F1065C", "F1550C")
            and hd["SCI"].header["BUNIT"] == "MJy/sr"
        )
        z = hd["SCI"].data.astype(float)
        z[(hd["DQ"].data & 3) != 0] = np.nan
        return z, hd["SCI"].header.copy(), hd[0].header.copy()


def calibrators(paths):
    stamps = []
    centers = []
    diagnostics = []
    for p in paths:
        data, sh, head = load(p)
        assert head["EXP_TYPE"] == "MIR_CORONCAL"
        z = fill(data)
        py, px = np.indices(data.shape)
        guess_x = sh["CRPIX1"] - 1 + head["XOFFSET"] / 0.11
        guess_y = sh["CRPIX2"] - 1 + head["YOFFSET"] / 0.11
        nearby = np.hypot(px - guess_x, py - guess_y) < 14
        y, x = np.unravel_index(
            np.nanargmax(np.where(nearby, data, np.nan)), data.shape
        )
        core = z[y - 2 : y + 3, x - 2 : x + 3]
        dy, dx = ndi.center_of_mass(core - np.min(core))
        y += dy - 2
        x += dx - 2
        yy, xx = np.indices(data.shape)
        r = np.hypot(xx - x, yy - y)
        z -= median(z[(r >= 18) & (r < 24) & np.isfinite(data)])
        flux = np.sum(z[r < 16])
        cy, cx = np.mgrid[-15:16, -15:16]
        band = head["FILTER"]
        std_jy = band_jy(p.parent / "hd167060_mod_006.fits", p.parent / (band + ".vot"))
        # Divide out the calibrator's own modeled transmission because the
        # trial location receives its absolute relative-to-clear map value.
        tmap = fits.getdata(next(p.parent.glob("JWST_MIRI_" + band + "*fits")), 0)
        tcal = float(ndi.map_coordinates(tmap, [[y], [x]], order=1)[0])
        total_units = std_jy / (sh["PIXAR_SR"] * 1e6)
        diagnostics.append(
            dict(
                filename=p.name,
                model_jy=std_jy,
                aperture16_over_model=float(flux / total_units),
                calibration_transmission=tcal,
            )
        )
        stamp = ndi.map_coordinates(z, [cy + y, cx + x], order=3, mode="nearest") / (
            total_units * tcal
        )
        # This finite support is outside the RDI fit for every injection.
        stamp[np.hypot(cx, cy) > 12] = 0.0
        stamps.append(stamp)
        centers.append([float(x), float(y)])
    return median(stamps[:2], axis=0), median(stamps[2:], axis=0), centers, diagnostics


def rdi(science, refs, radius):
    fit = (radius > 0.36) & (radius < 1.55)
    mat = np.stack([fill(x) for x in refs] + [np.ones_like(radius)], axis=-1)
    out = []
    for img in science:
        ok = fit & np.isfinite(img)
        # Positive RDI coefficients keep donor speckles and noise bounded.
        bounds = (np.array([0] * len(refs) + [-np.inf]), np.inf)
        for _ in range(3):
            fitresult = lsq_linear(mat[ok], img[ok], bounds=bounds)
            if not fitresult.success:
                raise RuntimeError(fitresult.message)
            diff = img - mat @ fitresult.x
            center = median(diff[ok])
            noise = 1.4826 * median(abs(diff[ok] - center))
            ok = fit & np.isfinite(img) & (abs(diff - center) < 4 * noise)
        out.append(diff)
    return np.array(out)


def levels(amp, rad, center, scale, usable, beam_px):
    cut = np.full_like(amp, np.nan)
    snr = np.full_like(amp, np.nan)
    stats = []
    for lo in np.arange(INNER, OUTER - 0.01, 0.5):
        hi = lo + 0.5
        rp = (lo + hi) / 2 / scale
        n = int(2 * np.pi * rp / beam_px)
        pa = np.arange(n) * 2 * np.pi / n
        x = np.rint(center[0] + rp * np.cos(pa)).astype(int)
        y = np.rint(center[1] + rp * np.sin(pa)).astype(int)
        blank = amp[y, x][usable[y, x] & np.isfinite(amp[y, x])]
        n = len(blank)
        mu = np.mean(blank)
        sd = np.std(blank, ddof=1) * np.sqrt(1 + 1 / n)
        threshold = mu + t.isf(norm.sf(5), n - 1) * sd
        selected = (rad >= lo) & (rad < hi) & usable
        cut[selected] = threshold
        snr[selected] = (amp[selected] - mu) / sd
        stats.append(
            dict(
                lo_arcsec=lo,
                hi_arcsec=hi,
                n_control=n,
                mean=mu,
                sigma=sd,
                threshold=threshold,
            )
        )
    return cut, snr, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("."))
    ap.add_argument("--host-jy", type=float, required=True)
    ap.add_argument("--band", choices=["F1065C", "F1550C"], default="F1065C")
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    specs = {
        "F1065C": ("38", "40", "35", "33", "25", 3.5),
        "F1550C": ("37", "39", "36", "34", "27", 4.8),
    }
    sci_obs, sci_sky, ref_obs, bkg_obs, calobs, beam_px = specs[args.band]
    prefix = args.band.lower() + "_"

    def exposures(prefix):
        paths = sorted(args.inputs.glob(prefix + "*fits"))
        return [load(p)[0] for p in paths]

    (sci_file,) = args.inputs.glob("jw090560" + sci_obs + "*fits")
    sci, h, h0 = load(sci_file)
    target_bkg = exposures("jw090560" + sci_sky)
    ref_bkg = exposures("jw090560" + bkg_obs)
    assert len(target_bkg) == len(ref_bkg) == 2
    sky = median(np.concatenate(target_bkg), axis=0)
    ref_sky = median(np.concatenate(ref_bkg), axis=0)
    sci = sci - sky
    refs = [
        median(frame, axis=0) - ref_sky for frame in exposures("jw090560" + ref_obs)
    ]
    assert len(refs) == 9 and len(sci) == 6
    cd = np.array([[h["CD1_1"], h["CD1_2"]], [h["CD2_1"], h["CD2_2"]]]) * 3600.0
    scale = np.sqrt(abs(np.linalg.det(cd)))
    center = np.array([h["CRPIX1"] - 1, h["CRPIX2"] - 1])
    yy, xx = np.indices(sci.shape[1:])
    radius = np.hypot(xx - center[0], yy - center[1]) * scale
    residual = rdi(sci, refs, radius)
    med = median(residual, axis=0)
    med[np.isfinite(residual).sum(axis=0) < 3] = np.nan
    primary_psf, inj_psf, calcent, caldata = calibrators(
        sorted(args.inputs.glob("jw015380" + calobs + "*fits"))
    )
    kernel = primary_psf[11:20, 11:20]
    trans = fits.getdata(
        next(args.inputs.glob("JWST_MIRI_" + args.band + "*fits")), 0
    ).astype(float)
    use = (ndi.minimum_filter(trans, size=5) > 0.85) & (trans >= TRANSMISSION_FLOOR)
    amp = matched(med, kernel)
    cut, snr, stats = levels(amp, radius, center, scale, use, beam_px)
    split = [
        levels(
            matched(median(part, axis=0), kernel), radius, center, scale, use, beam_px
        )[1]
        for part in [residual[::2], residual[1::2]]
    ]
    peaks = amp == ndi.maximum_filter(np.nan_to_num(amp, nan=-1e30), size=5)
    hits = []
    for y, x in zip(*np.where((amp > cut) & peaks)):
        east, north = cd @ ([x - center[0], y - center[1]])
        hits.append(
            dict(
                x0=x,
                y0=y,
                east_arcsec=east,
                north_arcsec=north,
                sep_arcsec=np.hypot(east, north),
                pa_deg=np.degrees(np.arctan2(east, north)) % 360,
                snr=snr[y, x],
                odd_snr=split[0][y, x],
                even_snr=split[1][y, x],
                flux_uJy=amp[y, x] * h["PIXAR_SR"] * 1e12 / trans[y, x],
            )
        )
    curves = []
    trials = []
    pix2ujy = h["PIXAR_SR"] * 1e12
    # PSF candidates are displaced across the quadrant phase plate. The
    # map gives relative intensity transmission at each trial coordinate.
    for rho, noise in zip(RADII, stats):
        nominal = noise["threshold"] * pix2ujy
        test = []
        for ang in (np.arange(48) + 0.37) * 2 * np.pi / 48:
            x, y = center + rho / scale * np.array([np.cos(ang), np.sin(ang)])
            ix, iy = round(x), round(y)
            cut_patch = cut[iy - 4 : iy + 5, ix - 4 : ix + 5]
            med_patch = med[iy - 8 : iy + 9, ix - 8 : ix + 9]
            shifted = ndi.shift(inj_psf, [y - iy, x - ix], order=3)
            # Sub-pixel interpolation is truncated to enforce disjoint support.
            shifted[np.hypot(*np.mgrid[-15:16, -15:16]) > 13] = 0.0
            assert rho - 13 * scale > 1.55
            transmission = ndi.map_coordinates(trans, [[y], [x]], order=1)[0]
            source = shifted[7:24, 7:24] * transmission
            blank = fit_at(med_patch, kernel)
            sensitivity = fit_at(
                np.where(np.isfinite(med_patch), source, np.nan), kernel
            )
            offsets = np.hypot(
                ix + np.arange(-4, 5)[None, :] - x, iy + np.arange(-4, 5)[:, None] - y
            )
            valid = np.isfinite(cut_patch) & (offsets <= 1.5)
            pre = np.any(
                (blank == ndi.maximum_filter(np.nan_to_num(blank, nan=-1e30), size=5))
                & valid
                & (blank > cut_patch)
            )
            eligible = bool(use[iy, ix] and np.isfinite(amp[iy, ix]))
            test.append(
                (np.degrees(ang), blank, sensitivity, cut_patch, valid, pre, eligible)
            )
        neligible = sum(q[-1] for q in test)
        frac = []
        grid = nominal * np.array(
            [0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5, 1.8, 2.2, 3.0, 4.0]
        )
        for uj in grid:
            success = 0
            for ang, base, unit, limit, valid, pre, eligible in test:
                field = base + uj / pix2ujy * unit
                imax = field == ndi.maximum_filter(
                    np.nan_to_num(field, nan=-1e30), size=5
                )
                recovered = bool(
                    eligible and not pre and np.any(imax & valid & (field > limit))
                )
                success += int(recovered)
                trials.append(
                    dict(
                        radius_arcsec=rho,
                        detector_angle_deg=round(ang, 4),
                        eligible=int(eligible),
                        flux_uJy=uj,
                        contrast=uj * 1e-6 / args.host_jy,
                        recovered=int(recovered),
                    )
                )
            frac.append(success / neligible)

        def lim(level):
            return next(
                (float(grid[i]) for i in range(len(grid)) if min(frac[i:]) >= level),
                None,
            )

        f50, f95 = lim(0.5), lim(0.95)
        curves.append(
            dict(
                radius_arcsec=rho,
                total_angles=48,
                eligible_angles=neligible,
                f50_uJy=f50,
                f95_uJy=f95,
                contrast50=None if f50 is None else f50 * 1e-6 / args.host_jy,
                contrast95=None if f95 is None else f95 * 1e-6 / args.host_jy,
            )
        )
    # One end-to-end injection tests local equivalence with the pre-fit pipeline.
    ix = int(center[0]) + 33
    iy = int(center[1]) + 26
    stim = sci.copy()
    sig = inj_psf * trans[iy, ix] * 30
    stim[:, iy - 15 : iy + 16, ix - 15 : ix + 16] += sig
    check = rdi(stim, refs, radius)
    check_im = median(check, axis=0)
    check_im[np.isfinite(check).sum(axis=0) < 3] = np.nan
    predicted = med.copy()
    predicted[iy - 15 : iy + 16, ix - 15 : ix + 16] += sig
    np.testing.assert_allclose(predicted, check_im, atol=1e-10, equal_nan=True)
    coverage = np.isfinite(snr[(radius >= INNER) & (radius < OUTER)]).mean()
    result = dict(
        band=args.band,
        target=h0["TARGNAME"],
        science_mjd=h0["EXPSTART"],
        date_obs=h0["DATE-OBS"],
        eff_time_s=h0["EFFEXPTM"],
        pipeline=h0["CAL_VER"],
        crds=h0["CRDS_CTX"],
        pixel_scale_arcsec=scale,
        center_xy0=center.tolist(),
        host_jy=args.host_jy,
        search_arcsec=[INNER, OUTER],
        annulus_coverage=coverage,
        n_peaks=len(hits),
        peak_score=float(np.nanmax(snr)),
        n_injection_trials=len(trials),
        calibrator_centers_xy=calcent,
        flux_calibration=caldata,
        negative_peaks=int(
            sum(
                np.sum(
                    (radius >= v["lo_arcsec"])
                    & (radius < v["hi_arcsec"])
                    & use
                    & (
                        (-amp)
                        == ndi.maximum_filter(np.nan_to_num(-amp, nan=-1e30), size=5)
                    )
                    & (amp < 2 * v["mean"] - v["threshold"])
                )
                for v in stats
            )
        ),
    )
    for filename, rows in [
        ("miri_noise.csv", stats),
        ("miri_contrast.csv", curves),
        ("miri_injections.csv", trials),
    ]:
        write_table(out / (prefix + filename), list(rows[0]), rows)
    write_table(
        out / (prefix + "miri_peaks.csv"),
        [
            "x0",
            "y0",
            "east_arcsec",
            "north_arcsec",
            "sep_arcsec",
            "pa_deg",
            "snr",
            "odd_snr",
            "even_snr",
            "flux_uJy",
        ],
        hits,
    )
    (out / (prefix + "miri_measurement.json")).write_text(
        json.dumps(result, indent=2) + "\n"
    )
    np.savez_compressed(
        out / (prefix + "miri_maps.npz"),
        residual_MJy_sr=med.astype("f4"),
        snr=snr.astype("f4"),
        threshold=cut.astype("f4"),
        amplitude=amp.astype("f4"),
        relative_transmission=trans.astype("f4"),
    )
    print(result, curves, hits)


if __name__ == "__main__":
    main()
