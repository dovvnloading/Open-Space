"""F1140C empirical injections at the separation of HD 202628's visit 007 peak."""

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from miri_rolls import read

LEVELS = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0])
OFFSETS = np.array([1.1779125369960428, 3.0645695071040415])
RADIUS = np.linalg.norm(OFFSETS)
PA = np.rad2deg(np.arctan2(*OFFSETS)) % 360


def empirical_psfs(cache, trans):
    from standard_photometry import band_jy
    from nircam_reduce import median, fill
    from search import refine
    from vet import aperture

    flux = band_jy(cache / "hd2811_mod_006.fits", cache / "F1140C.vot")
    stamps = {}
    growth = []
    for f in sorted(cache.glob("jw01536041001*mirimage_cal.fits")):
        with fits.open(f) as hd:
            im = hd["SCI"].data.astype(float)
            im[(hd["DQ"].data & 3) != 0] = np.nan
            p = hd[0].header
            h = hd["SCI"].header
            assert p["FILTER"] == "F1140C" and p["TARGNAME"] == "HD 2811"
            assert p["EXP_TYPE"] == "MIR_CORONCAL" and h["BUNIT"] == "MJy/sr"
            assert im.shape == trans.shape
            cy, cx = np.unravel_index(np.nanargmax(im), im.shape)
            x, y = refine(im, cx, cy)
            # 70% aperture / sky annulus and aperture correction from STScI
            # Table 2. A local sky is also removed from the PSF stamp.
            # Some far-southern standards have a clipped 27-pixel annulus.
            if min(x, y, im.shape[1] - 1 - x, im.shape[0] - 1 - y) < 31:
                continue
            ratio = (
                aperture(
                    fill(im),
                    x,
                    y,
                    dict(radius=10.68, skyin=19.10, skyout=27.06, apcorr=1.461),
                )
                * h["PIXAR_SR"]
                * 1e6
                / flux
            )
            if not np.isfinite(ratio):
                continue
            yy, xx = np.indices(im.shape)
            r = np.hypot(xx - x, yy - y)
            back = median(im[(r >= 19.1) & (r < 27.06)])
            dy, dx = np.mgrid[-13:14, -13:14]
            # Linear recentering, independently for each dither.
            stamp = ndimage.map_coordinates(
                im - back, [dy + y, dx + x], order=1, prefilter=False
            )
            stamp[np.hypot(dx, dy) > 13] = 0
            # DQ defects outside the core can be filled by other dithers.
            reftrans = ndimage.map_coordinates(trans, [[y], [x]], order=1)[0]
            quad = (x >= 119.75, y >= 112.24)
            stamps.setdefault(quad, []).append(stamp / (flux * 1e6 * reftrans))
            growth.append(
                dict(
                    filename=f.name,
                    x=x,
                    y=y,
                    model_jy=flux,
                    aperture_corrected_over_model=ratio,
                    relative_transmission=reftrans,
                    quadrant=str(quad),
                )
            )
    if set(stamps) != {(False, False), (False, True), (True, False), (True, True)}:
        raise ValueError("Empirical standard missing a quadrant")
    psf = {q: median(st, axis=0) for q, st in stamps.items()}
    for q, img in psf.items():
        if not np.isfinite(img).all():
            raise ValueError(f"Missing standard pixels {q}")
    return psf, growth


def search_map(stack, h):
    from nircam_reduce import median
    from search import smooth_response

    im = median(stack, axis=0)
    valid = np.isfinite(im)
    y, x = np.indices(im.shape)
    cx, cy = h["CRPIX1"] - 1, h["CRPIX2"] - 1
    scale = np.sqrt(h["PIXAR_SR"]) * 206264.806
    radius = np.hypot(x - cx, y - cy) * scale
    good = (
        ndimage.minimum_filter(valid, size=5, mode="constant", cval=0)
        & (abs(x - cx) * scale > 0.37)
        & (abs(y - cy) * scale > 0.37)
        & (radius >= 2.5)
        & (radius < 5.5)
    )
    filt = smooth_response(np.where(valid, im, median(im)))
    noise = np.ones_like(im)
    score = np.zeros_like(im)
    for lo in np.arange(2.5, 5.5, 0.5):
        ring = good & (radius >= lo) & (radius < lo + 0.5)
        m = median(filt[ring])
        sig = 1.4826 * median(abs(filt[ring] - m))
        noise[ring] = sig
        score[ring] = (filt[ring] - m) / sig
    return dict(noise=noise, score=score, good=good), im, radius


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    sys.path[:0] = [str(a.merged / "nircam"), str(a.merged / "miri")]
    from nircam_reduce import reduce_cube, write_table, median
    from search import injected_positions

    trans = fits.getdata(a.cache / "F1140C_transmission.fits", 0).astype(float)
    psf, cal = empirical_psfs(a.cache, trans)
    write_table(a.out / "F1140C_standard_checks.csv", list(cal[0]), cal)
    rows = []
    verify = []
    for path in sorted(a.cache.glob("jw02538*mirimage_c1017_psfsub.fits")):
        pipeline, h, primary = read(path)
        w = WCS(h, naxis=2)
        cx, cy = h["CRPIX1"] - 1, h["CRPIX2"] - 1
        center = w.pixel_to_world(cx, cy)
        scale = np.sqrt(h["PIXAR_SR"]) * 206264.806
        native, _, _ = read(a.cache / path.name.replace("_c1017_psfsub", "_calints"))
        yy, xx = np.indices(native.shape[1:])
        rad = np.hypot(xx - cx, yy - cy) * scale
        for prefix in ["pipeline", "jw02538009", "jw02538014"]:
            refs = []
            method = "pipeline_post_KLIP"
            if prefix == "pipeline":
                cube = pipeline
            else:
                files = sorted(a.cache.glob(prefix + "*mirimage_calints.fits"))
                method = "native_RDI_" + read(files[0])[2]["TARGNAME"].replace(" ", "")
                refs = [read(f)[0] for f in files]
                cube = reduce_cube(native, refs, rad)
            obs, im, r = search_map(cube, h)
            for n, angle in enumerate(np.r_[np.arange(72) * 5.0 + 0.83, PA]):
                theta = np.deg2rad(angle)
                en = RADIUS * np.array([np.sin(theta), np.cos(theta)])
                pos = np.asarray(
                    w.world_to_pixel(center.spherical_offsets_by(*(en * u.arcsec)))
                ).astype(float)
                x, y = pos
                ix, iy = np.rint(pos).astype(int)
                q = (x >= cx, y >= cy)
                t = ndimage.map_coordinates(trans, [[y], [x]], order=1)[0]
                stamp = psf[q] * t
                # The empirical template comes from far outside the mask. Restrict
                # to >=85% modeled transmission to avoid strong 4QPM shape changes.
                eligible = bool(
                    obs["good"][iy, ix]
                    and t >= 0.85
                    and np.isfinite(im[iy - 3 : iy + 4, ix - 3 : ix + 4]).all()
                )
                near = []
                if eligible:
                    shifted = ndimage.shift(
                        np.pad(stamp, 1), (y - iy, x - ix), order=1, prefilter=False
                    )
                    obs["psf"] = np.where(
                        np.isfinite(im[iy - 14 : iy + 15, ix - 14 : ix + 15]),
                        shifted,
                        0.0,
                    )
                    # Reuse the parent's before-filter injected-position search.
                    # Only use pixels fully sampled by the injected compact core.
                    positions = injected_positions(obs, np.array([ix, iy]), LEVELS)
                    near = [
                        bool(np.any(np.linalg.norm(pp - pos, axis=1) * scale <= 0.20))
                        for pp in positions
                    ]
                    if refs and n == 72:
                        # Refit all native integrations and compare to local linear
                        # injection in a region disjoint from the inner halo fit.
                        canvas = np.zeros_like(im)
                        canvas[iy - 13 : iy + 14, ix - 13 : ix + 14] = stamp
                        canvas = ndimage.shift(
                            canvas, (y - iy, x - ix), order=1, prefilter=False
                        )
                        assert not np.any(canvas[(r >= 0.45) & (r <= 0.95)])
                        re = median(reduce_cube(native + 5 * canvas, refs, r), axis=0)
                        err = np.nanmax(np.abs(re - (im + 5 * canvas)))
                        assert err < 1e-9
                        verify.append(
                            dict(
                                roll=primary["OBSERVTN"],
                                method=method,
                                flux_uJy=5.0,
                                max_difference_MJy_sr=err,
                            )
                        )
                for k, level in enumerate(LEVELS[1:], 1):
                    rows.append(
                        dict(
                            roll=primary["OBSERVTN"],
                            method=method,
                            trial=n,
                            observed_PA=int(n == 72),
                            radius_arcsec=RADIUS,
                            pa_deg=angle,
                            x=x,
                            y=y,
                            transmission=t,
                            eligible=int(eligible),
                            flux_uJy=level,
                            baseline_peak=int(near[0]) if eligible else "",
                            recovered=int(near[k] and not near[0]) if eligible else "",
                        )
                    )
    write_table(a.out / "F1140C_recovery_trials.csv", list(rows[0]), rows)
    write_table(a.out / "F1140C_native_replay.csv", list(verify[0]), verify)
    summary = []
    for method in sorted({r["method"] for r in rows}):
        for level in LEVELS[1:]:
            entries = [
                r
                for r in rows
                if r["method"] == method
                and r["flux_uJy"] == level
                and not r["observed_PA"]
            ]
            r1 = [r for r in entries if int(r["roll"]) == 7]
            r2 = [r for r in entries if int(r["roll"]) == 8]
            both = [(u, v) for u, v in zip(r1, r2) if u["eligible"] and v["eligible"]]
            for roll, label in [(r1, "007"), (r2, "008")]:
                valid = [r for r in roll if r["eligible"]]
                count = sum(r["recovered"] for r in valid)
                summary.append(
                    dict(
                        method=method,
                        roll=label,
                        flux_uJy=level,
                        radius_arcsec=RADIUS,
                        eligible=len(valid),
                        geometric=72,
                        recovered=count,
                        fraction=count / len(valid),
                    )
                )
            count = sum(u["recovered"] and v["recovered"] for u, v in both)
            summary.append(
                dict(
                    method=method,
                    roll="both",
                    flux_uJy=level,
                    radius_arcsec=RADIUS,
                    eligible=len(both),
                    geometric=72,
                    recovered=count,
                    fraction=count / len(both),
                )
            )
    write_table(a.out / "F1140C_recovery.csv", list(summary[0]), summary)


if __name__ == "__main__":
    main()
