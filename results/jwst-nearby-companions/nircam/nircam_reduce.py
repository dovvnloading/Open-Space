"""Classical reference-star differential imaging and injection/recovery for GO 5835 and 8826.

Fetch the inputs listed in the selected field manifest first. Input/output paths can be supplied for replay.
"""

import argparse
import csv
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi
from scipy.signal import correlate
from scipy.stats import norm, t

# The discovery search was restricted to the fully illuminated outer region.
INNER, OUTER = 2.5, 5.5
RADIAL_BIN = 0.30
BEAM_PX = 3.2
FALSE_ALARM = norm.sf(5)
CALIBRATOR = "jw01537028"
FIELDS = {
    "mcc124": ("jw05835", "jw05835059001_03103_00001_nrcalong_calints.fits", 8),
    "ross619": ("jw08826", "jw08826009001_05101_00001_nrcalong_calints.fits", 3),
    "l788_34": ("jw08826", "jw08826011001_05101_00001_nrcalong_calints.fits", 3),
    "v374peg": ("jw08826", "jw08826020001_05101_00001_nrcalong_calints.fits", 3),
    "g99_49": ("jw08826", "jw08826022001_03103_00001_nrcalong_calints.fits", 3),
}


def write_table(path, fields, rows):
    with path.open("w", newline="") as out:
        table = csv.DictWriter(out, fieldnames=fields)
        table.writeheader()
        table.writerows(rows)


def cube(path):
    with fits.open(path) as h:
        assert h[0].header["APERNAME"] == "NRCA5_MASK335R"
        assert h[0].header["FILTER"] == "F444W"
        assert h["SCI"].header["BUNIT"] == "MJy/sr"
        data = h["SCI"].data.astype(float)
        data[(h["DQ"].data & 3) != 0] = np.nan
        return data, h["SCI"].header.copy(), h[0].header.copy()


def median(data, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(data, axis=axis)


def fill(image):
    missing = ~np.isfinite(image)
    idx = ndi.distance_transform_edt(
        missing, return_distances=False, return_indices=True
    )
    return image[tuple(idx)]


def offaxis_psf(paths):
    # Normalize by the flux inside 30 px per dither before cropping the PSF.
    # 30 px corresponds to 1.88 arcsec. Lost wings are a calibration systematic.
    stamps = []
    growth = []
    y, x = np.indices((320, 320))
    for path in sorted(paths):
        image = fill(median(cube(path)[0], axis=0))
        iy, ix = np.unravel_index(image.argmax(), image.shape)
        core = image[iy - 1 : iy + 2, ix - 1 : ix + 2]
        cy, cx = ndi.center_of_mass(core - np.min(core))
        cy += iy - 1
        cx += ix - 1
        rad = np.hypot(x - cx, y - cy)
        background = median(image[(rad >= 35) & (rad < 43)])
        image -= background
        flux = image[rad < 30].sum()
        growth.append(float(image[rad < 25].sum() / flux))
        cen = 18
        stamp = ndi.map_coordinates(
            image, [y[:37, :37] + cy - cen, x[:37, :37] + cx - cen], order=3
        )
        stamps.append(stamp / flux)
    return median(stamps, axis=0), growth


def reduce_cube(science, references, radius, sky_bounds=(6.3, 7.5)):
    fit = (radius >= 0.45) & (radius <= 0.95)
    sky = (radius >= sky_bounds[0]) & (radius <= sky_bounds[1])
    frames = []
    for data in references:
        img = median(data, axis=0)
        img -= median(img[sky])
        frames.append(img / median(img[fit]))
    template = fill(median(frames, axis=0))
    # Small differential pointing and intensity offsets are fitted using only
    # the inner halo. The eight targets enter the reference with equal weight.
    gy, gx = np.gradient(template)
    model = np.stack([np.ones_like(template), template, gx, gy], axis=-1)
    residual = []
    for img in science:
        ok = fit & np.isfinite(img)
        for step in range(3):
            coeff = np.linalg.lstsq(model[ok], img[ok], rcond=None)[0]
            diff = img - model @ coeff
            scale = 1.4826 * median(np.abs(diff[ok] - median(diff[ok])))
            ok = fit & np.isfinite(img) & (np.abs(diff - median(diff[ok])) <= 4 * scale)
        residual.append(diff)
    return np.array(residual)


def matched(image, kernel):
    # Missing pixels are excluded in a local constant + empirical PSF fit.
    mask = np.isfinite(image)
    ones = np.ones_like(kernel)
    z = np.where(mask, image, 0.0)

    def corr(a, b):
        return correlate(a, b, mode="same", method="fft")

    n = np.maximum(1, corr(mask.astype(float), ones))
    k = corr(mask.astype(float), kernel)
    denom = corr(mask.astype(float), kernel**2) - k**2 / n
    numerator = corr(z, kernel) - k * corr(z, ones) / n
    good = (n >= 0.83 * kernel.size) & (
        denom > 0.6 * np.sum((kernel - kernel.mean()) ** 2)
    )
    return np.divide(numerator, denom, out=np.full_like(n, np.nan), where=good)


def thresholds(amplitude, radius, center, scale):
    threshold = np.full_like(amplitude, np.nan)
    standardized = np.full_like(amplitude, np.nan)
    rows = []
    for lo in np.arange(INNER, OUTER - 1e-5, RADIAL_BIN):
        hi = lo + RADIAL_BIN
        r = (lo + hi) / (2 * scale)
        n = int(2 * np.pi * r / BEAM_PX)
        theta = np.arange(n) * 2 * np.pi / n
        xx = np.rint(center[0] + r * np.cos(theta)).astype(int)
        yy = np.rint(center[1] + r * np.sin(theta)).astype(int)
        vals = amplitude[yy, xx]
        vals = vals[np.isfinite(vals)]
        n = len(vals)
        mean, sd = np.mean(vals), np.std(vals, ddof=1)
        noise = sd * np.sqrt(1 + 1 / n)
        # One-sided per-resolution-element probability. No trial-factor claim.
        cutoff = mean + t.isf(FALSE_ALARM, n - 1) * noise
        where = (radius >= lo) & (radius < hi)
        standardized[where] = (amplitude[where] - mean) / noise
        threshold[where] = cutoff
        rows.append(
            dict(
                lo_arcsec=round(lo, 2),
                hi_arcsec=round(hi, 2),
                n_noise=n,
                mean=float(mean),
                std=float(sd),
                cutoff=float(cutoff),
            )
        )
    return threshold, standardized, rows


def fit_at(patch, kernel):
    # Matched amplitude on a 9-pixel square for the blind 5-pixel peak criterion.
    out = np.full((9, 9), np.nan)
    for i, j in np.ndindex(out.shape):
        z = patch[i : i + kernel.shape[0], j : j + kernel.shape[1]]
        good = np.isfinite(z)
        k = kernel[good]
        if good.mean() < 0.83:
            continue
        kc = k - k.mean()
        if np.sum(kc**2) > 0.6 * np.sum((kernel - kernel.mean()) ** 2):
            out[i, j] = np.dot(z[good], kc) / np.sum(kc**2)
    return out


def injection_trials(
    residual_cube,
    kernel,
    psf,
    threshold,
    radius,
    center,
    scale,
    host_jy,
    jy_per_unit,
    trans,
):
    med = median(residual_cube, axis=0)
    rows, summary = [], []
    for rho in (2.65, 2.95, 3.55, 4.15, 4.75, 5.35):
        ang = (np.arange(36) + 0.371) * 2 * np.pi / 36
        trial = []
        for angle in ang:
            x = center[0] + rho / scale * np.cos(angle)
            y = center[1] + rho / scale * np.sin(angle)
            ix, iy = int(round(x)), int(round(y))
            sl = np.s_[iy - 8 : iy + 9, ix - 8 : ix + 9]
            adj = ndi.shift(psf, (y - iy, x - ix), order=3)
            source = adj[10:27, 10:27] * trans[sl]
            patch = med[sl]
            # Real source photons undergo RDI before stacking. The cropped PSF
            # is disjoint from the inner fit pixels. Translation equivariance
            # of a nanmedian then makes this local computation exactly linear.
            assert rho - 18 * np.sqrt(2) * scale > 0.95
            blank = fit_at(patch, kernel)
            unit = fit_at(np.where(np.isfinite(patch), source, np.nan), kernel)
            cut = threshold[iy - 4 : iy + 5, ix - 4 : ix + 5]
            offsets = np.hypot(
                ix + np.arange(-4, 5)[None, :] - x, iy + np.arange(-4, 5)[:, None] - y
            )
            valid = np.isfinite(blank) & np.isfinite(cut) & (offsets <= 1.5)
            baseline_peaks = blank == ndi.maximum_filter(
                np.nan_to_num(blank, nan=-1e30), size=5
            )
            baseline = np.any(baseline_peaks & valid & (blank > cut))
            # Require a local maximum within 1.5 pixels, with finite throughput.
            trial.append((angle, blank, unit, cut, valid, offsets, baseline))
        typical = np.nanmedian(threshold[(abs(radius - rho) < 0.03)]) * jy_per_unit
        levels = typical * np.array(
            [0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5, 1.8, 2.2, 3.0, 4.0]
        )
        fractions = []
        for flux in levels:
            detections = 0
            for angle, blank, unit, cut, valid, offsets, baseline in trial:
                field = blank + flux / jy_per_unit * unit
                peaks = field == ndi.maximum_filter(
                    np.nan_to_num(field, nan=-1e30), size=5
                )
                detected = bool(np.any(peaks & valid & (field > cut)) and not baseline)
                detections += detected
                rows.append(
                    dict(
                        radius_arcsec=rho,
                        detector_angle_deg=round(np.degrees(angle), 3),
                        flux_uJy=flux * 1e6,
                        contrast=flux / host_jy,
                        recovered=int(detected),
                        preexisting_crossing=int(baseline),
                    )
                )
            fractions.append(detections / len(trial))

        def limit(level):
            ind = next(
                (j for j in range(len(levels)) if min(fractions[j:]) >= level), None
            )
            return None if ind is None else float(levels[ind])

        f50, f95 = limit(0.5), limit(0.95)
        summary.append(
            dict(
                radius_arcsec=rho,
                n_angles=len(trial),
                nominal_cut_uJy=typical * 1e6,
                f50_uJy=None if f50 is None else f50 * 1e6,
                f95_uJy=None if f95 is None else f95 * 1e6,
                contrast50=None if f50 is None else f50 / host_jy,
                contrast95=None if f95 is None else f95 / host_jy,
                max_recovered=round(max(fractions) * len(trial)),
            )
        )
    return rows, summary


def verify_injection(science, refs, psf, kernel, trans, center, radius, out):
    # Compare the fast calculation with an explicit, independent full injection
    # before reference fitting, including DQ, integration stack and detection.
    ix, iy = int(center[0]) + 57, int(center[1]) + 13
    stamp = ndi.shift(psf, (0.25, 0.38), order=3)
    sl = np.s_[iy - 18 : iy + 19, ix - 18 : ix + 19]
    injected = science.copy()
    injected[:, sl[0], sl[1]] += 25 * stamp * trans[sl]
    reprocessed = reduce_cube(injected, refs, radius)
    frame = median(reprocessed, axis=0)
    frame[np.isfinite(reprocessed).sum(axis=0) < 4] = np.nan
    full = matched(frame, kernel)[iy - 4 : iy + 5, ix - 4 : ix + 5]
    raw = out[iy - 8 : iy + 9, ix - 8 : ix + 9]
    signal = stamp[10:27, 10:27] * trans[iy - 8 : iy + 9, ix - 8 : ix + 9]
    predicted = fit_at(raw, kernel) + 25 * fit_at(
        np.where(np.isfinite(raw), signal, np.nan), kernel
    )
    np.testing.assert_allclose(full, predicted, rtol=2e-8, atol=1e-8, equal_nan=True)
    return float(np.nanmax(np.abs(full - predicted)))


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--inputs", type=Path, default=Path("files"))
    arg.add_argument("--out", type=Path, default=Path("."))
    arg.add_argument("--field", choices=FIELDS, default="mcc124")
    arg.add_argument("--host-jy", type=float, required=True)
    args = arg.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    prefix, science_name, nref = FIELDS[args.field]

    def output(suffix):
        return args.out / (args.field + "_" + suffix)

    science, header, primary = cube(args.inputs / science_name)
    refs = [
        p
        for p in sorted(args.inputs.glob(prefix + "*long_calints.fits"))
        if p.name != science_name
    ]
    standards = sorted(args.inputs.glob(CALIBRATOR + "*long_calints.fits"))
    assert len(refs) == nref and len(standards) == 5
    detection_psf, _ = offaxis_psf(standards[:2])
    psf, growth = offaxis_psf(standards[2:])
    kernel = detection_psf[14:23, 14:23]
    cd = np.array(
        [[header["CD1_1"], header["CD1_2"]], [header["CD2_1"], header["CD2_2"]]]
    )
    scale = 3600 * np.sqrt(abs(np.linalg.det(cd)))
    # Adopt the coronagraph mask center as the unresolved photocenter.
    center = np.array([header["CRPIX1"] - 1, header["CRPIX2"] - 1])
    yy, xx = np.indices(science.shape[1:])
    radius = np.hypot(xx - center[0], yy - center[1]) * scale
    ref_data = [cube(p)[0] for p in refs]
    res = reduce_cube(science, ref_data, radius)
    med = median(res, axis=0)
    med[np.isfinite(res).sum(axis=0) < 4] = np.nan
    res[:, ~np.isfinite(med)] = np.nan
    trans = fits.getdata(args.inputs / "jwst_nircam_psfmask_0019.fits", "SCI").astype(
        float
    )
    # The flux standard is itself off axis on the clear substrate.
    trans /= np.median(trans[170:190, 220:240])
    replay_error = verify_injection(
        science, ref_data, psf, kernel, trans, center, radius, med
    )
    amp = matched(med, kernel)
    amp[trans < 0.95] = np.nan
    cut, snr, rings = thresholds(amp, radius, center, scale)
    jy_per_unit = header["PIXAR_SR"] * 1e6
    split = []
    for part in (res[::2], res[1::2]):
        amap = matched(median(part, axis=0), kernel)
        _, smap, _ = thresholds(amap, radius, center, scale)
        split.append(smap)
    maxima = np.nan_to_num(amp, nan=-1e30) == ndi.maximum_filter(
        np.nan_to_num(amp, nan=-1e30), size=5
    )
    hits = []
    for iy, ix in zip(*np.where((amp > cut) & maxima)):
        east, north = cd @ (np.array([ix, iy]) - center) * 3600
        hits.append(
            dict(
                x0=int(ix),
                y0=int(iy),
                rho_arcsec=float(np.hypot(east, north)),
                pa_deg=float(np.degrees(np.arctan2(east, north)) % 360),
                east_arcsec=float(east),
                north_arcsec=float(north),
                snr=float(snr[iy, ix]),
                odd_snr=float(split[0][iy, ix]),
                even_snr=float(split[1][iy, ix]),
                flux_uJy=float(amp[iy, ix] * jy_per_unit * 1e6 / trans[iy, ix]),
            )
        )
    trial, curves = injection_trials(
        res, kernel, psf, cut, radius, center, scale, args.host_jy, jy_per_unit, trans
    )
    write_table(output("injections.csv"), list(trial[0]), trial)
    write_table(output("contrast.csv"), list(curves[0]), curves)
    write_table(output("noise.csv"), list(rings[0]), rings)
    write_table(
        output("peaks.csv"),
        [
            "x0",
            "y0",
            "rho_arcsec",
            "pa_deg",
            "east_arcsec",
            "north_arcsec",
            "snr",
            "odd_snr",
            "even_snr",
            "flux_uJy",
        ],
        hits,
    )
    np.savez_compressed(
        output("rdi_maps.npz"),
        residual_MJy_sr=med.astype("f4"),
        amplitude=amp.astype("f4"),
        threshold=cut.astype("f4"),
        snr=snr.astype("f4"),
    )
    negative_maxima = (-amp) == ndi.maximum_filter(
        np.nan_to_num(-amp, nan=-1e30), size=5
    )
    negative_count = sum(
        int(
            np.sum(
                negative_maxima
                & (radius >= row["lo_arcsec"])
                & (radius < row["hi_arcsec"])
                & (amp < 2 * row["mean"] - row["cutoff"])
            )
        )
        for row in rings
    )
    result = dict(
        target=primary["TARGNAME"],
        host_jy=args.host_jy,
        host="unresolved combined light",
        n_science=len(science),
        n_references=len(refs),
        n_calibrators=len(standards),
        scale_arcsec=scale,
        explicit_injection_max_abs_difference=replay_error,
        peak_snr=float(np.nanmax(snr)),
        negative_peaks_past_cut=negative_count,
        fit_arcsec=[0.45, 0.95],
        search_arcsec=[INNER, OUTER],
        n_threshold_peaks=len(hits),
        coverage=float(np.isfinite(snr[(radius >= INNER) & (radius < OUTER)]).mean()),
        min_transmission_relative=float(
            np.min(trans[(radius >= INNER) & (radius < OUTER)])
        ),
        calibration_growth_25px_over_30px=growth,
        n_injections=len(trial),
        center_xy0=center.tolist(),
        cal_ver=primary.get("CAL_VER"),
        crds_ctx=primary.get("CRDS_CTX"),
        expstart_mjd=primary["EXPSTART"],
        effective_exposure_s=primary["EFFEXPTM"],
    )
    (output("measurement.json")).write_text(json.dumps(result, indent=2) + "\n")
    print(result)
    print(hits)
    print(curves)


if __name__ == "__main__":
    main()
