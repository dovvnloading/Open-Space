"""Measure the polar-source images and compare SED fits with and without LoDeSS.

Example: PYTHONPATH=../release:.. python j100005.py ../release/data raw evidence
"""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path
from urllib.parse import urlencode

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
from scipy.optimize import brentq, minimize
from acquire_lodess import digest, retrieve
from measure_vlssr import PARAMETERS, SERVICE
from model_spectra import fit, spectrum
from photometry import RadioImage, mad
from resolution_check import aperture

RA, DEC = 150.024629507167, 81.4507636411097
TARGET = SkyCoord(RA * u.deg, DEC * u.deg)
POSITEXT = "10:00:05.91 +81:27:02.7"


def targeted_archives(raw):
    raw.mkdir(exist_ok=True, parents=True)
    stamps = {}
    receipts = []
    for beam, suffix in [(6, ""), (20, "-low")]:
        url = (
            "https://lofar-surveys.org/dr3"
            + suffix
            + "-cutout.fits?"
            + urlencode({"pos": POSITEXT, "size": 12})
        )
        path = raw / f"lotss{beam}_1000.fits"
        if not path.exists():
            retrieve(url, path)
        stamps[beam] = path
        receipts.append(dict(file=path.name, sha256=digest(path), request=url))
    path = raw / "vlss_1000.fits"
    params = dict(PARAMETERS, RA="10 00 05.91", Dec="+81 27 02.7")
    if not path.exists():
        cmd = ["curl", "-fsSL", "--retry", "2", "--max-time", "120"]
        for k, v in params.items():
            cmd += ["--data-urlencode", f"{k}={v}"]
        subprocess.run([*cmd, SERVICE, "-o", str(path)], check=True)
    stamps[75] = path
    receipts.append(
        dict(file=path.name, sha256=digest(path), request=SERVICE, post_form=params)
    )
    return stamps, receipts


def low_photometry(stamps, data):
    low = []
    for field, quality in [("P139+79", "good"), ("P210+80", "mediocre")]:
        with RadioImage(data / "images" / f"{field}-I.fits") as im:
            hdr = im.hdus[0].header
            offset = TARGET.separation(
                SkyCoord(hdr["CRVAL1"] * u.deg, hdr["CRVAL2"] * u.deg)
            ).deg
            low.append(
                dict(
                    field=field,
                    quality=quality,
                    date_obs=hdr["DATE-OBS"],
                    offaxis_deg=offset,
                    **im.measure(RA, DEC),
                )
            )
    with tempfile.TemporaryDirectory() as tmp:
        with fits.open(stamps[75]) as h:
            hdr = h[0].header
            if abs(hdr["CRVAL3"] / 1e6 - 73.8) > 0.1:
                raise ValueError("Unexpected VLSSr frequency")
            # Official postage stamps omit the uniform 75-arcsec beam in the header.
            hdr["BMAJ"] = hdr["BMIN"] = 75 / 3600
            f = Path(tmp) / "beam.fits"
            h.writeto(f)
        with RadioImage(f) as im:
            result = im.measure(RA, DEC)
            x, y = map(float, im.wcs.world_to_pixel_values(RA, DEC))
            d = im.data
            yy, xx = np.ogrid[: d.shape[0], : d.shape[1]]

            def integrated(xc, yc, rad):
                r = im.pixel_arcsec * np.hypot(xx - xc, yy - yc)
                bg = np.nanmedian(d[(r >= 180) & (r <= 270)])
                return np.nansum(d[r <= rad] - bg) / (
                    im.beam_area * (1 - np.exp(-4 * np.log(2) * (rad / 75) ** 2))
                )

            apertures = []
            for radius in [60, 75, 90, 105, 112.5, 120, 150]:
                grid = [
                    integrated(x + i * 3 * im.fwhm, y + j * 3 * im.fwhm, radius)
                    for i in range(-4, 5)
                    for j in range(-4, 5)
                    if 2 <= np.hypot(i, j) <= 4
                ]
                apertures.append(
                    dict(
                        radius_arcsec=radius,
                        raw_jy=integrated(x, y, radius),
                        stat_jy=mad(np.array(grid)),
                        controls=len(grid),
                        clean_bias_jy=0.66 * result["pixel_noise_jy_beam"],
                    )
                )
    return pd.DataFrame(low), pd.DataFrame(apertures), result


def lotss_images(stamps):
    rows = []
    for beam in (6, 20):
        with fits.open(stamps[beam]) as h:
            hdr = h[0].header
            native = hdr["BMAJ"] * 3600
            scale = abs(hdr["CDELT1"]) * 3600
            if abs(native - beam) > 0.01:
                raise ValueError("Unexpected LoTSS beam")
            x, y = WCS(hdr).celestial.world_to_pixel(TARGET)
            data = np.squeeze(h[0].data).astype(float)
            for matched in [beam, 45, 75]:
                smoothed = data
                if matched > beam:
                    width = np.sqrt(matched**2 - beam**2) / (
                        np.sqrt(8 * np.log(2)) * scale
                    )
                    smoothed = (
                        gaussian_filter(data, width, mode="nearest")
                        * (matched / beam) ** 2
                    )
                for radius in [45, 60, 90]:
                    flux, bg, noise = aperture(smoothed, x, y, scale, matched, radius)
                    rows.append(
                        dict(
                            native_arcsec=beam,
                            beam_arcsec=matched,
                            radius_arcsec=radius,
                            flux_jy=flux,
                            background_jy_beam=bg,
                            local_rms_jy_beam=noise,
                        )
                    )
    return pd.DataFrame(rows)


def peak_interval(data, best, statistic):
    """Profile the log-parabolic peak with log-amplitude and curvature refitted."""
    x = np.log(data.freq_mhz.to_numpy() / 144)
    flux = data.flux_jy.to_numpy()
    err = data.stat_jy.to_numpy()
    gain = data.frac_scale.to_numpy()
    peak_x = -best[1] / (2 * best[2])
    peak_log_flux = best[0] + best[1] * peak_x + best[2] * peak_x**2

    def boundary(xp):
        def objective(pars):
            mu = np.exp(np.clip(pars[0] + pars[1] * (x - xp) ** 2, -100, 30))
            var = err**2 + (gain * mu) ** 2
            return np.sum((flux - mu) ** 2 / var + np.log(var))

        r = minimize(
            objective,
            [peak_log_flux, best[2]],
            method="Nelder-Mead",
            bounds=[(-8, 8), (-3, -1e-9)],
            options={"maxiter": 700},
        )
        if not r.success:
            raise RuntimeError("Peak profile did not converge")
        return r.fun - statistic - 3.841459

    ends = []
    for direction in [-1, 1]:
        for step in np.geomspace(0.02, 5, 32):
            test = peak_x + direction * step
            if boundary(test) > 0:
                root = brentq(boundary, *sorted((peak_x, test)))
                ends.append(144 * np.exp(root))
                break
        else:
            ends.append(np.nan)
    return ends


def compare_fits(rev, low, vlss, out):
    primary = low[low.quality == "good"].iloc[0]
    v = vlss[vlss.radius_arcsec == 90].iloc[0]
    points = [
        ("LoDeSS", 23, primary.psf_jy, primary.psf_noise_jy, 0.39),
        ("VLSSr", 73.7939, v.raw_jy + v.clean_bias_jy, v.stat_jy, 0.15),
        ("LoTSS", 144, rev.lotss_flux_jy, rev.lotss_err_jy, 0.07),
        # This catalogue's quoted integrated error already includes 10% calibration.
        ("TGSS", 150, rev.gmrtas150m_flux_jy, rev.gmrtas150m_err_jy, 0),
        ("WENSS", 352, rev.wenss_flux_jy, rev.wenss_err_jy, 0.08),
        ("NVSS", 1400, rev.nvss_flux_jy, rev.nvss_err_jy, 0.05),
    ]
    points = pd.DataFrame(
        points, columns=["survey", "freq_mhz", "flux_jy", "stat_jy", "frac_scale"]
    )
    points.to_csv(out / "j100005_photometry.csv", index=False)
    fits_data = []
    for choice in [
        "baseline",
        "no_23",
        "brighter_23",
        "bright_vlssr",
        "no_lofar",
        "extra_20pct_epoch",
        "miyun",
        "miyun_no_23",
        "miyun_20pct",
        "catalogue_325",
        "350plus",
    ]:
        p = points.copy()
        if choice.startswith("miyun"):
            # Miyun 1988-02-27 field 153, unresolved catalogue peak; no reported integral.
            # Its 232 MHz gain is treated as independent of the remaining instruments.
            gain = 0.20 if choice == "miyun_20pct" else 0.166
            p.loc[len(p)] = ["Miyun", 232, 2.13, 2.13 / 41.8, gain]
        if choice == "catalogue_325":
            p.loc[p.survey == "WENSS", "freq_mhz"] = 325.0
        if choice in ["no_23", "miyun_no_23"]:
            p = p[p.survey != "LoDeSS"]
        if choice == "no_lofar":
            p = p[~p.survey.isin(["LoTSS", "LoDeSS"])]
        if choice == "350plus":
            p = p[p.freq_mhz >= 350]
        if choice == "brighter_23":
            rep = low[low.quality == "mediocre"].iloc[0]
            p.loc[p.survey == "LoDeSS", ["flux_jy", "stat_jy"]] = [
                rep.psf_jy,
                rep.psf_noise_jy,
            ]
        if choice == "bright_vlssr":
            # Catalogue peak times Gaussian/beam solid-angle ratio, with the point-source CLEAN correction.
            p.loc[p.survey == "VLSSr", "flux_jy"] = (
                2.19054 * 88.43868 * 74.9988 / 75**2 + v.clean_bias_jy
            )
        if choice == "extra_20pct_epoch":
            p["frac_scale"] = np.hypot(p.frac_scale, 0.20)
        for model in ["power", "curved", "absorbed"] if len(p) > 3 else ["power"]:
            pars, stat, peak = fit(p, model)
            ci = peak_interval(p, pars, stat) if model == "curved" else [np.nan, np.nan]
            fits_data.append(
                dict(
                    sample=choice,
                    n=len(p),
                    model=model,
                    log_amplitude=pars[0],
                    alpha=pars[1],
                    shape=pars[2] if len(pars) == 3 else np.nan,
                    minus2logL_constant=stat,
                    AIC_constant=stat + 2 * len(pars),
                    peak_mhz=peak,
                    peak_low95=ci[0],
                    peak_high95=ci[1],
                    predict23_jy=spectrum(pars, np.log(23 / 144), model),
                    predict74_jy=spectrum(pars, np.log(73.7939 / 144), model),
                )
            )
    fits_data = pd.DataFrame(fits_data)
    fits_data.to_csv(out / "j100005_fits.csv", index=False)
    return points, fits_data


def plot(points, models, low, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    err = np.hypot(points.stat_jy, points.frac_scale * points.flux_jy)
    ax.errorbar(
        points.freq_mhz,
        points.flux_jy,
        yerr=err,
        fmt="o",
        ms=6,
        capsize=3,
        color="#125579",
        label="Independent survey measurements",
    )
    repeat = low[low.quality == "mediocre"].iloc[0]
    ax.errorbar(
        [23],
        [repeat.psf_jy],
        yerr=[np.hypot(repeat.psf_noise_jy, 0.39 * repeat.psf_jy)],
        fmt="s",
        mfc="none",
        ms=6,
        color="#B46D19",
        capsize=3,
        label="Second LoDeSS pointing",
    )
    x = np.geomspace(20, 1550, 200)
    for sample, model, color, style, label in [
        ("baseline", "curved", "#15644e", "-", "Curved fit"),
        ("no_23", "curved", "#695192", ":", "Curved fit with 23 MHz held out"),
        ("350plus", "power", "#777", "--", "352–1400 MHz extrapolation"),
    ]:
        m = models[(models["sample"] == sample) & (models.model == model)].iloc[0]
        ax.plot(
            x,
            spectrum([m.log_amplitude, m.alpha, m["shape"]], np.log(x / 144), model),
            ls=style,
            color=color,
            label=label,
        )
    ax.errorbar(
        [232],
        [2.13],
        yerr=[np.hypot(2.13 / 41.8, 0.166 * 2.13)],
        fmt="D",
        ms=5,
        color="#ad4f44",
        capsize=3,
        label="Miyun (sensitivity fit)",
    )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Observed frequency (MHz)",
        ylabel="Flux density (Jy)",
        title="NVSS J100005+812702",
    )
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.18, which="both")
    fig.savefig(out / "j100005_sed.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ("data", "raw", "out"):
        ap.add_argument(arg, type=Path)
    args = ap.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    rev = pd.read_csv(args.out / "reverse_counterparts.csv.gz")
    rec = rev[(rev.name == "VLSSr J100007.6+812706") & (rev.field == "P139+79")].iloc[0]
    stamps, receipts = targeted_archives(args.raw)
    (args.out / "j100005_stamps.json").write_text(json.dumps(receipts, indent=2) + "\n")
    low, vlss, psf = low_photometry(stamps, args.data)
    (args.out / "j100005_vlssr_psf.json").write_text(json.dumps(psf, indent=2) + "\n")
    low.to_csv(args.out / "j100005_pointings.csv", index=False)
    vlss.to_csv(args.out / "j100005_vlssr_apertures.csv", index=False)
    lotss_images(stamps).to_csv(args.out / "j100005_lotss_apertures.csv", index=False)
    p, m = compare_fits(rec, low, vlss, args.out)
    plot(p, m, low, args.out)
