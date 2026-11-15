#!/usr/bin/env python3
"""Adversarial, reproducible fixed-position tests against PR #41 frozen FITS data.

Run from any directory: python results/wise-j0458-heldout/recheck.py
Requires numpy and astropy. The FITS reader and beam kernel are independent
of PR #41's Python code, while using the same frozen input images.
This writes ONLY results/wise-j0458-heldout/measurements.json.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "wise-j0458-radio"

RA, DEC, MJD, PMRA, PMDEC = 74.7268, 64.5825, 60859.0, 218.0, 300.0


def predicted(date):
    """Linear ICRS propagation, matching the catalog-only PR #41 hypothesis."""
    years = (Time(date).mjd - MJD) / 365.25
    return (RA + PMRA * years / 3.6e6 / np.cos(np.deg2rad(DEC)),
            DEC + PMDEC * years / 3.6e6)


class Image:
    """Minimal independent reader and native restoring-beam evaluator."""
    def __init__(self, path):
        with fits.open(path) as hdus:
            self.header = hdus[0].header.copy()
            self.image = np.asarray(hdus[0].data.squeeze(), dtype=float) * 1e3
        self.wcs = WCS(self.header).celestial
        self.date = self.header["DATE-OBS"]
        self.yy, self.xx = np.indices(self.image.shape)
        self.jac = self.wcs.pixel_scale_matrix * 3600
        self.pa = np.deg2rad(self.header["BPA"])
        self.axes = np.array([self.header["BMAJ"], self.header["BMIN"]]) * 3600

    def kernel(self, dx, dy):
        east, north = self.jac @ np.array([dx, dy])
        u = (east * np.sin(self.pa) + north * np.cos(self.pa)) / self.axes[0]
        v = (east * np.cos(self.pa) - north * np.sin(self.pa)) / self.axes[1]
        return np.exp(-4 * np.log(2) * (u*u + v*v))


LOFAR_CATALOG = PARENT / "independent-radio/lotss-source-cone.csv"
POINTING_CATALOG = PARENT / "independent-radio/pointings-within-3deg.csv"


def robust(xs):
    xs = np.asarray(xs, dtype=float)
    median = float(np.median(xs))
    return median, float(1.4826 * np.median(np.abs(xs - median)))


def patch_fit(im, anchors, source_skies, shift=(0, 0), half=9):
    """Beam-template least squares, preserving a common patch for the nulls.

    The columns contain one Gaussian beam for each fixed sky coordinate,
    a constant and local x/y slopes. Frequencies, beam sizes and observed
    dates are always native to the input visit. Do not interpret residual
    sums of squares as independent-pixel chi-squared.
    """
    anchor_pix = [np.asarray(im.wcs.world_to_pixel_values(*p)) for p in anchors]
    mid = (anchor_pix[0] + anchor_pix[1]) / 2 + np.asarray(shift)
    dx, dy = im.xx - mid[0], im.yy - mid[1]
    take = (np.abs(dx) < half) & (np.abs(dy) < half) & np.isfinite(im.image)
    x, y = im.xx[take], im.yy[take]
    cols = []
    for sky in source_skies:
        px, py = im.wcs.world_to_pixel_values(*sky)
        px += shift[0]
        py += shift[1]
        cols.append(im.kernel(x - px, y - py))
    cols.extend([np.ones_like(x), dx[take], dy[take]])
    design = np.column_stack(cols)
    values = im.image[take]
    beta = np.linalg.lstsq(design, values, rcond=None)[0]
    sse = float(np.sum((values - design @ beta) ** 2))
    return {"amplitudes": [float(a) for a in beta[:len(source_skies)]],
            "sse_mjy2": sse, "n_pixels": len(values)}


def main():
    images = [Image(PARENT / f"inputs/VLASS{i}.fits") for i in range(1, 5)]
    moving = [predicted(im.date) for im in images]
    with LOFAR_CATALOG.open() as stream:
        lotss = next(row for row in csv.DictReader(stream)
                     if row["source_name"] == "ILTJ045854.31+643455.8")
    static = (float(lotss["ra"]), float(lotss["dec"]))
    offsets = [(dx, dy) for dy in range(-30, 31, 6)
               for dx in range(-30, 31, 6)
               if 17 < np.hypot(dx, dy) < 31]
    output = {"method": "fixed sky positions, native-beam OLS in common 19x19 patches",
              "units": "mJy/beam; SSE is summed squared per-pixel residual in (mJy/beam)^2",
              "parent": "results/wise-j0458-radio/",
              "lofar_static_j2000_deg": static,
              "field_controls": len(offsets), "epochs": []}
    contrast_control = {}
    for i, (im, star) in enumerate(zip(images, moving), 1):
        anchor = (star, static)
        moving_only = patch_fit(im, anchor, [star])
        static_only = patch_fit(im, anchor, [static])
        two_components = patch_fit(im, anchor, [star, static])
        control_both = [patch_fit(im, anchor, [star, static], shift=o)
                        for o in offsets]
        control_moving = [patch_fit(im, anchor, [star], shift=o) for o in offsets]
        control_static = [patch_fit(im, anchor, [static], shift=o) for o in offsets]
        ctrl_a = [x["amplitudes"][0] for x in control_both]
        ctrl_b = [x["amplitudes"][1] for x in control_both]
        ctrl_diff = np.array([s["sse_mjy2"] - m["sse_mjy2"]
                              for m, s in zip(control_moving, control_static)])
        diff = static_only["sse_mjy2"] - moving_only["sse_mjy2"]
        m0, sm = robust(ctrl_a)
        s0, ss = robust(ctrl_b)
        dm, ds = robust(ctrl_diff)
        contrast_control[i] = (diff, ctrl_diff, ds)
        sep = 3600 * math.hypot((star[0] - static[0]) * math.cos(math.radians(star[1])),
                               star[1] - static[1])
        row = {"epoch": i, "date": im.date, "star_static_sep_arcsec": sep,
               "moving": moving_only, "static_lotss": static_only,
               "moving_plus_static": two_components,
               "control_dual_moving": {"median": m0, "mad_scale": sm,
                                       "robust_z": (two_components["amplitudes"][0] - m0) / sm},
               "control_dual_static": {"median": s0, "mad_scale": ss,
                                       "robust_z": (two_components["amplitudes"][1] - s0) / ss},
               "sse_static_minus_moving": diff,
               "contrast_control": {"median": dm, "mad_scale": ds,
                                    "rank_ge_including_target": 1 + int(np.sum(ctrl_diff >= diff))},
               "window_sensitivity": {
                   str(2 * half + 1): patch_fit(im, anchor, [star, static], half=half)
                   for half in (7, 9, 11)}}
        # Independently lock two *stationary* sources to the endpoint IR
        # coordinates. Neither source is allowed to move, but each flux can
        # change independently each epoch.
        ends = (moving[0], moving[-1])
        row["two_stationary_endpoints"] = {
            "both": patch_fit(im, ends, list(ends)),
            "single_actual_moving": patch_fit(im, ends, [star])}
        output["epochs"].append(row)
    # Off-source pseudo-tracks use the same displacement in both endpoint maps.
    comparison = np.zeros(len(offsets))
    selected = 0
    for i in (1, 4):
        d, ctrl, scale = contrast_control[i]
        comparison += ctrl / scale
        selected += d / scale
    output["endpoint_joint"] = {
        "test": "2017+2025 static-LoTSS versus moving SSE contrast / control MAD",
        "target": float(selected),
        "control_median": robust(comparison)[0],
        "control_mad": robust(comparison)[1],
        "n_controls": len(offsets),
        "rank_ge_including_target": 1 + int(np.sum(comparison >= selected)),
        "interpretation": "local descriptive rank, not survey-wide false-alarm probability"}
    with POINTING_CATALOG.open() as stream:
        points = list(csv.DictReader(stream))
    dt_years = (max(float(p["datelastobs"]) for p in points) -
                min(float(p["datefirstobs"]) for p in points)) / 365.25
    pm = math.hypot(218, 300) / 1000  # arcsec/year, from original source receipt
    parallax = 108.3  # mas, original source receipt
    d_pc = 1000 / parallax
    major = float(lotss["dc_maj"])
    e_major = float(lotss["e_dc_maj"])
    angle = float(lotss["dc_pa"])
    pa_pm = math.degrees(math.atan2(218, 300))
    output["morphology"] = {
        "lofar_deconvolved_major_arcsec": major,
        "catalog_error_major_arcsec": e_major,
        "lofar_deconvolved_minor_arcsec": float(lotss["dc_min"]),
        "lofar_pa_deg_e_of_n": angle, "pm_pa_deg_e_of_n": pa_pm,
        "pa_abs_difference_deg": abs(((angle - pa_pm + 90) % 180) - 90),
        "pointing_earliest_mjd": min(float(p["datefirstobs"]) for p in points),
        "pointing_latest_mjd": max(float(p["datelastobs"]) for p in points),
        "local_pointings_n": len(points),
        "local_pointing_span_year": dt_years,
        "maximum_catalog_pm_smear_arcsec": pm * dt_years,
        "major_over_smear": major / (pm * dt_years),
        "major_projected_au_if_associated": major * d_pc,
        "major_over_beam": major / 6,
        "caveat": "unknown mosaic pointing weights and non-astrophysical smearing"}
    flux = np.array([.7361465797122816, -.15003158291554516,
                     .4446350765499496, .48658454374487514])
    sigma = np.array([.1143940974605618, .1479645936445649,
                      .11262726393400843, .15460080907016788])
    subsets = [(0, 1, 2, 3), (1, 2, 3), (0, 1, 2), (0, 2, 3)]
    outvar = {}
    for s in subsets:
        w = sigma[list(s)] ** -2
        mean = float(np.sum(flux[list(s)] * w) / np.sum(w))
        chi = float(np.sum(((flux[list(s)] - mean) / sigma[list(s)]) ** 2))
        outvar["-".join(str(i + 1) for i in s)] = {
            "weighted_mean_mjy": mean, "chi2": chi, "nominal_df": len(s) - 1}
    output["variability_constant_flux_diagnostic"] = outvar
    outpath = HERE / "measurements.json"
    outpath.write_text(json.dumps(output, indent=2) + "\n")
    print("Wrote", outpath, "with", len(offsets), "co-translated controls")
    for e in output["epochs"]:
        print(e["epoch"], "moving and static mJy",
              [round(x, 3) for x in e["moving_plus_static"]["amplitudes"]])
    print("joint local rank", output["endpoint_joint"]["rank_ge_including_target"])


if __name__ == "__main__":
    main()
