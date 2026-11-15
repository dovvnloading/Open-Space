"""Can an airmass- or seeing-dependent photometric error imitate the parallax signal?

Adds a difference-imaging-style systematic, g * (A - 1) * (X - X0), to the
no-parallax model (X is airmass or seeing), fitted linearly with the source and
blend fluxes, and compares the chi2 gain with that from parallax. Also bins the
residuals of magnified epochs by airmass and seeing.

Usage: python -I check/airmass.py <phot.dat> <fits.json> <stress.json> <out.json>
"""
import json
import os
import sys

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pspl as lens  # noqa: E402

LCO = EarthLocation(lat=-29.0146 * u.deg, lon=-70.7017 * u.deg, height=2380 * u.m)   # OGLE Warsaw telescope
NM = {"xatol": 1e-6, "fatol": 1e-5, "maxiter": 8000, "maxfev": 16000}


def solve(a, f, w, extra=None):
    cols = [a, np.ones_like(a)] + ([extra] if extra is not None else [])
    m = np.vstack(cols).T
    coef = np.linalg.lstsq(m * np.sqrt(w)[:, None], f * np.sqrt(w), rcond=None)[0]
    r = f - m @ coef
    return float((w * r * r).sum()), coef


def main():
    phot, fitjson, sjson, outp = sys.argv[1:5]
    t, mag, err, see, sky = np.loadtxt(phot).T
    prev = json.load(open(fitjson))
    k = json.load(open(sjson))["error_scale"]
    sols = {key: prev[key]["params"] for key in ("no_parallax", "u0+", "u0-")}
    f = 10 ** (-0.4 * (mag - lens.ZP))
    w = 1.0 / (f * np.log(10) / 2.5 * err * k) ** 2
    dn, de = lens.parallax_offsets(t)
    target = SkyCoord("17:52:38.09 -31:47:36.1", unit=(u.hourangle, u.deg))
    alt = target.transform_to(AltAz(obstime=Time(t, format="jd"), location=LCO)).alt.deg
    x = 1.0 / np.sin(np.radians(alt))
    s26 = np.where(t > 2461040.0)[0]
    g = int(np.argmax(np.diff(t[s26])))
    late = t > t[s26][g] + 0.5
    mag_ev = t > 2460676.0           # 2025-2026, while the source is magnified
    out = {"airmass_late_points": np.round(x[late], 2).tolist(),
           "airmass_2026_median_other": float(np.median(x[(t > 2461040.0) & ~late])),
           "airmass_quantiles_2025_2026": np.round(np.percentile(x[mag_ev], [5, 50, 95]), 2).tolist()}

    def fit(params0, parallax, extra_of):
        def chi(p):
            pp = list(p) if parallax else [p[0], p[1], p[2], 0.0, 0.0]
            if pp[2] <= 1:
                return 1e30
            a = lens.magnification(pp, t, dn, de)
            return solve(a, f, w, None if extra_of is None else extra_of(a))[0]
        x0 = list(params0) if parallax else list(params0[:3])
        r = minimize(chi, x0, method="Nelder-Mead", options=NM)
        r = minimize(chi, r.x, method="Nelder-Mead", options=NM)
        pp = list(r.x) if parallax else list(r.x) + [0.0, 0.0]
        a = lens.magnification(pp, t, dn, de)
        c, coef = solve(a, f, w, None if extra_of is None else extra_of(a))
        return pp, c, coef

    terms = {"none": None,
             "airmass": lambda a: (a - 1) * (x - 1.0),
             "seeing": lambda a: (a - 1) * (see - np.median(see)),
             "sky": lambda a: (a - 1) * (sky - np.median(sky)) / 100.0}
    res = {}
    for name, term in terms.items():
        p_np, c_np, coef_np = fit(sols["no_parallax"], False, term)
        p_pl, c_pl, coef_pl = fit(sols["u0+"], True, term)
        res[name] = {"chi2_no_parallax": c_np, "chi2_parallax_u0+": c_pl,
                     "no_parallax_tE": p_np[2], "parallax_piE": float(np.hypot(p_pl[3], p_pl[4])),
                     "parallax_tE": p_pl[2],
                     "term_over_fs_no_parallax": float(coef_np[2] / coef_np[0]) if term else 0.0,
                     "term_over_fs_parallax": float(coef_pl[2] / coef_pl[0]) if term else 0.0}
    base = res["none"]["chi2_parallax_u0+"]
    for name in res:
        res[name]["delta_vs_best_parallax"] = {
            "no_parallax": res[name]["chi2_no_parallax"] - base,
            "parallax": res[name]["chi2_parallax_u0+"] - base}
    out["systematic_terms"] = res

    # binned residuals of magnified epochs (late points excluded), in mmag
    def model_mag(p):
        a = lens.magnification(p, t, dn, de)
        _, coef = solve(a, f, w)
        return lens.ZP - 2.5 * np.log10(coef[0] * a + coef[1])
    r_np = mag - model_mag(sols["no_parallax"])
    r_pl = mag - model_mag(sols["u0+"])
    sel = mag_ev & ~late
    bins = {}
    for label, var, edges in (("airmass", x, [1.0, 1.1, 1.3, 1.6, 2.0, 4.0]),
                              ("seeing", see, [0, 4.5, 5.5, 6.5, 7.5, 20]),
                              ("sky", sky, [0, 300, 450, 600, 900, 1e5])):
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = sel & (var >= lo) & (var < hi)
            if m.sum() < 3:
                continue
            rows.append({"range": [lo, hi], "n": int(m.sum()),
                         "mean_resid_no_parallax_mmag": round(1000 * float(np.mean(r_np[m])), 1),
                         "mean_resid_parallax_mmag": round(1000 * float(np.mean(r_pl[m])), 1),
                         "sem_mmag": round(1000 * float(np.std(r_pl[m]) / np.sqrt(m.sum())), 1)})
        bins[label] = rows
    out["binned_residuals_2025_2026_without_late"] = bins
    out["late_points_resid_mmag"] = {"no_parallax": np.round(1000 * r_np[late], 0).tolist(),
                                     "parallax": np.round(1000 * r_pl[late], 0).tolist()}
    with open(outp, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
