"""Stress tests of the parallax signal in OGLE-2026-BLG-0058.

Usage: python -I check/stress.py <phot.dat> <fits.json> <out.json>

All chi2 values are in units of errors rescaled so the best unconstrained
parallax fit to the full data set has chi2/dof = 1.
"""
import json
import os
import sys

import numpy as np
from astropy.time import Time
from scipy.optimize import minimize
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pspl as lens  # noqa: E402

NM = {"xatol": 1e-6, "fatol": 1e-5, "maxiter": 8000, "maxfev": 16000}
FLOOR = -0.158      # flux of a "negative" I = 20 star at ZP 18, the floor used by the result under test
SOFT = 0.03


def to_flux(mag, err, floor=0.0):
    f = 10 ** (-0.4 * (mag - lens.ZP))
    return f, f * np.log(10) / 2.5 * np.hypot(err, floor)


class F:
    """chi2 of one data subset, with an optional blend-flux floor or zero blend."""

    def __init__(self, t, f, sf, floor=None, no_blend=False):
        self.t, self.f, self.w = t, f, 1.0 / sf ** 2
        self.dn, self.de = lens.parallax_offsets(t)
        self.floor, self.no_blend = floor, no_blend

    def fluxes(self, p):
        a = lens.magnification(p, self.t, self.dn, self.de)
        if self.no_blend:
            return a, (self.w * a * self.f).sum() / (self.w * a * a).sum(), 0.0
        fs, fb, _ = lens.linear_fluxes(a, self.f, self.w)
        return a, fs, fb

    def chi2(self, p):
        if p[2] <= 1.0 or p[2] > 5000:
            return 1e30
        a, fs, fb = self.fluxes(p)
        r = self.f - fs * a - fb
        c = float((self.w * r * r).sum())
        if self.floor is not None and fb < self.floor:
            c += ((fb - self.floor) / SOFT) ** 2
        return c

    def resid(self, p):
        a, fs, fb = self.fluxes(p)
        return (self.f - fs * a - fb) * np.sqrt(self.w)


def nm(fun, x0):
    r1 = minimize(fun, x0, method="Nelder-Mead", options=NM)
    r2 = minimize(fun, r1.x, method="Nelder-Mead", options=NM)
    return r2 if r2.fun < r1.fun else r1


def fit_set(fx, sols):
    """No-parallax fit plus the best parallax fit for each sign of u0."""
    out = {}
    r = nm(lambda p: fx.chi2([p[0], p[1], p[2], 0.0, 0.0]), list(sols["no_parallax"][:3]))
    out["no_parallax"] = (list(r.x) + [0.0, 0.0], float(r.fun))
    for key in ("u0+", "u0-"):
        best = None
        s = sols[key]
        for start in (s, [s[0], s[1], s[2], 0.0, 0.0], [s[0], s[1], s[2], -s[3], -s[4]]):
            r = nm(fx.chi2, list(start))
            if best is None or r.fun < best[1]:
                best = (list(r.x), float(r.fun))
        out[key] = best
    return out


def describe(res, fx, scale=1.0):
    d = {}
    for key, (p, c) in res.items():
        _, fs, fb = fx.fluxes(p)
        d[key] = {"t0": p[0], "u0": p[1], "tE": p[2], "piEN": p[3], "piEE": p[4],
                  "piE": float(np.hypot(p[3], p[4])), "chi2": c / scale,
                  "blend_share": float(fb / (fs + fb))}
    best = min(res["u0+"][1], res["u0-"][1])
    d["delta_chi2_parallax"] = (res["no_parallax"][1] - best) / scale
    d["n"] = int(len(fx.t))
    return d


def fit_fixed_pie(fx, pie, base):
    """Best chi2 at fixed |piE| with the direction, t0, u0 and tE free, over both signs of u0."""
    best = None
    t0b, u0b, teb = base
    for sign in (1, -1):
        for dt0, fu0, fte in ((0, 1.0, 1.0), (0, 1.6, 0.75), (-15, 0.7, 1.3), (15, 1.0, 0.9)):
            for phi0 in np.radians([0, 90, 180, 270]):
                def fun(q):
                    return fx.chi2([q[0], q[1], q[2], pie * np.cos(q[3]), pie * np.sin(q[3])])
                r = nm(fun, [t0b + dt0, sign * abs(u0b) * fu0, teb * fte, phi0])
                if best is None or r.fun < best[1]:
                    q = r.x
                    best = ([q[0], q[1], q[2], pie * np.cos(q[3]), pie * np.sin(q[3])], float(r.fun))
    p, c = best
    _, fs, fb = fx.fluxes(p)
    return {"chi2": c, "t0": p[0], "u0": p[1], "tE": p[2], "piEN": p[3], "piEE": p[4],
            "fb": float(fb), "blend_share": float(fb / (fs + fb))}


def iso(t):
    return Time(t, format="jd").iso[:10]


def main():
    phot, fitjson, outpath = sys.argv[1:4]
    t, mag, err, see, sky = np.loadtxt(phot).T
    prev = json.load(open(fitjson))
    sols = {k: prev[k]["params"] for k in ("no_parallax", "u0+", "u0-")}
    f, sf = to_flux(mag, err)
    raw = F(t, f, sf)
    c_best = min(raw.chi2(sols["u0+"]), raw.chi2(sols["u0-"]))
    k = float(np.sqrt(c_best / (len(t) - 7)))
    sfr = sf * k
    R = F(t, f, sfr)
    out = {"error_scale": k, "full": {key: R.chi2(p) for key, p in sols.items()}}

    # 1. points after the long gap at the end of the 2026 season
    s26 = np.where(t > 2461040.0)[0]
    g = int(np.argmax(np.diff(t[s26])))
    gap0, gap1 = t[s26][g], t[s26][g + 1]
    late = t > gap0 + 0.5
    rn, rp, rm = R.resid(sols["no_parallax"]), R.resid(sols["u0+"]), R.resid(sols["u0-"])
    in26 = t > 2461040.0
    out["late"] = {
        "gap_from": iso(gap0), "gap_to": iso(gap1), "gap_days": float(gap1 - gap0), "n_late": int(late.sum()),
        "dchi2_late_u0+": float((rn[late] ** 2 - rp[late] ** 2).sum()),
        "dchi2_late_u0-": float((rn[late] ** 2 - rm[late] ** 2).sum()),
        "dchi2_all_u0+": float((rn ** 2 - rp ** 2).sum()),
        "season26_median_seeing": float(np.median(see[in26 & ~late])),
        "season26_median_sky": float(np.median(sky[in26 & ~late])),
        "rows": [[iso(t[i]), float(mag[i]), float(err[i]), float(see[i]), float(sky[i]),
                  round(float(rn[i]), 2), round(float(rp[i]), 2)] for i in np.where(late)[0]]}
    # cumulative parallax gain by season
    yrs = Time(t, format="jd").decimalyear.astype(int)
    out["dchi2_by_year_u0+"] = {int(y): round(float((rn[yrs == y] ** 2 - rp[yrs == y] ** 2).sum()), 2)
                               for y in np.unique(yrs)}

    # 2. refit without the post-gap points
    keep = ~late
    fx = F(t[keep], f[keep], sfr[keep])
    out["without_late"] = describe(fit_set(fx, sols), fx)

    # 3. season jackknife
    jack = {}
    for y in np.unique(yrs):
        keep = yrs != y
        fx = F(t[keep], f[keep], sfr[keep])
        d = describe(fit_set(fx, sols), fx)
        jack[int(y)] = {"n_dropped": int((~keep).sum()), "delta_chi2_parallax": round(d["delta_chi2_parallax"], 2),
                        "u0+": [round(d["u0+"]["tE"], 1), round(d["u0+"]["piE"], 4)],
                        "u0-": [round(d["u0-"]["tE"], 1), round(d["u0-"]["piE"], 4)]}
    out["jackknife_year"] = jack

    # 4. the error model of the result under test: 0.004 mag floor and one 5-sigma robust clip
    f4, sf4 = to_flux(mag, err, floor=0.004)
    r4 = F(t, f4, sf4).resid(sols["u0+"])
    clip = np.abs(r4) < 5 * max(1.0, np.sqrt(np.median(r4 ** 2) / 0.4549))
    fx = F(t[clip], f4[clip], sf4[clip])
    res = fit_set(fx, sols)
    scale = min(res["u0+"][1], res["u0-"][1]) / (clip.sum() - 7)
    out["their_error_model"] = describe(res, fx, scale)
    out["their_error_model"]["n_clipped"] = int((~clip).sum())

    # 5. residuals against seeing and sky, 2025-2026
    ev = t > 2460676.0
    out["residual_correlations_2025_2026"] = {
        name: {"seeing": [round(x, 4) for x in spearmanr(r[ev], see[ev])],
               "sky": [round(x, 4) for x in spearmanr(r[ev], sky[ev])]}
        for name, r in (("no_parallax", rn), ("u0+", rp))}

    # 6. fixed |piE| with the direction free, with and without the blend floor
    scan = {}
    for label, floor in (("free_blend", None), ("blend_floor", FLOOR)):
        fx = F(t, f, sfr, floor=floor)
        rows = {}
        for pie in (0.08, 0.10, 0.15, 0.20, 0.30):
            r = fit_fixed_pie(fx, pie, sols["u0+"][:3])
            r["delta_chi2"] = r["chi2"] - out["full"]["u0+"]
            rows[str(pie)] = r
        scan[label] = rows
    out["fixed_piE_scan"] = scan

    # 7. no blend light at all
    fx = F(t, f, sfr, no_blend=True)
    out["zero_blend"] = describe(fit_set(fx, sols), fx)

    # 8. forecasts from the best fits
    when = {"2027-02-15": Time("2027-02-15").jd, "2027-03-01": Time("2027-03-01").jd,
            "2027-04-01": Time("2027-04-01").jd}
    tt = np.array(list(when.values()))
    dn, de = lens.parallax_offsets(tt)
    fc = {}
    for key, p in sols.items():
        _, fs, fb = R.fluxes(p)
        a = lens.magnification(p, tt, dn, de)
        fc[key] = dict(zip(when, [round(float(m), 3) for m in lens.ZP - 2.5 * np.log10(fs * a + fb)]))
    out["forecast_I"] = fc

    with open(outpath, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
