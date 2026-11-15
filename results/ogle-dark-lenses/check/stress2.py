"""Follow-up tests for OGLE-2026-BLG-0058.

(a) Holdout: fit only the data before the 71-day gap in July-September 2026,
    then predict the seven later points.
(b) The parallax upper limit without those points, with the blend floor.
(c) The other analysis's own likelihood code evaluated at this script's
    constrained optima, to separate optimiser failures from model differences.
(d) Forecast differences between the models through the 2027 season.

Usage: python -I check/stress2.py <phot.dat> <fits.json> <stress.json> <../lens.py> <out.json>
"""
import importlib.util
import json
import os
import sys

import numpy as np
from astropy.time import Time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pspl as lens  # noqa: E402
import stress as S  # noqa: E402


def main():
    phot, fitjson, sjson, their_path, outp = sys.argv[1:6]
    t, mag, err, _, _ = np.loadtxt(phot).T
    prev = json.load(open(fitjson))
    st = json.load(open(sjson))
    sols = {k: prev[k]["params"] for k in ("no_parallax", "u0+", "u0-")}
    f, sf = S.to_flux(mag, err)
    sfr = sf * st["error_scale"]
    s26 = np.where(t > 2461040.0)[0]
    g = int(np.argmax(np.diff(t[s26])))
    late = t > t[s26][g] + 0.5
    keep = ~late
    out = {"n_late": int(late.sum())}

    # (a) holdout prediction of the late points, with and without the blend floor
    dn_l, de_l = lens.parallax_offsets(t[late])
    for label, floor in (("free_blend", None), ("blend_floor", S.FLOOR)):
        fx = S.F(t[keep], f[keep], sfr[keep], floor=floor)
        res = S.fit_set(fx, sols)
        d = S.describe(res, fx)
        for key, (p, _) in res.items():
            _, fs, fb = fx.fluxes(p)
            model = fs * lens.magnification(p, t[late], dn_l, de_l) + fb
            z = (f[late] - model) / sfr[late]
            d[key]["holdout_chi2_late"] = float((z ** 2).sum())
            d[key]["holdout_mean_mag_offset"] = float(np.mean(mag[late] - (lens.ZP - 2.5 * np.log10(model))))
        out["pre_gap_" + label] = d
        if floor is not None:
            best = min(res["u0+"][1], res["u0-"][1])
            scan = {}
            for pie in (0.10, 0.15, 0.20):
                r = S.fit_fixed_pie(fx, pie, sols["u0+"][:3])
                r["delta_chi2"] = r["chi2"] - best
                scan[str(pie)] = r
            out["pre_gap_blend_floor_fixed_piE"] = scan

    # (c) the other analysis's likelihood at these optima
    spec = importlib.util.spec_from_file_location("their_lens", their_path)
    tl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tl)
    tt, mm, ee, _, _ = tl.read_phot(phot)
    flux, sigma = tl.to_flux(mm, ee)
    ra, dec = tl.sexagesimal("17:52:38.09", "-31:47:36.1")
    par = tl.Parallax(tt, ra, dec, lens.TPAR - 2450000.0)

    def chi2(theta, sig):
        return float(np.sum(tl.residuals(np.asarray(theta, float), tt, flux, sig, par) ** 2))

    def conv(p):
        return [p[0] - 2450000.0, p[1], p[2], p[3], p[4]]

    fits = {key: tl.fit(np.array(conv(sols[key])), tt, flux, sigma, par) for key in ("u0+", "u0-")}
    sig = sigma * np.sqrt(min(v[2] for v in fits.values()) / (len(tt) - 7))
    c0 = min(chi2(v[0], sig) for v in fits.values())
    theirs = {"n_points": int(len(tt)), "best_chi2_rescaled": c0,
              "no_parallax_delta": chi2(tl.fit(np.array(conv(sols["no_parallax"])[:3]), tt, flux, sig)[0][:3].tolist()
                                         + [0.0, 0.0], sig) - c0}
    for key, m in st["fixed_piE_scan"]["blend_floor"].items():
        theirs["delta_at_this_piE_" + key + "_optimum"] = chi2(
            [m["t0"] - 2450000.0, m["u0"], m["tE"], m["piEN"], m["piEE"]], sig) - c0
    out["their_likelihood"] = theirs

    # (d) forecast differences, I(parallax) - I(no parallax)
    dates = [f"2027-{mo:02d}-{dy:02d}" for mo in range(2, 11) for dy in (1, 15)]
    jd = Time(dates).jd
    dn, de = lens.parallax_offsets(jd)
    full = S.F(t, f, sfr)
    curve = {}
    for key, p in sols.items():
        _, fs, fb = full.fluxes(p)
        curve[key] = lens.ZP - 2.5 * np.log10(fs * lens.magnification(p, jd, dn, de) + fb)
    out["forecast"] = {d: {"no_parallax_I": round(float(curve["no_parallax"][i]), 3),
                           "u0+_minus_np": round(float(curve["u0+"][i] - curve["no_parallax"][i]), 3),
                           "u0-_minus_np": round(float(curve["u0-"][i] - curve["no_parallax"][i]), 3)}
                       for i, d in enumerate(dates)}
    with open(outp, "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
