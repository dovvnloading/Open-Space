"""Independent point-lens fits, with and without annual parallax, for OGLE-2026-BLG-0058.

Written separately from ../lens.py: the Sun's position comes from astropy
rather than ERFA and the source and blend fluxes are solved in closed form.
The sign conventions of u0 and the parallax vector match MulensModel
(validate.py checks this).

Usage: python -I check/pspl.py <phot.dat> <out.json>
"""
import json
import sys

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time
from scipy.optimize import minimize

RA_DEC = "17:52:38.09 -31:47:36.1"
TPAR = 2461213.0          # parallax reference time, near the peak
ZP = 18.0                 # arbitrary flux zero point (I mag)
KAPPA = 8.144             # mas / Msun

coord = SkyCoord(RA_DEC, unit=(u.hourangle, u.deg))
ra, dec = coord.ra.rad, coord.dec.rad
EAST = np.array([-np.sin(ra), np.cos(ra), 0.0])
NORTH = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])


def sun_projected(times):
    tt = Time(times, format="jd", scale="tdb")
    e = get_body_barycentric("earth", tt).xyz.to(u.au).value
    s = get_body_barycentric("sun", tt).xyz.to(u.au).value
    vec = s - e                      # Sun as seen from Earth
    return NORTH @ vec, EAST @ vec


def parallax_offsets(times, tpar=TPAR):
    s_n, s_e = sun_projected(times)
    ref_n, ref_e = sun_projected(np.array([tpar - 1.0, tpar, tpar + 1.0]))
    v_n = (ref_n[2] - ref_n[0]) / 2.0
    v_e = (ref_e[2] - ref_e[0]) / 2.0
    return (s_n - ref_n[1] - (times - tpar) * v_n,
            s_e - ref_e[1] - (times - tpar) * v_e)


def magnification(params, t, ds_n, ds_e):
    t0, u0, te, pn, pe = params
    tau = (t - t0) / te + pn * ds_n + pe * ds_e
    beta = -u0 - pn * ds_e + pe * ds_n      # sign chosen to match MulensModel's u0 convention
    uu = np.sqrt(tau ** 2 + beta ** 2)
    return (uu ** 2 + 2) / (uu * np.sqrt(uu ** 2 + 4))


def linear_fluxes(a, f, w):
    # weighted least squares for f = fs * a + fb
    s_w, s_a, s_f = w.sum(), (w * a).sum(), (w * f).sum()
    s_aa, s_af = (w * a * a).sum(), (w * a * f).sum()
    det = s_w * s_aa - s_a ** 2
    fs = (s_w * s_af - s_a * s_f) / det
    fb = (s_aa * s_f - s_a * s_af) / det
    cov = np.array([[s_w, -s_a], [-s_a, s_aa]]) / det   # var(fs), cov, var(fb)
    return fs, fb, cov


class Fitter:
    def __init__(self, t, f, sf):
        self.t, self.f, self.w = t, f, 1.0 / sf ** 2
        self.ds_n, self.ds_e = parallax_offsets(t)

    def chi2(self, params):
        if params[2] <= 1.0:
            return 1e30
        a = magnification(params, self.t, self.ds_n, self.ds_e)
        fs, fb, _ = linear_fluxes(a, self.f, self.w)
        r = self.f - fs * a - fb
        return float((self.w * r * r).sum())

    def fit(self, start, free_parallax=True, fixed_pi=(0.0, 0.0)):
        if free_parallax:
            fun = self.chi2
            x0 = list(start)
        else:
            fun = lambda p: self.chi2([p[0], p[1], p[2], fixed_pi[0], fixed_pi[1]])
            x0 = list(start[:3])
        best = None
        for scale in (1.0, 0.5):
            res = minimize(fun, x0, method="Nelder-Mead",
                           options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000, "maxfev": 40000})
            x0 = res.x
            if best is None or res.fun < best.fun:
                best = res
        x = list(best.x) if free_parallax else list(best.x) + list(fixed_pi)
        return x, best.fun

    def details(self, params):
        a = magnification(params, self.t, self.ds_n, self.ds_e)
        fs, fb, cov = linear_fluxes(a, self.f, self.w)
        return a, fs, fb, cov


def mag_of(flux):
    return ZP - 2.5 * np.log10(flux) if flux > 0 else float("nan")


def main():
    data = np.loadtxt(sys.argv[1])
    t, mag, err = data[:, 0], data[:, 1], data[:, 2]
    f = 10 ** (-0.4 * (mag - ZP))
    sf = f * err * np.log(10) / 2.5
    out = {"n_points": int(len(t)), "t_first": float(t.min()), "t_last": float(t.max())}

    fit = Fitter(t, f, sf)
    # no-parallax fit
    x_np, c_np = fit.fit([2461213.4, 0.5, 229.0, 0.0, 0.0], free_parallax=False)
    out["no_parallax"] = {"params": x_np, "chi2": c_np}
    # parallax fits, both signs of u0, several starts
    sols = {}
    for sign in (+1, -1):
        best = None
        for pn0, pe0 in [(0.05, 0.0), (-0.05, 0.0), (0.0, 0.05), (0.0, -0.05), (0.1, 0.1), (-0.1, -0.1)]:
            x, c = fit.fit([x_np[0], sign * abs(x_np[1]), x_np[2], pn0, pe0])
            if best is None or c < best[1]:
                best = (x, c)
        sols["u0+" if sign > 0 else "u0-"] = best
    dof = len(t) - 7
    for key, (x, c) in sols.items():
        a, fs, fb, cov = fit.details(x)
        pie = float(np.hypot(x[3], x[4]))
        out[key] = {"params": x, "chi2": c, "piE": pie, "fs": fs, "fb": fb,
                    "I_source": mag_of(fs), "I_blend": mag_of(fb),
                    "fb_sigma": float(np.sqrt(cov[1, 1])),
                    "blend_fraction_of_baseline": fb / (fs + fb)}
    a, fs, fb, cov = fit.details(x_np)
    out["no_parallax"].update({"fs": fs, "fb": fb, "I_source": mag_of(fs), "I_blend": mag_of(fb),
                               "fb_sigma": float(np.sqrt(cov[1, 1]))})
    best_key = min(sols, key=lambda k: sols[k][1])
    c_best = sols[best_key][1]
    out["best"] = best_key
    out["error_scale"] = float(np.sqrt(c_best / dof))
    out["delta_chi2_parallax_raw"] = c_np - c_best
    out["delta_chi2_parallax_rescaled"] = (c_np - c_best) / (c_best / dof)

    with open(sys.argv[2], "w", newline="\n") as fh:
        json.dump(out, fh, indent=1, default=float)

    print(json.dumps({k: out[k] for k in ("n_points", "best", "error_scale",
                                          "delta_chi2_parallax_raw", "delta_chi2_parallax_rescaled")},
                     default=float))
    for key in ("no_parallax", "u0+", "u0-"):
        o = out[key]
        p = o["params"]
        print(key, "t0=%.2f u0=%.4f tE=%.1f piEN=%.4f piEE=%.4f chi2=%.1f I_s=%.3f I_b=%s fb=%.4f+-%.4f"
              % (p[0], p[1], p[2], p[3], p[4], o["chi2"], o["I_source"],
                 ("%.2f" % o["I_blend"]) if o["I_blend"] == o["I_blend"] else "neg", o["fb"], o["fb_sigma"]))


if __name__ == "__main__":
    main()
