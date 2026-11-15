"""Parallax profile of one OGLE event: the best fit at a fixed parallax magnitude.

Usage: python profile.py DATADIR EVENT OUT.json [--free-blend] [--tmax HJD-2450000]

At each |piE| the direction of the parallax vector, t0, u0 and tE are refitted
from 40 starts (both signs of u0, four start geometries, five start directions)
with the likelihood of candidate.py: the same clipping, uncertainties rescaled
to chi2/dof = 1 at the best parallax fit, and the blend floor of lens.py.
`--free-blend` removes the floor; `--tmax` drops later photometry. Delta chi2
is measured from the best parallax fit of the same data.
"""
import argparse
import json
import pathlib

import numpy as np
from scipy.optimize import minimize

import candidate
import lens

PIE = (0.0, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30)
OPTIONS = {"xatol": 1e-6, "fatol": 1e-5, "maxiter": 8000, "maxfev": 16000}


def minimise(fun, x0):
    first = minimize(fun, x0, method="Nelder-Mead", options=OPTIONS)
    second = minimize(fun, first.x, method="Nelder-Mead", options=OPTIONS)
    return second if second.fun < first.fun else first


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir")
    parser.add_argument("event")
    parser.add_argument("out")
    parser.add_argument("--free-blend", action="store_true")
    parser.add_argument("--tmax", type=float, default=None)
    args = parser.parse_args()
    data = candidate.prepare(pathlib.Path(args.datadir), args.event, args.tmax)
    t, flux, sigma, par = data["t"], data["flux"], data["sigma"], data["par"]
    if args.free_blend:
        lens.BLEND_FLOOR = -1e30

    def chi2(theta):
        return float(np.sum(lens.residuals(np.asarray(theta, float), t, flux, sigma, par) ** 2))

    best = {sign: data["fits"][sign][0] for sign in (1, -1)}
    c_best = min(chi2(best[sign]) for sign in best)
    rows = []
    for pie in PIE:
        found = None
        for sign in (1, -1):
            t0, u0, te, pn, pe = best[sign]
            phi_best = np.arctan2(pe, pn)
            for dt0, fu0, fte in ((0, 1.0, 1.0), (0, 1.6, 0.75), (-15, 0.7, 1.3), (15, 1.0, 0.9)):
                start = [t0 + dt0, u0 * fu0, te * fte]
                if pie == 0.0:
                    r = minimise(lambda q: chi2([q[0], q[1], q[2], 0.0, 0.0]), start)
                    theta = [*r.x, 0.0, 0.0]
                    if found is None or r.fun < found[1]:
                        found = (theta, r.fun)
                    continue
                for phi0 in (phi_best, 0.0, np.pi / 2, np.pi, -np.pi / 2):
                    def fun(q):
                        return chi2([q[0], q[1], q[2], pie * np.cos(q[3]), pie * np.sin(q[3])])
                    r = minimise(fun, start + [phi0])
                    if found is None or r.fun < found[1]:
                        q = r.x
                        found = ([q[0], q[1], q[2], pie * np.cos(q[3]), pie * np.sin(q[3])], r.fun)
        theta, value = found
        _, fs, fb = lens.model(np.asarray(theta), t, flux, sigma, par)
        rows.append({"pie": pie, "delta_chi2": round(value - c_best, 2), "t0": round(theta[0], 3),
                     "u0": round(theta[1], 4), "te": round(theta[2], 1), "pi_n": round(theta[3], 4),
                     "pi_e": round(theta[4], 4), "blend_fraction": round(float(fb / (fs + fb)), 3)})
        print(rows[-1])
    out = {"event": args.event, "n_points": len(t), "tmax": args.tmax, "blend_floor": not args.free_blend,
           "chi2_best_parallax": c_best, "profile": rows}
    with open(args.out, "w", newline="\n") as handle:
        json.dump(out, handle, indent=1)


if __name__ == "__main__":
    main()
