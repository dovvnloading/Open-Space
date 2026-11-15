"""Independent emcee posterior for OGLE-2026-BLG-0058, both signs of u0.

Errors are rescaled so the best parallax fit has chi2/dof = 1. A soft floor
keeps the blend flux above that of a negative I = 20 star (a physical prior,
the same one the result under test uses); pass --free-blend to drop it.

Usage: python -I check/mcmc.py <phot.dat> <fits.json> <out.npz> [--free-blend] [--tmax=HJD] [--inflate=X]

--tmax keeps photometry up to that HJD; --inflate multiplies the rescaled
uncertainties. A JSON summary is written next to the chains.
"""
import json
import os
import sys

import emcee
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pspl as lens  # noqa: E402

FLOOR, SOFT = -0.158, 0.03
WALKERS, STEPS, BURN, THIN = 32, 9000, 3000, 10


def main():
    phot, fitjson, outp = sys.argv[1:4]
    free = "--free-blend" in sys.argv
    opts = dict(a.split("=") for a in sys.argv[4:] if "=" in a)
    t, mag, err = np.loadtxt(phot, usecols=(0, 1, 2)).T
    if "--tmax" in opts:
        keep = t <= float(opts["--tmax"])
        t, mag, err = t[keep], mag[keep], err[keep]
    inflate = float(opts.get("--inflate", 1.0))
    f = 10 ** (-0.4 * (mag - lens.ZP))
    sf = f * np.log(10) / 2.5 * err
    prev = json.load(open(fitjson))
    k2 = min(prev["u0+"]["chi2"], prev["u0-"]["chi2"]) / (len(t) - 7)
    w = 1.0 / (sf ** 2 * k2 * inflate ** 2)
    dn, de = lens.parallax_offsets(t)

    def log_prob(p, sign):
        if p[2] <= 1 or p[2] > 5000 or p[1] * sign <= 0 or abs(p[1]) > 3 or abs(p[3]) > 5 or abs(p[4]) > 5:
            return -np.inf
        a = lens.magnification(p, t, dn, de)
        fs, fb, _ = lens.linear_fluxes(a, f, w)
        r = f - fs * a - fb
        c = float((w * r * r).sum())
        if not free and fb < FLOOR:
            c += ((fb - FLOOR) / SOFT) ** 2
        return -0.5 * c

    out, summary = {}, {"free_blend": free, "options": opts}
    for key, sign in (("plus", 1), ("minus", -1)):
        best = np.array(prev["u0+" if sign > 0 else "u0-"]["params"])
        rng = np.random.default_rng(20 + sign)
        p0 = best + np.array([1.0, 0.02, 10.0, 0.005, 0.005]) * rng.standard_normal((WALKERS, 5))
        p0[:, 1] = sign * np.abs(p0[:, 1])
        sampler = emcee.EnsembleSampler(WALKERS, 5, log_prob, args=(sign,))
        sampler._random = np.random.RandomState(300 + sign)
        sampler.run_mcmc(p0, STEPS, progress=False)
        ch = sampler.get_chain(discard=BURN, thin=THIN, flat=True)
        lp = sampler.get_log_prob(discard=BURN, thin=THIN, flat=True)
        fs, fb = np.empty(len(ch)), np.empty(len(ch))
        for i, p in enumerate(ch):
            fs[i], fb[i], _ = lens.linear_fluxes(lens.magnification(p, t, dn, de), f, w)
        tau = sampler.get_autocorr_time(tol=0)
        out.update({key: ch, key + "_logp": lp, key + "_fs": fs, key + "_fb": fb, key + "_tau": tau})
        pie = np.hypot(ch[:, 3], ch[:, 4])

        def q(x):
            return np.round(np.percentile(x, [16, 50, 84]), 4).tolist()
        summary[key] = {
            "tE": q(ch[:, 2]), "u0": q(ch[:, 1]), "piEN": q(ch[:, 3]), "piEE": q(ch[:, 4]), "piE": q(pie),
            "piE_95": round(float(np.percentile(pie, 95)), 4), "piE_99": round(float(np.percentile(pie, 99)), 4),
            "I_source": q(lens.ZP - 2.5 * np.log10(np.maximum(fs, 1e-9))), "blend_share": q(fb / (fs + fb)),
            "P(blend<0)": round(float(np.mean(fb < 0)), 3),
            "mass_per_masyr": q(ch[:, 2] / 365.25 / (lens.KAPPA * pie)),
            "chi2_min": round(float(-2 * lp.max()), 2), "autocorr": np.round(tau, 0).tolist(),
            "steps_over_tau": round(float(STEPS / tau.max()), 1)}
        print(key, json.dumps(summary[key]))
    np.savez_compressed(outp, **out)
    with open(outp.rsplit(".", 1)[0] + ".json", "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


if __name__ == "__main__":
    main()
