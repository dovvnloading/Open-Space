"""Posterior sampling and lens-mass inference for one OGLE event.

Usage: python candidate.py DATADIR EVENT OUTDIR [--ai A_I] [--steps N] [--tmax HJD-2450000]

Samples (t0, u0, tE, piE_N, piE_E) with emcee for both signs of u0, with flat
priors, linear source and blend fluxes, 5-sigma clipping against the best
parallax fit and uncertainties rescaled to chi2/dof = 1 (then multiplied by
`--inflate`). `--tmax` discards
later photometry so that a fit can be frozen and tested against newer data.
The samples are then weighted with the Galactic prior in galaxy.py.
"""
import argparse
import json
import pathlib

import emcee
import numpy as np

import galaxy
import lens
import screen


def prepare(datadir, event, tmax=None, inflate=1.0):
    year, kind, number = event.split("-")
    row = next(r for r in screen.table(datadir, year) if r["event"] == event)
    t, mag, err, _, _ = lens.read_phot(datadir / year / f"{kind.lower()}-{number}.dat")
    if tmax is not None:
        keep = t <= tmax
        t, mag, err = t[keep], mag[keep], err[keep]
    flux, sigma = lens.to_flux(mag, err)
    ra, dec = lens.sexagesimal(row["ra"], row["dec"])
    first = screen.best_start(row, t, flux, sigma)
    theta0 = first[0].copy()
    theta0[1] = abs(theta0[1])
    par = lens.Parallax(t, ra, dec, theta0[0])
    start = np.concatenate([theta0, [0.0, 0.0]])
    trial = lens.fit(start, t, flux, sigma, par)
    resid = lens.residuals(trial[0], t, flux, sigma, par)[:-1]
    scale = max(1.0, np.sqrt(np.median(resid ** 2) / 0.4549))
    keep = np.abs(resid) < 5 * scale
    t, flux, sigma, mag = t[keep], flux[keep], sigma[keep], mag[keep]
    par = lens.Parallax(t, ra, dec, theta0[0])
    fits = {}
    for sign in (1, -1):
        guess = start.copy()
        guess[1] *= sign
        fits[sign] = lens.fit(guess, t, flux, sigma, par)
    best = min(fits.values(), key=lambda f: f[2])
    sigma = sigma * np.sqrt(best[2] / (len(t) - 7)) * inflate
    fits = {sign: lens.fit(fits[sign][0], t, flux, sigma, par) for sign in fits}
    base = lens.fit(theta0, t, flux, sigma)
    return {"row": row, "t": t, "flux": flux, "sigma": sigma, "ra": ra, "dec": dec, "par": par,
            "t0par": float(theta0[0]), "fits": fits, "base": base, "n_clipped": int((~keep).sum())}


def log_prob(theta, t, flux, sigma, par, sign):
    if theta[2] <= 1 or theta[2] > 5000 or theta[1] * sign < 0 or abs(theta[1]) > 3 or abs(theta[3]) > 5 \
            or abs(theta[4]) > 5:
        return -np.inf
    return -0.5 * np.sum(lens.residuals(theta, t, flux, sigma, par) ** 2)


def sample(data, sign, steps, seed):
    best, err = data["fits"][sign][0], data["fits"][sign][1]
    rng = np.random.default_rng(seed)
    err = np.where(np.isfinite(err) & (err > 0), err, np.abs(best) * 0.01 + 1e-3)
    walkers = 40
    start = best + 0.1 * err * rng.standard_normal((walkers, 5))
    start[:, 1] = np.abs(start[:, 1]) * sign
    start[:, 2] = np.abs(start[:, 2])
    sampler = emcee.EnsembleSampler(walkers, 5, log_prob,
                                    args=(data["t"], data["flux"], data["sigma"], data["par"], sign))
    sampler._random = np.random.RandomState(1000 + seed)      # fixed stream: reruns are identical
    sampler.run_mcmc(start, steps, progress=False)
    chain = sampler.get_chain(discard=steps // 2, thin=10, flat=True)
    logp = sampler.get_log_prob(discard=steps // 2, thin=10, flat=True)
    try:
        tau = float(np.max(sampler.get_autocorr_time(tol=0)))
    except Exception:
        tau = float("nan")
    return chain, logp, tau


def summarize(values, weights=None):
    if weights is None:
        weights = np.ones(len(values))
    low, mid, high = galaxy.quantiles(values, weights)
    return [float(mid), float(mid - low), float(high - mid)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir")
    parser.add_argument("event")
    parser.add_argument("outdir")
    parser.add_argument("--ai", type=float, default=1.5)
    parser.add_argument("--steps", type=int, default=6000)
    parser.add_argument("--tmax", type=float, default=None)
    parser.add_argument("--inflate", type=float, default=1.0,
                        help="extra factor on the rescaled uncertainties (correlated-noise allowance)")
    args = parser.parse_args()
    data = prepare(pathlib.Path(args.datadir), args.event, args.tmax, args.inflate)
    t, flux, sigma, par = data["t"], data["flux"], data["sigma"], data["par"]
    _, _, vn, ve = lens.sun_projected(np.array([data["t0par"]]), data["ra"], data["dec"])
    earth_velocity = (-vn[0] * 365.25, -ve[0] * 365.25)
    out = {"event": args.event, "ra": data["row"]["ra"], "dec": data["row"]["dec"], "n_points": len(t),
           "n_clipped": data["n_clipped"], "t0par": data["t0par"], "tmax": args.tmax,
           "first_hjd": float(t.min()), "last_hjd": float(t.max()),
           "chi2_no_parallax": float(data["base"][2]), "te_no_parallax": float(data["base"][0][2]),
           "a_i_total_assumed": args.ai, "error_inflation": args.inflate, "solutions": {}}
    chains = {}
    for sign in (1, -1):
        chain, logp, tau = sample(data, sign, args.steps, 10 + sign)
        fs, fb = np.empty(len(chain)), np.empty(len(chain))
        for i, theta in enumerate(chain):
            _, fs[i], fb[i] = lens.model(theta, t, flux, sigma, par)
        chains[sign] = (chain, logp, fs, fb)
        best = chain[np.argmax(logp)]
        solution = {"chi2": float(-2 * logp.max()), "autocorr_steps": tau, "best": [float(x) for x in best],
                    "t0": summarize(chain[:, 0]), "u0": summarize(chain[:, 1]), "te": summarize(chain[:, 2]),
                    "pi_n": summarize(chain[:, 3]), "pi_e": summarize(chain[:, 4]),
                    "pie": summarize(np.hypot(chain[:, 3], chain[:, 4])),
                    "pie_95_upper": float(np.quantile(np.hypot(chain[:, 3], chain[:, 4]), 0.95)),
                    "source_i": summarize(lens.ZP - 2.5 * np.log10(np.maximum(fs, 1e-9))),
                    "blend_fraction": summarize(fb / (fs + fb))}
        for slope in (-1.0, -1.75, -2.35):
            post = galaxy.lens_posterior(chain[:, 2], chain[:, 3], chain[:, 4], fb, data["ra"], data["dec"],
                                         earth_velocity, a_i_total=args.ai, mass_slope=slope, seed=5)
            w = post["weight"]
            dark = post["ms_lens_i"] < post["blend_i"]       # a main-sequence lens would outshine the blend
            solution[f"prior_slope_{slope}"] = {
                "mass": summarize(post["mass"], w), "d_l": summarize(post["d_l"], w),
                "mu_hel": summarize(post["mu_hel"], w), "theta_e": summarize(post["theta_e"], w),
                "p_mass_gt_2": float(np.sum(w[post["mass"] > 2]) / w.sum()),
                "p_mass_gt_3": float(np.sum(w[post["mass"] > 3]) / w.sum()),
                "p_ms_brighter_than_blend": float(np.sum(w[dark]) / w.sum()),
                "p_dark_and_gt_2": float(np.sum(w[dark & (post["mass"] > 2)]) / w.sum()),
                "p_dark_and_gt_3": float(np.sum(w[dark & (post["mass"] > 3)]) / w.sum()),
                "disc_fraction": float(np.sum(w * post["disc_fraction"]) / w.sum())}
        out["solutions"]["u0_plus" if sign > 0 else "u0_minus"] = solution
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.tmax is None else f"-to{int(args.tmax)}"
    (outdir / f"{args.event}{suffix}.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(outdir / f"{args.event}{suffix}-chains.npz",
                        plus=chains[1][0][::5], minus=chains[-1][0][::5],
                        plus_logp=chains[1][1][::5], minus_logp=chains[-1][1][::5])
    for name, solution in out["solutions"].items():
        flat, steep = solution["prior_slope_-1.0"], solution["prior_slope_-2.35"]
        print(args.event, name, "chi2", round(solution["chi2"], 1), "(no parallax", round(out["chi2_no_parallax"], 1),
              ") tE", np.round(solution["te"], 1), "piE", np.round(solution["pie"], 4),
              "blend", np.round(solution["blend_fraction"], 2),
              "| M(flat logM)", np.round(flat["mass"], 2), "P(dark&>3)", round(flat["p_dark_and_gt_3"], 3),
              "| M(Salpeter)", np.round(steep["mass"], 2), "P(dark&>3)", round(steep["p_dark_and_gt_3"], 3))


if __name__ == "__main__":
    main()
