"""Fit every OGLE-IV EWS event with a point lens, with and without annual parallax.

Usage: python screen.py DATADIR OUT.csv [YEAR ...]

For each event the script fits (t0, u0, tE) and then adds the parallax vector
for both signs of u0. Source and blend fluxes are solved linearly. Points more
than five scaled sigma from the first parallax fit are removed once and the
uncertainties are rescaled so that the final best fit has chi2/dof = 1.
"""
import csv
import multiprocessing as mp
import pathlib
import sys

import numpy as np

import lens

FIELDS = ["event", "ra", "dec", "n", "n_clipped", "span_days", "last_hjd", "ews_tau", "ews_umin",
          "errscale", "chi2_pspl", "t0_0", "u0_0", "te_0", "fs_0", "fb_0",
          "chi2_par", "t0", "u0", "te", "pi_n", "pi_e", "e_t0", "e_u0", "e_te", "e_pi_n", "e_pi_e",
          "fs", "fb", "chi2_alt", "te_alt", "u0_alt", "pi_n_alt", "pi_e_alt", "e_pi_n_alt", "e_pi_e_alt",
          "fs_alt", "fb_alt", "dchi2", "i_base", "status"]


def table(datadir, year):
    for line in (datadir / year / "lenses.par").read_text().splitlines():
        p = line.split()
        if len(p) >= 14 and p[0][:4].isdigit():
            yield {"event": p[0], "ra": p[3], "dec": p[4], "tmax": float(p[5]) - 2450000,
                   "tau": float(p[7]), "umin": float(p[8])}


def best_start(row, t, flux, sigma):
    """Try the EWS solution and a small grid, return the best no-parallax fit."""
    best = None
    peak = t[np.argmax(flux)]
    for t0 in (row["tmax"], peak):
        for te in (row["tau"], 20.0, 80.0, 250.0):
            for u0 in (max(min(row["umin"], 1.5), 0.01), 0.1, 0.6):
                if not (0.5 < te < 3000):
                    continue
                try:
                    out = lens.fit(np.array([t0, u0, te]), t, flux, sigma)
                except Exception:
                    continue
                if np.isfinite(out[2]) and out[0][2] > 0 and (best is None or out[2] < best[2]):
                    best = out
    return best


def one(args):
    datadir, row = args
    year, kind, number = row["event"].split("-")
    out = {"event": row["event"], "ra": row["ra"], "dec": row["dec"], "ews_tau": row["tau"],
           "ews_umin": row["umin"], "status": "ok"}
    try:
        t, mag, err, _, _ = lens.read_phot(datadir / year / f"{kind.lower()}-{number}.dat")
        if len(t) < 60:
            out["status"] = "few points"
            return out
        flux, sigma = lens.to_flux(mag, err)
        ra, dec = lens.sexagesimal(row["ra"], row["dec"])
        first = best_start(row, t, flux, sigma)
        if first is None:
            out["status"] = "no fit"
            return out
        theta0 = first[0]
        theta0[1] = abs(theta0[1])
        par = lens.Parallax(t, ra, dec, theta0[0])
        start = np.concatenate([theta0, [0.0, 0.0]])
        trial = lens.fit(start, t, flux, sigma, par)
        resid = lens.residuals(trial[0], t, flux, sigma, par)[:-1]
        scale = max(1.0, np.sqrt(np.median(resid ** 2) / 0.4549))
        keep = np.abs(resid) < 5 * scale
        out["n"], out["n_clipped"] = int(keep.sum()), int((~keep).sum())
        t, flux, sigma = t[keep], flux[keep], sigma[keep]
        par = lens.Parallax(t, ra, dec, theta0[0])
        plus = lens.fit(start, t, flux, sigma, par)
        minus_start = start.copy()
        minus_start[1] *= -1
        minus = lens.fit(minus_start, t, flux, sigma, par)
        best, alt = (plus, minus) if plus[2] <= minus[2] else (minus, plus)
        k2 = best[2] / (len(t) - 7)
        sigma = sigma * np.sqrt(k2)
        base = lens.fit(theta0, t, flux, sigma)
        _, fs0, fb0 = lens.model(base[0], t, flux, sigma)
        _, fs, fb = lens.model(best[0], t, flux, sigma, par)
        _, fsa, fba = lens.model(alt[0], t, flux, sigma, par)
        e = best[1] * np.sqrt(k2)
        ea = alt[1] * np.sqrt(k2)
        out.update(span_days=float(np.ptp(t)), last_hjd=float(t.max()), errscale=float(np.sqrt(k2)),
                   chi2_pspl=base[2], t0_0=base[0][0], u0_0=base[0][1], te_0=base[0][2], fs_0=fs0, fb_0=fb0,
                   chi2_par=best[2] / k2, t0=best[0][0], u0=best[0][1], te=best[0][2], pi_n=best[0][3],
                   pi_e=best[0][4], e_t0=e[0], e_u0=e[1], e_te=e[2], e_pi_n=e[3], e_pi_e=e[4], fs=fs, fb=fb,
                   chi2_alt=alt[2] / k2, te_alt=alt[0][2], u0_alt=alt[0][1], pi_n_alt=alt[0][3],
                   pi_e_alt=alt[0][4], e_pi_n_alt=ea[3], e_pi_e_alt=ea[4], fs_alt=fsa, fb_alt=fba,
                   dchi2=base[2] - best[2] / k2,
                   i_base=float(lens.ZP - 2.5 * np.log10(max(fs + fb, 1e-9))))
    except Exception as error:
        out["status"] = f"error {type(error).__name__}: {error}"
    return out


def main():
    datadir = pathlib.Path(sys.argv[1])
    years = sys.argv[3:] or ["2022", "2023", "2024", "2025", "2026"]
    rows = [(datadir, r) for y in years for r in table(datadir, y)]
    with mp.Pool(10) as pool, open(sys.argv[2], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for i, out in enumerate(pool.imap(one, rows, chunksize=4)):
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in out.items()})
            if i % 500 == 0:
                print(i, flush=True)


if __name__ == "__main__":
    main()
