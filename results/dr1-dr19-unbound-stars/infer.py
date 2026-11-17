"""Posterior phase space and escape fractions. Distances in kpc, speeds in km/s."""

import hashlib
from pathlib import Path

import astropy.units as u
import numpy as np
from gala.potential import MilkyWayPotential
from galpy.potential.interpRZPotential import eval_potential_c
from galpy.potential.mwpotentials import McMillan17
from scipy.integrate import cumulative_trapezoid

from make_queue import number, rows, save
from screen_desi import cartesian

GALAXY = MilkyWayPotential(version="v2")
N = 4096
NAMES = ("ra", "dec", "parallax", "pmra", "pmdec")


def draws(g, s, model, n=N):
    sid = s["source_id"]
    seed = int.from_bytes(
        hashlib.sha256((sid + "/" + model).encode()).digest()[:8], "big"
    )
    rng = np.random.default_rng(seed)
    err = np.array([number(g[k + "_error"]) for k in NAMES])
    cov = np.diag(err**2)
    for i, k in enumerate(NAMES):
        for j, l in enumerate(NAMES[i + 1 :], i + 1):
            cov[i, j] = cov[j, i] = err[i] * err[j] * number(g[k + "_" + l + "_corr"])
    cov += np.diag([0, 0, 0.02**2, 0.02**2, 0.02**2])
    p = number(g["parallax"]) - number(s["zp"])
    ep = np.sqrt(cov[2, 2])
    spec = number(s.get("spec_d"))
    es = number(s.get("spec_e"))
    conditional = True
    if model in ("bj_geo", "bj_photogeo", "starhorse", "minesweeper", "gspphot"):
        if model in ("starhorse", "minesweeper", "gspphot"):
            med, lo, hi = [number(s.get(model + "_" + q)) for q in ("med", "lo", "hi")]
        else:
            typ = model[3:]
            med, lo, hi = [
                number(s.get("r_" + q + "_" + typ)) / 1000 for q in ("med", "lo", "hi")
            ]
        if not (0 < lo < med < hi):
            return None
        z = rng.normal(size=n)
        d = med * np.exp(z * np.where(z < 0, np.log(med / lo), np.log(hi / med)))
    else:
        grid = np.linspace(np.log(0.002), np.log(500.0), 12289)
        dist = np.exp(grid)
        if model in ("spec_only", "spec_plx", "spdist_only", "dwarf_only"):
            if model == "dwarf_only":
                spec, es = number(s.get("dwarf_d")), number(s.get("dwarf_e"))
            if model == "spdist_only":
                spec, es = number(s.get("sp_d")), number(s.get("sp_e"))
            if not (spec > 0 and es > 0):
                return None
            logsig = np.sqrt(np.log1p((es / spec) ** 2))
            lp = -0.5 * ((grid - np.log(spec)) / logsig) ** 2
            conditional = model == "spec_plx"
        elif model == "edsd_plx":
            lp = 3 * grid - dist / 1.35
        else:
            raise ValueError(model)
        if conditional:
            lp -= 0.5 * ((p - 1 / dist) / ep) ** 2
        cdf = cumulative_trapezoid(np.exp(lp - lp.max()), grid, initial=0)
        d = np.exp(np.interp(rng.random(n) * cdf[-1], cdf, grid))
    keep = [0, 1, 3, 4]
    cov4 = cov[np.ix_(keep, keep)]
    delta = np.zeros((n, 4))
    if conditional:
        cross = cov[keep, 2]
        cov4 -= np.outer(cross, cross) / cov[2, 2]
        delta += (1 / d - p)[:, None] * cross / cov[2, 2]
    delta += rng.multivariate_normal(np.zeros(4), cov4, n, check_valid="raise")
    ra = number(g["ra"]) + delta[:, 0] / (3.6e6 * np.cos(np.deg2rad(number(g["dec"]))))
    dec = number(g["dec"]) + delta[:, 1] / 3.6e6
    rv = rng.normal(
        number(s["rv"]),
        np.hypot(number(s["e_rv"]), number(s.get("rv_floor") or 1.0)),
        n,
    )
    x, v = cartesian(
        ra,
        dec,
        d,
        number(g["pmra"]) + delta[:, 2],
        number(g["pmdec"]) + delta[:, 3],
        rv,
    )
    return d, x, v


def escape(x):
    R = np.hypot(x[0], x[1])
    z = x[2]
    pot = McMillan17
    phi, flag = eval_potential_c(pot, R / pot._ro, z / pot._ro)
    if flag:
        raise RuntimeError("galpy failed")
    # Both mass models use zero potential at infinity. Work in (km/s)^2.
    mc = np.sqrt(-2 * phi) * pot._vo
    gala = np.sqrt(-2 * GALAXY.energy(x * u.kpc).to_value((u.km / u.s) ** 2))
    return mc, gala


def summarize(s, g):
    out = []
    for model in (
        "spec_only",
        "spec_plx",
        "spdist_only",
        "dwarf_only",
        "edsd_plx",
        "bj_geo",
        "bj_photogeo",
        "starhorse",
        "minesweeper",
        "gspphot",
    ):
        sample = draws(g, s, model)
        if sample is None:
            continue
        d, x, v = sample
        if not np.isfinite(x).all() or not np.isfinite(v).all():
            raise ValueError((s["source_id"], model, "non-finite state"))
        speed = np.linalg.norm(v, axis=0)
        mc, ga = escape(x)
        line = dict(
            source_id=s["source_id"],
            channel=s["channel"],
            model=model,
            n=len(d),
            ra=g["ra"],
            dec=g["dec"],
            rv=s["rv"],
            e_rv=s["e_rv"],
        )
        for k, a in zip(
            ("d", "x", "y", "z", "vx", "vy", "vz", "v_gal", "vesc_mc17", "vesc_gala22"),
            (d, *x, *v, speed, mc, ga),
        ):
            for q, val in zip(("16", "50", "84"), np.percentile(a, [16, 50, 84])):
                line[k + "_" + q] = round(float(val), 5)
        for k, a in [
            ("p_mc17", speed > mc),
            ("p_gala22", speed > ga),
            ("p_v1000", speed > 1000),
        ]:
            line[k] = float(np.mean(a))
        out.append(line)
    return out


def load_inputs():
    stars = rows("desi_queue.csv")
    gs = rows("gaia_queue.csv")
    if Path("sdss_queue.csv").exists():
        stars += rows("sdss_queue.csv")
        gs += rows("sdss/gaia_tail.csv")
    if Path("moderate/queue.csv").exists():
        import glob

        focus = {"4435527217525195264", "4550108286326235520"}
        stars += [s for s in rows("moderate/queue.csv") if s["source_id"] in focus]
        gs += [g for path in glob.glob("moderate/gaia/rows_*.csv") for g in rows(path)]
    dwarfs = rows("literature/followup/dwarf_distances.csv")
    sh = {r["source_id"]: r for r in rows("literature/followup/starhorse.csv")}
    xp = {r["source_id"]: r for r in rows("desi/los/gspphot.csv")}
    mines = {r["source_id"]: r for r in rows("sdss/audit/minesweeper.csv")}
    for s in stars:
        sid = s["source_id"]
        if sid in sh:
            s.update(
                zip(
                    ("starhorse_lo", "starhorse_med", "starhorse_hi"),
                    [sh[sid]["dist" + q] for q in ("16", "50", "84")],
                )
            )
        if sid in xp:
            s.update(
                zip(
                    ("gspphot_med", "gspphot_lo", "gspphot_hi"),
                    [
                        number(xp[sid][k]) / 1000
                        for k in (
                            "distance_gspphot",
                            "distance_gspphot_lower",
                            "distance_gspphot_upper",
                        )
                    ],
                )
            )
        if sid in mines:
            m = mines[sid]
            d = number(m["Dist"])
            s.update(
                minesweeper_med=d,
                minesweeper_lo=d - number(m["Dist_lerr"]),
                minesweeper_hi=d + number(m["Dist_uerr"]),
            )
        # Prefer direct SDSS ugriz when both it and PS1 transformations are available.
        phot = [r for r in dwarfs if r["source_id"] == sid]
        if phot:
            fit = next((r for r in phot if r["survey"] == "SDSS"), phot[0])
            d = number(fit["distance_single_kpc"])
            s["dwarf_d"] = d
            s["dwarf_e"] = d * (0.25 if number(fit["g_i_dered"]) < 0.4 else 0.20)
    return stars, {g["source_id"]: g for g in gs}


if __name__ == "__main__":
    stars, gs = load_inputs()
    allrows = []
    for s in stars:
        allrows.extend(summarize(s, gs[s["source_id"]]))
    save("inference.csv", allrows)
