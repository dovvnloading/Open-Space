"""Lens class and mass for OGLE-2026-BLG-0058 from posterior samples.

Written independently of the result's galaxy.py. Choices that differ from it:
- an old population's present-day mass function with stellar remnants: Kroupa
  IMF from 0.08 to 120 Msun, main-sequence stars up to a 1.0 Msun turn-off,
  white dwarfs from 1-8 Msun (Kalirai et al. 2008 relation), neutron stars at
  1.35 +- 0.07 Msun from 8-20 Msun, black holes log-normal about 9 Msun from
  20-120 Msun ("pdmf"); or a power law M^slope for comparison ("powerlaw");
- natal kicks: 265 km/s per axis for neutron stars, BH_KICK for black holes;
- a bar rotating as a solid body (50 km/s/kpc, capped at 120 km/s), angle 27 deg;
- source distances weighted by rho(DS) DS^2.
The event-rate weight per unit ln M and DS is g(M) M * M DL^4 mu^4 * sum over
lens and source populations of rho_L f(mu) rho_S DS^2, which follows from
Batista et al. (2011) eq. 18 after changing variables from mu to M.
Main-sequence lenses brighter than the blend of their own posterior sample are
rejected; white dwarfs, neutron stars and black holes are dark.

Usage: python -I check/galaxy_check.py <chains.npz> <out.json> [quick]
"""
import json
import sys

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric_posvel
from astropy.time import Time

KAPPA, R0, V0, KV = 8.144, 8.2, 233.0, 4.74047
TPAR = 2461213.0
TARGET = SkyCoord("17:52:38.09 -31:47:36.1", unit=(u.hourangle, u.deg))
GAL = TARGET.galactic
L, B = GAL.l.rad, GAL.b.rad
PA_L = TARGET.position_angle(SkyCoord(l=GAL.l + 1 * u.arcmin, b=GAL.b, frame="galactic")).rad
PA_B = TARGET.position_angle(SkyCoord(l=GAL.l, b=GAL.b + 1 * u.arcmin, frame="galactic")).rad
_ra, _dec = TARGET.ra.rad, TARGET.dec.rad
EAST = np.array([-np.sin(_ra), np.cos(_ra), 0.0])
NORTH = np.array([-np.sin(_dec) * np.cos(_ra), -np.sin(_dec) * np.sin(_ra), np.cos(_dec)])
_v = get_body_barycentric_posvel("earth", Time(TPAR, format="jd", scale="tdb"))[1].xyz.to(u.au / u.yr).value
V_EARTH = (float(NORTH @ _v), float(EAST @ _v))     # au/yr, projected north and east

SIG = {"disc": (35.0, 25.0), "bar": (115.0, 100.0)}   # km/s along l and b
MS_MASS = np.array([0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00])
MS_MI = np.array([13.5, 12.4, 10.3, 9.3, 8.5, 7.7, 6.9, 6.2, 5.6, 5.0, 4.3])


def xyz(d):
    return d * np.cos(B) * np.cos(L) - R0, d * np.cos(B) * np.sin(L), d * np.sin(B)


def density(d):
    x, y, z = xyz(d)
    disc = 0.05 * np.exp(-(np.hypot(x, y) - R0) / 2.6 - np.abs(z) / 0.3)
    a = np.radians(27.0)
    xb, yb = -x * np.cos(a) + y * np.sin(a), x * np.sin(a) + y * np.cos(a)
    rs = (((xb / 1.58) ** 2 + (yb / 0.62) ** 2) ** 2 + (z / 0.43) ** 4) ** 0.25
    return {"disc": disc, "bar": 1.23 * np.exp(-0.5 * rs ** 2)}


def mean_mu(d, pop):
    """Mean proper motion (mas/yr, along l and b) of a population at distance d, relative to the Sun."""
    x, y, _ = xyz(d)
    if pop == "disc":
        r = np.maximum(np.hypot(x, y), 1e-3)
        vx, vy = (V0 - 10.0) * y / r, -(V0 - 10.0) * x / r
    else:
        vx, vy = 50.0 * y, -50.0 * x
        cap = np.minimum(1.0, 120.0 / np.maximum(np.hypot(vx, vy), 1e-9))
        vx, vy = vx * cap, vy * cap
    vx, vy, vz = vx - 11.1, vy - (V0 + 12.24), -7.25
    vl = -vx * np.sin(L) + vy * np.cos(L)
    vb = -vx * np.sin(B) * np.cos(L) - vy * np.sin(B) * np.sin(L) + vz * np.cos(B)
    return vl / (KV * d), vb / (KV * d)


def kroupa(m):
    return np.where(m < 0.5, m ** -1.3, 0.5 * m ** -2.3) * ((m >= 0.08) & (m <= 120.0))


def mass_function(m, mode):
    if mode.startswith("powerlaw"):
        slope = float(mode.split(":")[1])
        return {"ANY": m ** slope}
    mi = (m - 0.394) / 0.109
    n_ns = 0.5 * (8.0 ** -1.3 - 20.0 ** -1.3) / 1.3
    n_bh = 0.5 * (20.0 ** -1.3 - 120.0 ** -1.3) / 1.3
    return {"MS": np.where(m <= 1.0, kroupa(m), 0.0),
            "WD": np.where((mi > 1.0) & (mi <= 8.0), kroupa(mi) / 0.109, 0.0),
            "NS": n_ns * np.exp(-0.5 * ((m - 1.35) / 0.07) ** 2) / (0.07 * np.sqrt(2 * np.pi)),
            "BH": n_bh * np.exp(-0.5 * (np.log(m / 9.0) / 0.35) ** 2) / (m * 0.35 * np.sqrt(2 * np.pi))}


def analyse(chain, fb, mode="pdmf", a_i=2.0, bh_kick=50.0, n_keep=1500, seed=3, light_cut=True):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(chain), n_keep, replace=False)
    te, pn, pe, fbs = chain[idx, 2] / 365.25, chain[idx, 3], chain[idx, 4], fb[idx]
    pie = np.hypot(pn, pe)
    lnm = np.linspace(np.log(0.08), np.log(80.0), 800)
    m = np.exp(lnm)
    ds = np.linspace(4.5, 11.5, 29)
    g = mass_function(m, mode)
    kick = {"MS": 0.0, "WD": 0.0, "NS": 265.0, "BH": bh_kick, "ANY": 0.0}
    rho_s = density(ds)
    src_w = {p: rho_s[p] * ds ** 2 for p in rho_s}
    mu_s = {p: mean_mu(ds, p) for p in rho_s}
    w = {c: np.zeros((n_keep, len(m))) for c in g}
    for lo in range(0, n_keep, 25):
        sl = slice(lo, lo + 25)
        th = KAPPA * m[None, :, None] * pie[sl, None, None]
        pirel = th * pie[sl, None, None]
        dl = 1.0 / (pirel + 1.0 / ds[None, None, :])
        mug = th / te[sl, None, None]
        mun = mug * (pn[sl] / pie[sl])[:, None, None] + V_EARTH[0] * pirel
        mue = mug * (pe[sl] / pie[sl])[:, None, None] + V_EARTH[1] * pirel
        mul = mun * np.cos(PA_L) + mue * np.sin(PA_L)
        mub = mun * np.cos(PA_B) + mue * np.sin(PA_B)
        rho_l = density(dl)
        base = m[None, :, None] * dl ** 4 * mug ** 4
        if light_cut:
            share = (1 - np.exp(-dl * abs(np.sin(B)) / 0.1)) / (1 - np.exp(-ds[None, None, :] * abs(np.sin(B)) / 0.1))
            ms_i = np.interp(np.log10(np.clip(m, 0.08, 1.0)), np.log10(MS_MASS), MS_MI)[None, :, None] \
                + 5 * np.log10(100 * dl) + a_i * share
            ms_ok = 10 ** (-0.4 * (ms_i - 18.0)) <= np.maximum(fbs[sl], 0.0)[:, None, None]
        for c in g:
            f = 0.0
            for pl in rho_l:
                ml_l, ml_b = mean_mu(dl, pl)
                for ps in src_w:
                    vl = ((SIG[pl][0] ** 2 + kick[c] ** 2) / dl ** 2 + SIG[ps][0] ** 2 / ds ** 2) / KV ** 2
                    vb = ((SIG[pl][1] ** 2 + kick[c] ** 2) / dl ** 2 + SIG[ps][1] ** 2 / ds ** 2) / KV ** 2
                    dmu_l, dmu_b = mul - (ml_l - mu_s[ps][0]), mub - (ml_b - mu_s[ps][1])
                    f = f + rho_l[pl] * src_w[ps] * np.exp(-0.5 * (dmu_l ** 2 / vl + dmu_b ** 2 / vb)) \
                        / (2 * np.pi * np.sqrt(vl * vb))
            ww = (g[c] * m)[None, :, None] * base * f
            if light_cut and c == "MS":
                ww = np.where(ms_ok, ww, 0.0)
            w[c][sl] = ww.sum(axis=2)
    total = sum(w[c] for c in w)
    z = float(total.sum() / n_keep)
    per_m = total.sum(axis=0)
    cdf = np.cumsum(per_m) / per_m.sum()
    q = np.interp([0.16, 0.5, 0.84], cdf, m)
    res = {"Z": z, "mass_16_50_84": np.round(q, 2).tolist(),
           "P(M>3)": round(float(per_m[m > 3].sum() / per_m.sum()), 3),
           "P(class)": {c: round(float(w[c].sum() / total.sum()), 4) for c in w}}
    if "BH" in w:
        bh = w["BH"].sum(axis=0)
        cb = np.cumsum(bh) / bh.sum()
        res["BH_mass_16_50_84"] = np.round(np.interp([0.16, 0.5, 0.84], cb, m), 1).tolist()
    return res


def main():
    data = np.load(sys.argv[1])
    out = {"l_deg": float(np.degrees(L)), "b_deg": float(np.degrees(B)),
           "PA_increasing_l_deg": float(np.degrees(PA_L)), "PA_north_galactic_deg": float(np.degrees(PA_B)),
           "earth_velocity_ne_au_per_yr": V_EARTH, "runs": {}}
    runs = [("powerlaw:-1.0", 2.0, 50.0, False), ("powerlaw:-2.35", 2.0, 50.0, False),
            ("pdmf", 1.5, 50.0, True), ("pdmf", 2.5, 50.0, True),
            ("pdmf", 2.0, 0.0, True), ("pdmf", 2.0, 100.0, True), ("pdmf", 2.0, 50.0, False)]
    if len(sys.argv) > 3 and sys.argv[3] == "quick":
        runs = [("powerlaw:-2.35", 2.0, 50.0, False), ("pdmf", 2.0, 50.0, True)]
    for mode, a_i, kick, cut in runs:
        key = f"{mode} A_I={a_i} BHkick={kick} lightcut={cut}"
        out["runs"][key] = {}
        for sol in ("plus", "minus"):
            out["runs"][key][sol] = analyse(data[sol], data[sol + "_fb"], mode, a_i, kick, light_cut=cut)
        zp, zm = out["runs"][key]["plus"]["Z"], out["runs"][key]["minus"]["Z"]
        out["runs"][key]["prior_odds_plus_over_minus"] = zp / zm
        print(key, json.dumps(out["runs"][key]))
    with open(sys.argv[2], "w", newline="\n") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in ("l_deg", "b_deg", "PA_increasing_l_deg", "PA_north_galactic_deg",
                                          "earth_velocity_ne_au_per_yr")}))


if __name__ == "__main__":
    main()
