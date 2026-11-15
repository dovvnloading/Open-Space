"""Reproduce target measurements, Fourier fits, bootstrap and cadence controls."""
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from scipy.optimize import minimize_scalar
from screen import PHASE, URL, fetch, model, photo

NAME = "OGLE-LMC-RRLYR-09022"
SOURCE = "4652104345088779264"
TAP = "https://gaia.aip.de/tap/sync?"
RA, DEC = 15 * (5 + 15/60 + 25.99/3600), -(70 + 29/60 + 51.4/3600)
SEED = 9022
PERIODS = np.array([.360343, .485272])
RANGES = [(2000, 4600), (5200, 5750), (5750, 6120), (6120, 6485),
          (6485, 6850), (6850, 7215), (7215, 7580)]


def gaia(query, cache):
    url = TAP + urlencode(dict(REQUEST="doQuery", LANG="ADQL", QUERY=query))
    raw = fetch(url, cache)
    root = ET.fromstring(raw)
    fields = [x.get("name") for x in root.findall(".//{*}FIELD")]
    rows = [dict(zip(fields, [v.text for v in row])) for row in root.findall(".//{*}TR")]
    return rows, url


def fit(data, p, harmonics=2):
    t, y, e = data.T
    X = model(t, p, harmonics)
    w = 1 / np.maximum(.015, e)
    design = X * w[:, None]
    coeff = np.linalg.lstsq(design, y * w, rcond=None)[0]
    prediction = X @ coeff
    return dict(chi2=float(np.sum(((y - prediction)*w)**2)), prediction=prediction,
                amp=np.hypot(coeff[1::2*harmonics], coeff[2::2*harmonics]))


def refine(data, centers, halfwidth=.0012):
    """Resolve narrow local aliases by searching the whole bracket before refinement."""
    p = np.array(centers, dtype=float)
    low, high = p - halfwidth, p + halfwidth
    for _ in range(2):
        for index in range(len(p)):
            def objective(period):
                test = p.copy()
                test[index] = period
                return fit(data, test)["chi2"]
            grid = np.linspace(low[index], high[index], 101)
            scores = [objective(v) for v in grid]
            loc = int(np.argmin(scores))
            left, right = grid[max(loc-1, 0)], grid[min(loc+1, len(grid)-1)]
            p[index] = minimize_scalar(objective, bounds=(left, right), method="bounded",
                                       options=dict(xatol=1e-10)).x
    return p


def bootstrap(data, periods, count, rng):
    # Resample whole observing nights. Periods are conditional on this window's fit.
    groups = np.unique(np.floor(data[:, 0]), return_inverse=True)[1]
    blocks = [data[groups == i] for i in np.unique(groups)]
    trial = np.array([fit(np.vstack([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))]),
                          periods)["amp"] for _ in range(count)])
    return np.percentile(trial, [2.5, 97.5], axis=0).T.tolist()


def summarize(data, rng, count, fixed=None):
    p = refine(data, PERIODS) if fixed is None else np.array(fixed)
    sol = fit(data, p)
    return dict(n=len(data), first=data[:, 0].min(), last=data[:, 0].max(), periods=p.tolist(),
                semiamp=sol["amp"].tolist(), night_bootstrap_95=bootstrap(data, p, count, rng),
                median_err=np.median(data[:, 2]))


def residual_search(data, primary, grid):
    """Full-frequency scan after subtraction of the independently fitted dominant mode."""
    from scipy.signal import lombscargle
    base = fit(data, [primary])
    t, y, _ = data.T
    residual = y - base["prediction"]
    power = lombscargle(t - t.mean(), residual - residual.mean(), 2*np.pi*grid, normalize=True)
    k = int(np.argmax(power))
    return dict(best_frequency=grid[k], best_period=1/grid[k], peak_power=power[k])


def main(cache, count):
    rng = np.random.default_rng(SEED)
    urls = [URL.format(PHASE[p], "lmc", f"phot/{band}/{NAME}.dat")
            for p in (3, 4) for band in ("I", "V")]
    data = {f"{p}{b}": photo(p, NAME, b, cache) for p in (3, 4) for b in ("I", "V")}
    cone = (f"SELECT source_id,ra,dec,phot_g_mean_mag,phot_variable_flag,has_epoch_photometry "
            f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
            f"CIRCLE('ICRS',{RA},{DEC},0.0017))")
    queries = [cone, f"SELECT * FROM gaiadr3.vari_rrlyrae WHERE source_id={SOURCE}",
               f"SELECT * FROM gaiadr3.epoch_photometry WHERE source_id={SOURCE}"]
    with ThreadPoolExecutor(3) as pool:
        replies = list(pool.map(lambda q: gaia(q, cache), queries))
    urls.extend(u for rows, u in replies)
    neighbors, meta, epochs = [rows for rows, u in replies]
    field = epochs[0]
    def col(name):
        return np.array(field[name].split(), dtype=float)
    ok = np.array(field["variability_flag_g_reject"].split()) == "False"
    g = np.c_[col("g_transit_time")+5197.5, col("g_transit_mag"),
              (2.5 / np.log(10)) / col("g_transit_flux_over_error")]
    g = g[ok & np.isfinite(g).all(axis=1)]
    g = g[g[:, 0].argsort()]
    for row in neighbors:
        row["arcsec"] = np.hypot((float(row["ra"])-RA)*np.cos(np.deg2rad(DEC)),
                                 float(row["dec"])-DEC)*3600
    neighbors.sort(key=lambda r:r["arcsec"])

    windows = {}
    for tag, subset in data.items():
        for start, stop in RANGES:
            chosen = subset[(subset[:,0]>=start) & (subset[:,0]<stop)]
            if len(chosen) >= 30:
                windows[f"{tag}_{start}_{stop}"] = summarize(chosen, rng, count)
    early, late = data["3I"], data["4I"]
    tests = {}
    p_early = np.array(windows["3I_2000_4600"]["periods"])
    p_late = np.array(windows["4I_6850_7215"]["periods"])
    early_fit = fit(early, p_early)
    late_fit = fit(late[(late[:,0]>=6850)&(late[:,0]<7215)], p_late)
    # Does the early window suppress a late-sized fundamental? Keep residuals,
    # observation times and errors, remove both fitted modes, and re-inject at random phases.
    residual = early[:,1] - early_fit["prediction"]
    samples = []
    for phase in rng.uniform(0, 2*np.pi, count):
        injected = early.copy()
        injected[:,1] = (fit(early,[p_early[0]])["prediction"] + residual
                         + late_fit["amp"][1]*np.sin(2*np.pi*early[:,0]/p_late[1]+phase))
        samples.append(fit(injected, [p_early[0],p_late[1]])["amp"][1])
    tests["early_injected_amplitude"] = late_fit["amp"][1]
    tests["early_injection_recovery_p025_p50_p975"] = np.percentile(samples,[2.5,50,97.5]).tolist()
    # Free search of the early and late windows in a broad frequency interval independent of
    # which catalog mode is expected. A one-mode fit does not use the other mode's period.
    grid = np.linspace(1.2, 4.5, 13000)
    for key, d, subtract in [
            ("earlyI_after_overtone", early, p_early[0]),
            ("lateI_after_fundamental", late[(late[:,0]>=6850)&(late[:,0]<7215)], p_late[1]),
            ("middleI_after_overtone", late[(late[:,0]>=5750)&(late[:,0]<6120)],
             windows["4I_5750_6120"]["periods"][0])]:
        tests[key] = residual_search(d, subtract, grid)

    # Independent Gaia holdout at the published Gaia fundamental period. Only one
    # fundamental is fitted; the ability to fit a weak secondary in 39 visits is limited.
    pf = float(meta[0]["pf"])
    gf = fit(g,[pf])
    t,y,e = g.T
    loo, const = [], []
    for i in range(len(g)):
        mask = np.arange(len(t)) != i
        X = model(t, [pf], 2)
        w = 1/np.maximum(e,.015)
        beta = np.linalg.lstsq(X[mask]*w[mask,None],y[mask]*w[mask],rcond=None)[0]
        loo.append(y[i]-X[i]@beta)
        const.append(y[i]-np.average(y[mask],weights=w[mask]**2))
    # Same astrometric source and period as OGLE, with its own independent detector.
    gaia_result = dict(source=SOURCE, n=len(g), first=t.min(), last=t.max(),
                       catalog_type=meta[0]["best_classification"], catalog_period=pf,
                       catalog_period_error=float(meta[0]["pf_error"]),
                       semiamp=gf["amp"].item(), loo_rms=float(np.sqrt(np.mean(np.square(loo)))),
                       constant_loo_rms=float(np.sqrt(np.mean(np.square(const)))),
                       cone=neighbors)
    # Use scientific-row hashes for TAP: dynamic job metadata in response headers may vary.
    row_hashes = [hashlib.sha256(json.dumps(r,sort_keys=True).encode()).hexdigest()
                  for r in [sorted([dict((k,v) for k,v in r.items() if k!="arcsec") for r in neighbors],
                                   key=lambda r:r["source_id"]), meta, epochs]]
    sources = {u: hashlib.sha256(fetch(u,cache)).hexdigest() for u in urls[:4]}
    return dict(target=NAME, windows=windows, checks=tests, gaia=gaia_result,
                ogle_sha256=sources, tap_queries=urls[4:], tap_scientific_rows_sha256=row_hashes,
                bootstrap_trials=count, seed=SEED)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("raw-cache"))
    parser.add_argument("--out", type=Path, default=Path("measurements.json"))
    parser.add_argument("--trials", type=int, default=1000)
    args = parser.parse_args()
    args.out.write_text(json.dumps(main(args.cache,args.trials), indent=2, allow_nan=False) + "\n")
