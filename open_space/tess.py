"""Search a MAST two-minute TESS sector and predict the same ephemeris in another."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

from .tess_archive import (DEFAULT_PLAN, SEED, atomic_bytes, cached_product,
                           fetch_plan, make_plan)


MIN_PERIOD = 0.75
MAX_PERIOD = 8.0
DURATIONS = (1/24, 2/24, 4/24)
BIN_DAYS = 10 / 1440
ALGORITHM = "spoc-two-sector-bls-v1"


def load_product(product, cache, tic):
    import numpy as np
    from astropy.io import fits

    path, record = cached_product(product, cache)
    with fits.open(path, memmap=False) as hdus:
        head = hdus[0].header
        table = hdus["LIGHTCURVE"]
        info = table.header
        if (int(head["TICID"]) != tic or info["TIMESYS"] != "TDB" or info["TIMEUNIT"] != "d"
                or info["BJDREFI"] + info["BJDREFF"] != 2457000):
            raise ValueError("Target or time convention differs from planned two-minute SPOC product")
        time = np.asarray(table.data["TIME"], dtype=float)
        flux = np.asarray(table.data["PDCSAP_FLUX"], dtype=float)
        err = np.asarray(table.data["PDCSAP_FLUX_ERR"], dtype=float)
        good = (table.data["QUALITY"] == 0) & np.isfinite(time) & np.isfinite(flux) & np.isfinite(err) & (err > 0)
        if not good.any():
            raise ValueError("No valid PDCSAP cadences")
        t, f, e = time[good], flux[good], err[good]
        if not np.all(np.diff(t) > 0):
            raise ValueError("Expected increasing time samples")
        record.update(tic_id=tic, ra_deg=float(head["RA_OBJ"]), dec_deg=float(head["DEC_OBJ"]),
                      camera=int(head["CAMERA"]), ccd=int(head["CCD"]), start_btjd=float(t[0]),
                      stop_btjd=float(t[-1]), valid_cadences=int(good.sum()), total_cadences=len(time),
                      timescale="BJD_TDB - 2457000")
        for key in ("CROWDSAP", "FLFRCSAP"):
            val = info.get(key)
            record[key.lower()] = float(val) if val is not None and np.isfinite(val) else None
    return (t, f, e), record


def flatten(curve):
    """Bin once, then model a continuously varying baseline per uninterrupted run."""
    import numpy as np

    t, f, err = curve
    labels = np.floor((t-t[0])/BIN_DAYS).astype(int)
    starts = np.r_[0, np.flatnonzero(np.diff(labels)) + 1]
    n = np.diff(np.r_[starts, len(t)])
    w = 1/err**2
    sw = np.add.reduceat(w, starts)
    keep = n >= 3
    t = (np.add.reduceat(w*t, starts) / sw)[keep]
    f = (np.add.reduceat(w*f, starts) / sw)[keep]
    err = (1/np.sqrt(sw))[keep]
    ts, fs, es = [], [], []
    for ix in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > 0.35) + 1):
        if len(ix) < 60 or t[ix[-1]] - t[ix[0]] < 0.8:
            continue
        tx, fx, ex = t[ix], f[ix], err[ix]
        anchors = np.linspace(tx[0], tx[-1], max(2, int(np.ceil((tx[-1]-tx[0])/0.2))+1))
        meds = []
        for a in anchors:
            lo, hi = np.searchsorted(tx, [a-0.65, a+0.65])
            meds.append(np.median(fx[lo:hi]))
        baseline = np.interp(tx, anchors, meds)
        valid = np.isfinite(baseline) & (baseline > 0)
        tx, rel, ex = tx[valid], fx[valid]/baseline[valid]-1, ex[valid]/baseline[valid]
        if len(tx) < 2:
            continue
        noise = 1.4826*np.median(np.abs(np.diff(rel)-np.median(np.diff(rel))))/math.sqrt(2)
        ts.extend(tx)
        fs.extend(rel)
        es.extend(np.maximum(ex, noise))
    return tuple(np.asarray(a) for a in (ts, fs, es))


def enough(curve):
    return len(curve[0]) >= 200 and curve[0][-1] - curve[0][0] >= 8


def blind(curve):
    import numpy as np
    from astropy.timeseries import BoxLeastSquares

    t, f, err = curve
    span = t[-1]-t[0]
    dlog = min(DURATIONS)/(2.5*span)
    grid = np.geomspace(MIN_PERIOD, MAX_PERIOD, int(math.ceil(math.log(MAX_PERIOD/MIN_PERIOD)/dlog))+1)
    bls = BoxLeastSquares(t, f, dy=err)
    fit = bls.power(grid, DURATIONS, objective="snr", oversample=8)
    order = np.argsort(fit.power)[::-1]
    best = int(order[0])
    alternatives = []
    accepted = [float(grid[best])]
    for j in order[1:]:
        if all(abs(1/grid[j] - 1/p) * span > 0.5 for p in accepted):
            alternatives.append({"period_days": float(grid[j]), "bls_snr": float(fit.power[j])})
            accepted.append(float(grid[j]))
        if len(alternatives) == 3:
            break
    return {"period_days": float(fit.period[best]), "epoch_btjd": float(fit.transit_time[best]),
            "duration_days": float(fit.duration[best]), "bls_snr": float(fit.power[best]),
            "grid_trials": len(grid), "alternative_peaks": alternatives}


def _depth(values, errors):
    import numpy as np
    w = 1 / errors**2
    return float(np.sum(w*values)/np.sum(w)), float(1/math.sqrt(np.sum(w)))


def assess(curve, ephem):
    import numpy as np

    t, flux, err = curve
    p, t0, dur = (ephem[k] for k in ("period_days", "epoch_btjd", "duration_days"))
    events = []
    # Use each event's local baseline. A gap suppresses coverage instead of
    # interpolating the photometry. Uncertainty includes both local estimates.
    for epoch in range(math.ceil((t[0]-t0)/p-0.5), math.floor((t[-1]-t0)/p+0.5)+1):
        center = t0+epoch*p
        outer = min(2.5*dur, 0.42*p)
        lo, hi = np.searchsorted(t, [center-outer, center+outer])
        dx = t[lo:hi]-center
        inside = abs(dx) <= 0.5*dur
        outside = abs(dx) >= 0.78*dur
        if inside.sum() < max(3, math.ceil(0.50*dur/BIN_DAYS)) or outside.sum() < 5:
            continue
        if not ((dx[outside] < 0).any() and (dx[outside] > 0).any()):
            continue
        fm, fe = _depth(flux[lo:hi][inside], err[lo:hi][inside])
        bm, be = _depth(flux[lo:hi][outside], err[lo:hi][outside])
        events.append((epoch, bm-fm, math.hypot(fe,be)))
    if not events:
        return {"depth": None, "depth_err": None, "snr": None, "covered_events": 0,
                "supported_events": 0, "odd_even_sigma": None, "events": []}
    epoch, values, errors = np.array(events).T
    def aggregate(mask):
        if not mask.any():
            return None, None
        m, e = _depth(values[mask], errors[mask])
        # Inflate formal uncertainties for variable event depths.
        if mask.sum() > 1:
            e *= max(1, math.sqrt(float(np.sum(((values[mask]-m)/errors[mask])**2)/(mask.sum()-1))))
        return m, e
    all_events = np.ones(len(events), dtype=bool)
    depth, error = aggregate(all_events)
    a, ea = aggregate(epoch % 2 == 0)
    b, eb = aggregate(epoch % 2 == 1)
    odd_even = abs(a-b)/math.hypot(ea,eb) if a is not None and b is not None else None
    return {"depth": depth, "depth_err": error, "snr": depth/error,
            "covered_events": len(events), "supported_events": int(np.sum(values/errors >= 2)),
            "odd_even_sigma": odd_even,
            "events": [{"cycle": int(k), "depth": float(d), "depth_err": float(e)} for k,d,e in events]}


def search_pair(train, test):
    if not enough(train) or not enough(test):
        return {"disposition": "insufficient_coverage", "flags": ["insufficient_baseline_or_bins"]}
    ephem = blind(train)
    first, other = assess(train, ephem), assess(test, ephem)
    secondary_ephem = {**ephem, "epoch_btjd": ephem["epoch_btjd"] + ephem["period_days"]/2}
    secondary = [assess(curve, secondary_ephem) for curve in (train, test)]
    flags = []
    for label, primary, sec in zip(("discovery","validation"), (first,other), secondary):
        if (primary["odd_even_sigma"] or 0) > 4:
            flags.append(f"{label}_odd_even")
        if ((sec["snr"] or 0) >= 5 and (sec["depth"] or 0) > 0.2*(primary["depth"] or 0)):
            flags.append(f"{label}_secondary")
        if (primary["depth"] or 0) > 0.05:
            flags.append(f"{label}_deep")
    ratio = other["depth"]/first["depth"] if first["depth"] and other["depth"] is not None else None
    if ratio is not None and not 0.4 <= ratio <= 2.5:
        flags.append("sector_depth_inconsistent")
    two_sector = ((first["snr"] or 0) >= 7 and (other["snr"] or 0) >= 5 and
                  min(first["supported_events"], other["supported_events"]) >= 3)
    if two_sector:
        disposition = "periodic_dip_with_veto" if flags else "two_sector_periodic_dip"
    else:
        disposition = "below_two_sector_threshold"
    for sec in secondary:
        sec.pop("events")
    return {"ephemeris": ephem, "discovery": first, "validation": other,
            "secondary_by_sector": secondary, "validation_depth_ratio": ratio,
            "flags": flags, "disposition": disposition}


def run_target(target, cache_dir):
    """Persistent-worker adapter. Returns measured and rejected hypotheses alike."""
    products = target["products"]
    tic = int(target["tic_id"])
    if len(products) != 2 or products[0]["sector"] >= products[1]["sector"]:
        raise ValueError("Expected discovery and later validation sector")
    inputs = [load_product(p, cache_dir, tic) for p in products]
    flattened = [flatten(pair[0]) for pair in inputs]
    sources = [pair[1] for pair in inputs]
    result = search_pair(*flattened)
    if enough(flattened[0]) and enough(flattened[1]):
        inverted = [(t, -f, e) for t,f,e in flattened]
        inv = search_pair(*inverted)
        result["inverted_flux_control"] = {"disposition": inv["disposition"],
            "period_days": inv["ephemeris"]["period_days"], "flags": inv["flags"],
            "discovery_snr": inv["discovery"]["snr"], "validation_snr": inv["validation"]["snr"]}
    if target["role"] == "known_control":
        result["benchmark"] = {"name": "WASP-18 b", "catalog_period_days": 0.9414522,
                                "fractional_period_error": abs(result["ephemeris"]["period_days"]/0.9414522-1)
                                if "ephemeris" in result else None}
    return {"schema_version": 1, "track": "tess", "algorithm": ALGORITHM, "role": target["role"],
            "tic_id": tic, "ra_deg": sources[0]["ra_deg"], "dec_deg": sources[0]["dec_deg"],
            "binned_samples": [len(c[0]) for c in flattened], "sources": sources, **result}


def injection_experiment(target, cache):
    curves = [load_product(p, cache, target["tic_id"])[0] for p in target["products"]]
    trials = []
    for p in (1.7, 3.1, 5.3):
        for depth in (0.001, 0.005, 0.01):
            t0 = float(curves[0][0][0] + 0.63*p)
            modified = []
            for t, f, e in curves:
                phase = (t-t0+p/2)%p-p/2
                in_transit = abs(phase) <= 1/24
                modified.append(flatten((t, f*(1-depth*in_transit), e)))
            fit = search_pair(*modified)
            same_phase = None
            found = False
            if "ephemeris" in fit:
                meas = fit["ephemeris"]
                diff = (meas["epoch_btjd"]-t0+p/2)%p-p/2
                same_phase = float(diff)
                found = (abs(meas["period_days"]/p-1) < 0.01 and abs(diff) < 1/24 and
                         fit["disposition"] == "two_sector_periodic_dip")
            trials.append({"injected_period_days": p, "injected_depth": depth, "epoch_btjd": t0,
                           "duration_days": 2/24, "measured_period_days": fit.get("ephemeris",{}).get("period_days"),
                           "disposition": fit["disposition"], "epoch_offset_days": same_phase,
                           "recovered": bool(found)})
    return {"background_tic": target["tic_id"], "type": "injected_control", "trials": trials}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    plan = sub.add_parser("plan")
    plan.add_argument("--sectors", type=int, nargs=2, default=[1, 2])
    plan.add_argument("--offset", type=int, default=0)
    plan.add_argument("--limit", type=int, default=24)
    plan.add_argument("--seed", default=SEED)
    plan.add_argument("--include-control", action="store_true")
    plan.add_argument("--inventories", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    fetch = sub.add_parser("fetch")
    fetch.add_argument("--plan", type=Path, required=True)
    fetch.add_argument("--cache", type=Path, required=True)
    fetch.add_argument("--output", type=Path, required=True)
    fetch.add_argument("--workers", type=int, default=3)
    run = sub.add_parser("run")
    run.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    run.add_argument("--cache", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--injection-controls", action="store_true")
    run.add_argument("--workers", type=int, default=1, help="Parallel target processes (1–32)")
    args = parser.parse_args(argv)
    if args.action == "plan":
        result = make_plan(*args.sectors, args.limit, args.offset, args.seed, args.inventories, args.include_control)
    else:
        plan = json.loads(args.plan.read_bytes())
        if args.action == "fetch":
            result = fetch_plan(plan, args.cache, args.workers)
        else:
            if not 1 <= args.workers <= 32:
                raise ValueError("workers must be between 1 and 32")
            if args.workers == 1:
                targets = [run_target(t, args.cache) for t in plan["targets"]]
            else:
                with ProcessPoolExecutor(max_workers=args.workers) as pool:
                    targets = list(pool.map(partial(run_target, cache_dir=args.cache), plan["targets"]))
            result = {"schema_version": 1, "algorithm": ALGORITHM, "track":"tess",
                      "plan_sha256": hashlib.sha256(args.plan.read_bytes()).hexdigest(),
                      "targets": targets, "dispositions": dict(Counter(t["disposition"] for t in targets))}
            if args.injection_controls:
                first = next(t for t in plan["targets"] if t["role"] == "exploration")
                result["injection_controls"] = injection_experiment(first, args.cache)
    atomic_bytes(args.output, (json.dumps(result, indent=2, allow_nan=False)+"\n").encode())


if __name__ == "__main__":
    main()
