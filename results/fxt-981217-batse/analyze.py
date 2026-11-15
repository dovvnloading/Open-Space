"""Recompute the WFC–BATSE crossmatch and fixed-window photometry from the extracts."""
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DAY = 86400.0
END_OF_BATSE = 51699.0  # All retained WFC targets precede CGRO shutdown, 2000-06-04.


def rows(name):
    with (ROOT / "inputs" / name).open() as handle:
        return list(csv.DictReader(handle))


def vector(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.stack([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)], axis=-1)


def crossmatch():
    wfc = rows("wfc_catalog.csv")
    stern = rows("stern_catalog.csv")
    subset = [r for r in wfc if r["class"] == "GRB?" and float(r["mjd"]) < END_OF_BATSE]
    times = np.array([float(r["mjd"]) for r in stern])
    coords = vector([float(r["ra"]) for r in stern], [float(r["dec"]) for r in stern])
    radii = np.array([float(r["r68_deg"]) for r in stern])
    rng = np.random.default_rng(981217)
    limits = []
    for _ in range(50000):
        # Reflect positive and negative shifts. Uniform in elapsed days, not burst density.
        limits.append(rng.uniform(1, 180) * rng.choice([-1, 1]))
    shifts = np.array([0.0] + limits)
    families = [[r for r in wfc if float(r["mjd"]) < END_OF_BATSE], subset]
    result = {"wfc_rows": len(wfc), "ambiguous_wfc_rows": sum(r["class"] == "GRB?" for r in wfc),
              "stern_rows": len(stern), "trials": len(limits)}
    for label, family in zip(["all_classes", "ambiguous"], families):
        # Retain original discovery screen and offer a stringent one-radius sensitivity test.
        for gate in ["wide", "r68"]:
            at_least_one = np.zeros(shifts.size, bool)
            matches = []
            for w in family:
                time = float(w["mjd"])
                pos = vector(float(w["ra"]), float(w["dec"]))
                sep = np.rad2deg(np.arccos(np.clip(coords @ pos, -1, 1)))
                window = np.maximum(20, 2*radii) if gate == "wide" else radii
                ix = np.nonzero(sep <= window)[0]
                ix = ix[np.argsort(times[ix])]
                tt = times[ix]
                width = float(w["t90_s"])
                # Same screening interval for data and time shifts.
                low = time + shifts - 120/DAY
                high = time + shifts + (width + 120)/DAY
                a, b = np.searchsorted(tt, low), np.searchsorted(tt, high, side="right")
                at_least_one |= b > a
                for k in ix[a[0]:b[0]]:
                    matches.append({"wfc": w["id"], "stern": stern[k]["id"],
                                    "delta_s": round((times[k]-time)*DAY, 4),
                                    "separation_deg": round(float(sep[k]), 5),
                                    "radius_deg": radii[k], "wfc_batse_id": w["batse_id"]})
            count = int(at_least_one[1:].sum())
            result[f"{label}_{gate}"] = {"targets": len(family), "zero_lag_matches": matches,
                                          "shifts_with_any": count,
                                          "fraction_plus_one": (1+count)/shifts.size}
    w = next(r for r in wfc if r["id"] == "981217")
    gusbad = rows("gusbad_excerpt.csv")[0]
    sep_g = math.degrees(math.acos(np.clip(vector(float(w["ra"]), float(w["dec"])) @
                                           vector(float(gusbad["ra"]), float(gusbad["dec"])), -1, 1)))
    radius_w = math.radians(float(w["r68_arcmin"])/60)
    radius_b = math.radians(next(radii[k] for k, r in enumerate(stern) if r["id"] == "STK01 11164a"))
    result["gusbad_crosscheck"] = {"id":gusbad["name"], "delta_s":round((float(gusbad["time"])-float(w["mjd"]))*DAY,4),
                                  "separation_deg":round(sep_g,5),"radius_deg":float(gusbad["error_radius"])}
    result["reported_circle_area_ratio"] = round((math.sin(radius_b/2)/math.sin(radius_w/2))**2, 1)
    return result


def measure():
    raw = np.genfromtxt(ROOT / "inputs/batse_lc.csv", delimiter=",", names=True)
    wfc = np.genfromtxt(ROOT / "inputs/wfc_lc.csv", delimiter=",", names=True)
    t0 = 51164.079960
    t = (raw["mjd"]-t0)*DAY
    t90 = 111.
    channels = np.stack([sum(raw[f"d{det}_ch{ch}"] for det in range(8)) for ch in range(1, 7)], axis=1)
    channel23 = channels[:, 1:3].sum(axis=1)
    y = np.column_stack((channels, channel23))
    source = (t >= 0) & (t < t90)
    results = {}
    for degree, edge in [(1,210), (2,210), (1,150), (2,150)]:
        off = ((t >= -edge) & (t < -30)) | ((t >= t90+30) & (t < t90+edge))
        X = np.polynomial.polynomial.polyvander((t-t90/2)/100, degree)
        coef = np.linalg.lstsq(X[off], y[off], rcond=None)[0]
        resid = y - X@coef
        sum_x = X[source].sum(axis=0)
        factor = source.sum() + sum_x@np.linalg.solve(X[off].T@X[off],sum_x)
        error = np.sqrt(factor * (resid[off]**2).sum(axis=0)/(off.sum()-degree-1))
        results[f"degree{degree}_outer{edge}"] = {"source_bins":int(source.sum()), "background_bins":int(off.sum()),
            "net_counts": np.round(resid[source].sum(axis=0), 1).tolist(),
            "residual_se_counts": np.round(error, 1).tolist()}
        if degree == 1 and edge == 210:
            primary = resid[:, -1]
    bt = (wfc["mjd_begin"]-t0)*DAY
    et = (wfc["mjd_end"]-t0)*DAY
    mid = (bt+et)/2
    lc = []
    for lo, hi, flux, err in zip(bt, et, wfc["flux"], wfc["flux_error"]):
        use = (t >= lo) & (t < hi)
        rate = primary[use].mean()/1.024 if use.sum()>=8 else None
        lc.append({"begin_s":round(lo,3), "end_s":round(hi,3), "wfc_flux":float(flux), "wfc_error":float(err),
                   "batse_excess_rate":round(float(rate),2) if rate is not None else None,"batse_bins":int(use.sum())})
    ig = np.nonzero((mid >= 0) & (mid <= t90))[0]
    peak_w = ig[np.argmax(wfc["flux"][ig])]
    peak_b = ig[np.argmax([lc[i]["batse_excess_rate"] for i in ig])]
    # LAD 7/3 come from GUSBAD's own catalog. All eight LADs enter the primary measurement.
    x = np.column_stack((np.ones(t.size), (t-t90/2)/100))
    bg = ((t >= -210) & (t < -30)) | ((t >= t90+30) & (t < t90+210))
    det_y = np.stack([raw[f"d{det}_ch2"]+raw[f"d{det}_ch3"] for det in range(8)], axis=1)
    detector_residuals = det_y - x@np.linalg.lstsq(x[bg],det_y[bg],rcond=None)[0]
    results["all_lad_channel_23_net_counts"] = np.round(detector_residuals[source].sum(axis=0),1).tolist()
    results["peak_bins"] = {"wfc":lc[peak_w], "batse":lc[peak_b]}
    # No lag fit. The WFC peak is broad; use equal early/late thirds for a less bin-sensitive diagnostic.
    for name, lo, hi in [("first_third",0,37),("last_third",74,111)]:
        g = (mid>=lo)&(mid<hi)
        bins = [lc[i] for i in np.nonzero(g)[0]]
        results[name] = {"wfc_bins": len(bins), "wfc_mean_flux":round(float(wfc["flux"][g].mean()),4),
                         "wfc_se": round(float(np.sqrt(np.sum(wfc["flux_error"][g]**2))/g.sum()),4),
                         "batse_mean_rate":round(float(np.mean([b["batse_excess_rate"] for b in bins])),2)}
    results["flags_on_source"] = [[list(map(int, fb)), int(n)] for fb,n in
        zip(*np.unique(np.stack([raw[f"flag{i}"][source] for i in range(1,5)], axis=1), axis=0, return_counts=True))]
    return results, lc


if __name__ == "__main__":
    result = crossmatch()
    result["photometry"], lightcurve = measure()
    (ROOT / "measurements.json").write_text(json.dumps(result, indent=2)+"\n")
    with (ROOT / "joint_lightcurve.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(lightcurve[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(lightcurve)
