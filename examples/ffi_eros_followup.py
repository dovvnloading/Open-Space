"""Check image-period candidates with star-resolved EROS-2 and nearby OGLE light curves."""

import argparse
from collections import defaultdict
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode

import numpy as np
from astropy.timeseries import LombScargle

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from open_space.tess_archive import atomic_bytes, download

VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
LC_COLUMNS = "EROS,HJD,Rmag,e_Rmag,Bmag,e_Bmag"
# Pinned known neighbors of this particular image peak, independently identified by OGLE.
OGLE = {
    "OGLE-LMC-RRLYR-05589": ("rrlyr", 0.2476053, "HJD-2450000"),
    "OGLE-LMC-ELL-01519": ("ecl", 220.4366851, "HJD"),
}


def fetch(url, cache):
    key = hashlib.sha256(url.encode()).hexdigest()
    path = cache / (key + ".raw")
    raw = path.read_bytes() if path.exists() else download(url, 15_000_000)
    return raw, path, {"url": url, "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def vizier(table, ra, dec, limit, cache, all_columns=False):
    pars = {"-source": table}
    if all_columns:
        pars.update({"-c": f"{ra:.7f} {dec:.7f}", "-c.rs": "30", "-out.all": "1", "-out.max": limit})
    else:
        pars.update({"-out": LC_COLUMNS, "-c": f"{ra:.7f} {dec:.7f}", "-c.rs": "30", "-out.max": limit})
    raw, path, receipt = fetch(VIZIER + "?" + urlencode(pars), cache)
    lines = raw.decode().splitlines()
    if "#Name: " + table not in lines or any(l.startswith(("#ERROR", "#INFO\tError")) for l in lines):
        raise ValueError("Incomplete EROS query")
    data = [line for line in lines if line.strip() and not line.startswith("#")]
    if len(data) < 4 or not data[2].startswith("---"):
        raise ValueError("Missing EROS data")
    header = data[0].split("\t")
    rows = [dict(zip(header, (s.strip() for s in l.split("\t")), strict=True)) for l in data[3:]]
    if len(rows) >= limit:
        raise ValueError("Possibly truncated EROS response")
    if not path.exists():
        atomic_bytes(path, raw)
    return rows, {**receipt, "rows": len(rows)}


def clean(array):
    t, y, e = array.T
    ok = np.isfinite(t + y + e) & (e > 0) & (e <= 0.3)
    t, y, e = t[ok], y[ok], e[ok]
    if not len(t):
        return t, y, e
    order = np.argsort(t, kind="stable")
    t, y, e = t[order], y[order], np.sqrt(e[order] ** 2 + 0.01**2)
    # Independent per-season zeropoints suppress slow transparency and long-period leakage.
    seasons = np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > 75) + 1)
    for indices in seasons:
        y[indices] -= np.median(y[indices])
    mad = 1.4826 * np.median(abs(y - np.median(y)))
    ok = abs(y) <= 6 * np.maximum(mad, e)
    return t[ok], y[ok], e[ok]


def peak(curve, fmin, fmax):
    t, y, e = curve
    freq = np.arange(fmin, fmax, 1 / (5 * np.ptp(t)))
    power = LombScargle(t, y, e).power(freq, method="fast")
    i = int(np.argmax(power))
    return float(freq[i]), float(power[i])


def fit(curve, frequency, tref):
    t, y, e = curve
    phase = 2 * np.pi * (t - tref) * frequency
    design = np.column_stack((np.ones(len(t)), np.cos(phase), np.sin(phase)))
    beta = np.linalg.lstsq(design / e[:, None], y / e, rcond=None)[0]
    residual = y - design @ beta
    total = np.sum(((y - np.average(y, weights=1 / e**2)) / e) ** 2)
    return {
        "variance_fraction": float(1 - np.sum((residual / e) ** 2) / total),
        "semiamplitude_mag": float(np.hypot(*beta[1:])),
        "phase_cycles": float(math.atan2(beta[2], beta[1]) / (2 * math.pi)),
    }


def time_slice(curve, mask):
    return tuple(x[mask] for x in curve)


def survey_star(rows, band, window):
    curve = clean(
        np.array(
            [
                [float(row[k] or "nan") for k in ("HJD", band + "mag", "e_" + band + "mag")]
                for row in rows
            ]
        )
    )
    t, y, e = curve
    if len(t) < 120 or np.ptp(t) < 365:
        return {"usable": len(t), "status": "insufficient"}, None
    tmid = np.median(t)
    train, test = time_slice(curve, t <= tmid), time_slice(curve, t > tmid)
    frequency, train_power = peak(train, *window)
    locked = fit(test, frequency, t[0])
    train_fit = fit(train, frequency, t[0])
    full_frequency, full_power = peak(curve, *window)
    return {
        "status": "measured",
        "usable": len(t),
        "start_eros_hjd": float(t[0]),
        "end_eros_hjd": float(t[-1]),
        "median_error_mag_with_floor": float(np.median(e)),
        "all_peak_period_days": 1 / full_frequency,
        "all_peak_variance_fraction": full_power,
        "train_peak_period_days": 1 / frequency,
        "train_peak_variance_fraction": train_power,
        "validation_at_train_frequency": locked,
        "train_to_validation_phase_drift": abs(
            math.remainder(locked["phase_cycles"] - train_fit["phase_cycles"], 1)
        ),
        "train_count": len(train[0]),
        "validation_count": len(test[0]),
    }, curve


def daily_alias(curve, center):
    """Power near the first solar and sidereal sampling aliases in this ground survey."""
    out = []
    for frequency in (1.0, 1.00273791):
        for sign in (-1, 1):
            alias = abs(center + sign * frequency)
            f, power = peak(curve, max(0.025, alias - 1 / 27), alias + 1 / 27)
            out.append(
                {
                    "sampling_frequency_per_day": frequency,
                    "sign": sign,
                    "period_days": 1 / f,
                    "variance_fraction": power,
                }
            )
    return out


def recovery(cr, cb, window, phase, amplitude):
    """Use real noise, times and errors, with one injected phase per source in both bands."""
    center = np.mean(window)

    def injected(curve):
        t, y, e = curve
        return t, y + amplitude * np.sin(2 * np.pi * ((t - cr[0][0]) * center + phase)), e

    red, blue = injected(cr), injected(cb)
    median = np.median(red[0])
    tr = time_slice(red, red[0] <= median)
    vr = time_slice(red, red[0] > median)
    vb = time_slice(blue, blue[0] > median)
    frequency, power = peak(tr, *window)
    train, rv, bv = (fit(c, frequency, cr[0][0]) for c in (tr, vr, vb))
    drift = max(abs(math.remainder(v["phase_cycles"] - train["phase_cycles"], 1)) for v in (rv, bv))
    return (
        abs(frequency - center) <= 1 / np.ptp(tr[0])
        and min(power, rv["variance_fraction"], bv["variance_fraction"]) >= 0.05
        and drift <= 0.2
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--evidence", type=Path, default=Path("results/ffi-lmc-grid.json"))
    args = p.parse_args()
    source = next(
        s
        for f in json.loads(args.evidence.read_bytes())["fields"]
        for s in f["sources"]
        if s["source_id"] == "ICRS-77.05502-69.23710"
    )
    ra, dec = source["ra_deg"], source["dec_deg"]
    periods = [source[k]["period_days"] for k in ("discovery", "validation")]
    center = np.mean(1 / np.array(periods))
    window = (center - 1 / 27, center + 1 / 27)
    catalog, catalog_receipt = vizier("II/390/xcat", ra, dec, 500, args.cache, True)
    rows, lc_receipt = vizier("II/390/xtime", ra, dec, 100000, args.cache)
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["EROS"]].append(row)
    if set(grouped) != {c["EROS"] for c in catalog}:
        raise ValueError("Catalog roster differs from light curve roster")
    results = []
    for c in sorted(catalog, key=lambda c: float(c["_r"])):
        key = c["EROS"]
        if len(grouped[key]) != int(c["Nlc"]):
            raise ValueError(f"Incomplete source {key}")
        result = {
            "eros_id": key,
            "ra_j2000_deg": float(c["RAJ2000"]),
            "dec_j2000_deg": float(c["DEJ2000"]),
            "separation_arcsec": float(c["_r"]),
            "catalog_red_magnitude": float(c["Rmag"]),
            "exposure_rows": len(grouped[key]),
        }
        red, cr = survey_star(grouped[key], "R", window)
        blue, cb = survey_star(grouped[key], "B", window)
        result.update({"red": red, "blue": blue})
        if cr is not None and cb is not None:
            # Reserve all blue points at late epochs from red frequency selection.
            frequency = 1 / red["train_peak_period_days"]
            early_red = time_slice(cr, cr[0] <= np.median(cr[0]))
            late_blue = time_slice(cb, cb[0] > np.median(cr[0]))
            if len(late_blue[0]) >= 60:
                br = fit(late_blue, frequency, cr[0][0])
                rr = fit(early_red, frequency, cr[0][0])
                result["late_blue_at_early_red"] = br
                phase = int.from_bytes(hashlib.sha256(key.encode()).digest()[:4], "big") / 2**32
                result["injected_recovery"] = {
                    str(amplitude): bool(recovery(cr, cb, window, phase, amplitude))
                    for amplitude in (0.01, 0.03, 0.1)
                }
                result["red_to_blue_phase_drift"] = abs(
                    math.remainder(br["phase_cycles"] - rr["phase_cycles"], 1)
                )
            for band, curve in (("red", cr), ("blue", cb)):
                f, power = peak(curve, 0.1, 5.0)
                result[band]["wide_peak_period_days"] = 1 / f
                result[band]["wide_peak_variance_fraction"] = power
                result[band]["daily_aliases"] = daily_alias(curve, center)
        results.append(result)
    ogle = []
    for star, (group, known, scale) in OGLE.items():
        url = f"https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/{group}/phot/I/{star}.dat"
        raw, path, receipt = fetch(url, args.cache)
        arr = np.loadtxt(io.BytesIO(raw))
        if arr.shape[1] != 3:
            raise ValueError("Unexpected OGLE table")
        if not path.exists():
            atomic_bytes(path, raw)
        curve = clean(arr)
        narrow, power = peak(curve, *window)
        ogle.append(
            {
                "id": star,
                "receipt": receipt,
                "time_scale": scale,
                "points": len(arr),
                "usable": len(curve[0]),
                "baseline_days": float(np.ptp(curve[0])),
                "known_period_days": known,
                "at_known_period": fit(curve, 1 / known, curve[0][0]),
                "tess_window_best_days": 1 / narrow,
                "tess_window_best_power": power,
            }
        )
    out = {
        "source_id": source["source_id"],
        "center_deg": [ra, dec],
        "radius_arcsec": 30,
        "tess_periods_days": periods,
        "frequency_window_per_day": list(window),
        "wide_frequency_per_day": [0.1, 5.0],
        "source_rows": len(results),
        "injections": {
            "amplitudes_mag": [0.01, 0.03, 0.1],
            "injected_period_days": 1 / center,
            "paired_stars": sum("injected_recovery" in r for r in results),
            "recovered_counts": {
                str(a): sum(r.get("injected_recovery", {}).get(str(a), False) for r in results)
                for a in (0.01, 0.03, 0.1)
            },
            "rule": "both holdout bands and early red >=5% variance, <=0.2-cycle drift and early frequency within one training Rayleigh width; same sinusoidal magnitude amplitude in both EROS bands",
        },
        "exposure_rows": len(rows),
        "input": catalog_receipt,
        "exposures": lc_receipt,
        "sources": results,
        "ogle": ogle,
    }
    atomic_bytes(args.output, (json.dumps(out, indent=2, allow_nan=False) + "\n").encode())


if __name__ == "__main__":
    main()
