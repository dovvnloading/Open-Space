"""Search raw stellar variability on two TESS sectors, independent of the transit filter."""

import argparse
import csv
import io
import json
import math
from pathlib import Path
from urllib.parse import urlencode

from .archive_index import retrieve, record_response
from .tess_archive import cached_product


def binned(path):
    import numpy as np
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdus:
        tab = hdus[1].data
        t = np.asarray(tab["TIME"], dtype=float)
        f = np.asarray(tab["PDCSAP_FLUX"], dtype=float)
        err = np.asarray(tab["PDCSAP_FLUX_ERR"], dtype=float)
        ok = (
            (tab["QUALITY"] == 0)
            & np.isfinite(t)
            & np.isfinite(f)
            & np.isfinite(err)
            & (err > 0)
            & (f > 0)
        )
        if ok.sum() < 200:
            return None
        t, f = t[ok], f[ok]
        if not np.all(np.diff(t) > 0):
            raise ValueError("Nonmonotonic time in stellar light curve")
        scale = np.median(f)
        bins = ((t - t[0]) * 48).astype(int)
        _, starts, n = np.unique(bins, return_index=True, return_counts=True)
        out = [
            (np.median(t[i : i + k]), np.median(f[i : i + k]) / scale - 1)
            for i, k in zip(starts, n)
            if k >= 4
        ]
        if len(out) < 180:
            return None
        time, val = np.asarray(out).T
        return time, val


def search(curves):
    import numpy as np
    from astropy.timeseries import LombScargle

    if any(c is None for c in curves):
        return {"status": "insufficient_coverage"}
    (t, y), (tv, yv) = curves
    span = t[-1] - t[0]
    if min(span, tv[-1] - tv[0]) < 15:
        return {"status": "insufficient_baseline"}
    grid = np.arange(0.2, 8.0, 1 / (8 * span))

    def periodogram(times, vals):
        return LombScargle(times, vals, center_data=True, fit_mean=True).power(
            grid, method="fast"
        )

    if min(np.ptp(y), np.ptp(yv)) == 0:
        return {"status": "flat_light_curve"}
    train, valid = periodogram(t, y), periodogram(tv, yv)
    i, j = int(np.argmax(train)), int(np.argmax(valid))
    freq = float(grid[i])

    def fit(times, vals):
        phase = (times - t[0]) * freq * 2 * np.pi
        d = np.column_stack((np.cos(phase), np.sin(phase), np.ones(len(phase))))
        a, b, c = np.linalg.lstsq(d, vals, rcond=None)[0]
        model = a * d[:, 0] + b * d[:, 1] + c
        r2 = 1 - np.sum((vals - model) ** 2) / np.sum((vals - np.mean(vals)) ** 2)
        return float(np.hypot(a, b)), float(np.arctan2(b, a)), float(r2)

    amp, ph, r1 = fit(t, y)
    av, pv, r2 = fit(tv, yv)
    offset = abs(float(grid[j]) - freq)
    drift = abs(math.remainder(pv - ph, 2 * math.pi)) / (2 * math.pi)
    ratio = av / amp if amp > 0 else None
    passed = (
        r1 >= 0.15
        and r2 >= 0.10
        and offset <= 1 / max(span, tv[-1] - tv[0])
        and drift <= 0.15
        and ratio is not None
        and 0.4 <= ratio <= 2.5
    )
    return {
        "status": (
            "two_sector_coherent_variability"
            if passed
            else "below_variability_threshold"
        ),
        "period_days": 1 / freq,
        "epoch_min_btjd": float(t[0] + ((ph / (2 * math.pi) + 0.5) % 1) / freq),
        "validation_peak_period_days": 1 / float(grid[j]),
        "phase_offset_cycles": drift,
        "semiamplitude_ppt": [amp * 1000, av * 1000],
        "variance_fraction": [r1, r2],
        "bins": [len(t), len(tv)],
        "frequency_trials": len(grid),
        "method": "sinusoid-v1",
    }


def run_target(target, cache):
    import numpy as np
    from astropy.io import fits

    products = target["products"]
    if len(products) != 2 or products[1]["sector"] <= products[0]["sector"]:
        raise ValueError("Two increasing sectors required")
    checked = [cached_product(p, cache) for p in products]
    if any(int(p[1]["filename"].split("-")[2]) != target["tic_id"] for p in checked):
        raise ValueError("Star differs from requested TIC")
    curves = [binned(path) for path, _ in checked]
    base = search(curves)
    # The validation-site time frequencies are preserved. Permuting nightly flux
    # blocks provides a diagnostic of whether coherent phase survives block reassignment.
    if base["status"] == "two_sector_coherent_variability":
        t, y = curves[1]
        blocks = np.array_split(np.arange(len(t)), max(3, round((t[-1] - t[0]) / 1.0)))
        rng = np.random.default_rng(int(target["tic_id"]))
        n = 32
        survived = 0
        for _ in range(n):
            shuffled = np.concatenate(
                [y[blocks[k]] for k in rng.permutation(len(blocks))]
            )
            survived += (
                search([curves[0], (t, shuffled)])["status"]
                == "two_sector_coherent_variability"
            )
        base["block_control"] = {
            "survived": survived,
            "trials": n,
            "seed": int(target["tic_id"]),
        }
    with fits.open(checked[0][0], memmap=False) as hd:
        ra, dec = map(float, (hd[0].header["RA_OBJ"], hd[0].header["DEC_OBJ"]))
    return {
        **base,
        "tic_id": target["tic_id"],
        "ra_deg": ra,
        "dec_deg": dec,
        "role": target["role"],
        "sources": [record for _, record in checked],
    }


def simbad(ra, dec, cache):
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError("Invalid sky position")
    query = (
        f"SELECT TOP 50 main_id,otype,ra,dec FROM basic WHERE "
        f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.0125))"
    )
    url = "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?" + urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    raw, path = retrieve(url, cache, "text/csv")
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if set(reader.fieldnames or []) != {"main_id", "otype", "ra", "dec"}:
        raise ValueError("SIMBAD lookup incomplete")

    def sep(row):
        ar, dr = math.radians(ra), math.radians(dec)
        br, er = math.radians(float(row["ra"])), math.radians(float(row["dec"]))
        val = (
            math.sin((dr - er) / 2) ** 2
            + math.cos(dr) * math.cos(er) * math.sin((ar - br) / 2) ** 2
        )
        return math.degrees(2 * math.asin(math.sqrt(min(1, val)))) * 3600

    rows = sorted(
        ({**r, "separation_arcsec": sep(r)} for r in reader),
        key=lambda r: r["separation_arcsec"],
    )
    return {
        "records": rows,
        "radius_arcsec": 45,
        "limit": 50,
        "potential_truncation": len(rows) >= 50,
        "receipt": record_response(raw, path, url),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("plan", type=Path)
    p.add_argument("--cache", type=Path, required=True)
    args = p.parse_args()
    rows = [
        run_target(t, args.cache) for t in json.loads(args.plan.read_bytes())["targets"]
    ]
    print(
        json.dumps(
            {"track": "stellar_variability", "targets": rows}, indent=2, allow_nan=False
        )
    )


if __name__ == "__main__":
    main()
