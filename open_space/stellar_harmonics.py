"""Check an measured stellar period and its double across independent TESS sectors."""

import argparse
import json
import math
from pathlib import Path

from .stellar import binned
from .tess_archive import atomic_bytes, cached_product, inventory


def analyze_sector(times, flux, period):
    import numpy as np
    from astropy.timeseries import LombScargle

    span = float(np.ptp(times))
    if span < 15 or not math.isfinite(period) or period <= 0:
        raise ValueError("Need a positive period and 15 days of observations")
    grid = np.arange(0.08, 8.0, 1 / (10 * span))
    ls = LombScargle(times, flux, fit_mean=True, center_data=True)
    power = ls.power(grid)
    i = int(np.argmax(power))
    return {
        "start_btjd": float(times[0]),
        "end_btjd": float(times[-1]),
        "bins": len(times),
        "baseline_days": span,
        "frequency_resolution_per_day": 1 / span,
        "frequency_trials": len(grid),
        "wide_peak_period_days": float(1 / grid[i]),
        "wide_peak_variance_fraction": float(power[i]),
        "power_at_measured_period": float(ls.power(1 / period)),
        "power_at_twice_measured_period": float(ls.power(1 / (2 * period))),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tic", required=True, type=int)
    p.add_argument("--period", required=True, type=float)
    p.add_argument("--sectors", required=True, nargs="+", type=int)
    p.add_argument("--discovery-sectors", required=True, nargs="+", type=int)
    p.add_argument("--catalog-dir", required=True, type=Path)
    p.add_argument("--cache", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    sectors = sorted(set(args.sectors))
    if not set(args.discovery_sectors).issubset(sectors):
        raise ValueError("Discovery sectors must be included in --sectors")
    args.catalog_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for sector in sectors:
        products, source = inventory(sector, args.catalog_dir)
        name = products[args.tic]
        product = {"sector": sector, "filename": name, "uri": "mast:TESS/product/" + name}
        path, receipt = cached_product(product, args.cache)
        curve = binned(path)
        if curve is None:
            raise ValueError(f"Insufficient data in sector {sector}")
        results.append(
            {
                "role": "discovery" if sector in args.discovery_sectors else "holdout",
                "inventory": source,
                "product": receipt,
                **analyze_sector(*curve, args.period),
            }
        )
    out = {
        "tic_id": args.tic,
        "measured_period_days": args.period,
        "search_frequency_per_day": [0.08, 8],
        "method": "30-minute PDCSAP medians, floating-mean one-term Lomb-Scargle, 10 samples per resolution element",
        "sectors": results,
    }
    atomic_bytes(args.output, (json.dumps(out, indent=2, allow_nan=False) + "\n").encode())


if __name__ == "__main__":
    main()
