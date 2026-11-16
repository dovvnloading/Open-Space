"""Compare USC_SEAPRISM with older, quality-assured daily field radiances at matched bands."""
import argparse
import hashlib
from pathlib import Path
import subprocess
from urllib.parse import urlencode

import numpy as np
import pandas as pd

YEARS = range(2010, 2026)
BASE = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"


def run(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    all_days, manifest = [], []
    for year in YEARS:
        query = urlencode(dict(site="USC_SEAPRISM", year=year, month=2, day=1, year2=year,
                               month2=5, day2=31, LWN20=1, AVG=20, if_no_html=1))
        url = f"{BASE}?{query}"
        path = cache / f"usc_aoc_v3_l2_{year}.csv"
        if not path.exists():
            raw = subprocess.check_output(["curl", "-fsSL", "--retry", "2", "--max-time", "90", url])
            if not raw.startswith(b"AERONET Data Download"):
                raise ValueError(f"Unrecognized response for {year}")
            path.write_bytes(raw)
        raw = path.read_bytes()
        if not raw.startswith(b"AERONET Data Download"):
            raise ValueError(f"Unrecognized input {path}")
        receipt = {"year": year, "url": url, "sha256": hashlib.sha256(raw).hexdigest(), "days": 0}
        manifest.append(receipt)
        if b"Date(dd-mm-yyyy)" not in raw:
            if b"No records" not in raw and raw.strip() != b"AERONET Data Download (Version 3 SEAPRISM Processing)":
                raise ValueError(f"Unrecognized empty release for {year}")
            continue
        d = pd.read_csv(path, skiprows=5, na_values=-999)
        receipt["days"] = len(d)
        day = pd.DataFrame({"date": pd.to_datetime(d["Date(dd-mm-yyyy)"], format="%d:%m:%Y")})
        # Use the same nominal channels at every epoch for the historical comparison.
        # Delta log ratios cancel fixed per-channel conversion constants (e.g. solar irradiance).
        for red in (560, 667):
            bands = d[[f"Lwn_f/Q[{w}nm]" for w in [443, 490, red]]]
            logs = np.log(bands.where(bands > 0))
            day[f"C{red}"] = logs.iloc[:, 1] - ((red - 490) * logs.iloc[:, 0] + 47 * logs.iloc[:, 2]) / (red - 443)
        all_days.append(day)
    data = pd.concat(all_days).sort_values("date").set_index("date")
    pd.DataFrame(manifest).to_csv(out / "field_archive_manifest.csv", index=False)
    results = []
    controls = []
    for start, end in [("2024-03-11", "2024-03-15"), ("2024-04-10", "2024-04-15")]:
        first, last = pd.Timestamp(start), pd.Timestamp(end)
        days = (last - first).days
        # A common leap calendar preserves month/day across leap and nonleap years.
        phase = data.index.map(lambda d: d.replace(year=2000))
        past = data.loc[(data.index.year < last.year) & (abs((phase-last.replace(year=2000)).days) <= 30)]
        delta = data["C667"] - data["C667"].shift(freq=f"{days}D")
        before = delta.reindex(past.index).dropna()
        observed = data.at[last, "C667"] - data.at[first, "C667"]
        ctr = pd.DataFrame({"control_end": before.index, "control_delta_443_490_667": before.values})
        ctr["event_end"] = end
        controls.append(ctr)
        results.append({"start": start, "end": end, "days": days, "daily_delta_443_490_560":
                        data.at[last, "C560"] - data.at[first, "C560"],
                        "daily_delta_443_490_667": observed, "historical_endpoint_days": past.C667.count(),
                        "historical_endpoint_years": past.C667.dropna().index.year.nunique(),
                        "historical_exact_lag_pairs": len(before), "historical_delta_min": before.min(),
                        "historical_delta_median": before.median(), "historical_as_or_more_negative": int((before <= observed).sum())})
    pd.DataFrame(results).to_csv(out / "field_history_summary.csv", index=False, float_format="%.9g")
    pd.concat(controls).to_csv(out / "field_control_intervals.csv", index=False, float_format="%.9g")
    data.loc["2024-03-07":"2024-04-19"].dropna().to_csv(out / "field_daily_shape.csv", float_format="%.9g")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    run(args.cache, args.out)
