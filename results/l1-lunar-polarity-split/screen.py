"""Re-run the exploratory hourly ACE/Wind event-selection screen (numpy, pandas)."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

FIELDS = {"ace": ("AC_H2_MFI", "BGSM,SC_pos_GSE"), "wind": ("WI_H0_MFI@2", "B1GSM,P1GSE")}
ROOT = Path(__file__).resolve().parent


def fetch_or_read(job):
    mission, year = job
    dataset, fields = FIELDS[mission]
    query = {"id": dataset, "parameters": fields, "format": "csv", "time.min": f"{year}-01-01T00:00:00Z",
             "time.max": f"{year+1}-01-01T00:00:00Z"}
    path = args.cache / f"{mission}{year}.csv"
    if not path.exists():
        with urllib.request.urlopen("https://cdaweb.gsfc.nasa.gov/hapi/data?" + urllib.parse.urlencode(query), timeout=150) as response:
            content = response.read()
        if not content.startswith(str(year).encode()) and b"no data for time range" not in content:
            raise ValueError(content[:200])
        path.write_bytes(content)
    content = path.read_bytes()
    if content.startswith(b"{") and b"no data for time range" in content:
        return mission, year, pd.Series(dtype=float), hashlib.sha256(content).hexdigest()
    t = pd.read_csv(path, header=None, usecols=[0, 3], names=["utc", mission])
    t.index = pd.to_datetime(t.utc, utc=True).dt.floor("h")
    return mission, year, t[mission].where(t[mission].abs() < 200), hashlib.sha256(content).hexdigest()


def screen():
    jobs = [(mission, year) for mission in FIELDS for year in range(2015, 2026)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(fetch_or_read, jobs))
    series = [pd.concat([s for m, _, s, _ in rows if m == key and not s.empty]).rename(key) for key in FIELDS]
    paired = pd.concat(series, axis=1).sort_index().asfreq("h")
    summary = {"screen_pairs_with_both_bz":int(paired.dropna().shape[0]),
               "input_rows_and_sha256":{f"{m}_{year}":[int(s.size),sha] for m, year, s, sha in rows}}
    scores = []
    for hours in [1,2,3]:
        a = paired.ace.rolling(hours).mean()
        w = paired.wind.rolling(hours).mean()
        for direction in [1,-1]:
            by_shift = pd.concat([direction*(a-w.shift(i)) for i in [-1,0,1]],axis=1)
            # Retain only a discrepancy with the same orientation in all three clock alignments.
            v = by_shift.min(axis=1,skipna=False)
            for utc, value in v.dropna().nlargest(1).items():
                scores.append((hours,utc,float(value),direction,int(v.count())))
    summary["top_by_window"] = []
    for hours in [1,2,3]:
        _,utc,score,direction,count = max((x for x in scores if x[0]==hours), key=lambda x:x[2])
        summary["top_by_window"].append(dict(hours=hours,utc=str(utc),score_nT=round(score,3),
                                             ace_minus_wind_sign=direction,valid_windows=count))
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,required=True)
    args = p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    print(json.dumps(screen(),indent=2))
