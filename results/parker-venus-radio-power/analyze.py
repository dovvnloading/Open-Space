"""Measure two-sided log-power contrasts at Parker's closest approaches to Venus."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from itertools import permutations

HERE = Path(__file__).resolve().parent
BANDS = ["db_1.3_2", "db_2_4", "db_4_8", "db_8_16"]
ANTENNAS = ["ch0_V1V2", "ch1_V3V4"]


def contrast(arr, center, half_width, near, far):
    # Integer minute label denotes the beginning of a one-minute bin.
    offset = np.arange(arr.shape[0]) + 0.5 - center
    masks = [np.abs(offset) <= half_width, (offset >= -far) & (offset < -near),
             (offset > near) & (offset <= far)]
    meds = []
    for mask in masks:
        window = arr[mask]
        if len(window) < 5 or np.any(np.isfinite(window).sum(axis=0) < 0.80 * mask.sum()):
            return None
        meds.append(np.median(window, axis=0) if np.isfinite(window).all() else np.nanmedian(window, axis=0))
    return (meds[0] - (meds[1] + meds[2]) / 2, meds[0] - meds[1], meds[0] - meds[2])


def load_minutes(path):
    series = {}
    with gzip.open(path, "rt") as stream:
        for row in csv.DictReader(stream):
            key = row["date"], row["antenna"], int(row["nsum"])
            if key not in series:
                series[key] = np.full((1440, 4), np.nan)
            series[key][int(row["utc_minute"])] = [float(row[b]) if row[b] else np.nan for b in BANDS]
    return series


def rounded(vals):
    return [round(float(n), 4) for n in vals]


def analyze(series, manifest):
    results = []
    for event in manifest:
        day = event["date"]
        center = event["maximum_radius_utc_minute"]
        rec = {"date":day, "center_minute": center, "radius_deg":event["maximum_radius_deg"]}
        ants = {}
        for ant in ANTENNAS:
            for nsum in [40, 80]:
                arr = series.get((day, ant, nsum))
                if arr is None:
                    continue
                val = contrast(arr, center, 3, 20, 50)
                if val is None:
                    continue
                # Half-hour spaced control centers on the same UTC day and same receiver setup.
                controls = []
                for m in range(60, 1381, 30):
                    if abs(m-center) <= 100:
                        continue
                    ctrl = contrast(arr, m, 3, 20, 50)
                    if ctrl is not None:
                        controls.append(ctrl[0])
                ctl = np.array(controls).reshape(-1, 4)
                sensitivity = {}
                for hw, nr, fr in [(5,20,50),(3,15,30),(3,30,60)]:
                    alt = contrast(arr, center, hw, nr, fr)
                    sensitivity[f"{hw}_{nr}_{fr}"] = rounded(alt[0]) if alt is not None else None
                ants[ant] = dict(nsum=nsum, contrast_db=rounded(val[0]), vs_before_db=rounded(val[1]),
                                 vs_after_db=rounded(val[2]), matched_control_centers=len(controls),
                                 controls_as_low=(ctl <= val[0]).sum(axis=0).tolist(), sensitivity_db=sensitivity)
                break
        rec["antennas"] = ants
        results.append(rec)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minutes", default=HERE/"minutes.csv.gz",type=Path)
    args=parser.parse_args()
    results=analyze(load_minutes(args.minutes),json.loads((HERE/"manifest.json").read_text()))
    valid=[r for r in results if len(r["antennas"])==2]
    radius=[r["radius_deg"] for r in valid]
    # 4–8 MHz; mean of the simultaneous dipoles is one event, not two replicates.
    y=np.array([np.mean([x["contrast_db"][2] for x in r["antennas"].values()]) for r in valid])
    rho=spearmanr(radius,y).statistic
    perm=[spearmanr(radius, s).statistic for s in permutations(y)]
    summary={"events":results,"four_to_eight_mhz":{"usable_flybys":len(y),"spearman_rho":round(float(rho),5),
             "two_sided_exact_permutations_as_extreme":sum(bool(abs(s)+1e-10>=abs(rho)) for s in perm),"permutations":len(perm)}}
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
