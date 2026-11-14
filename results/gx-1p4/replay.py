"""Replay timing and flux comparisons from the included mission tables."""
import csv
from datetime import date
from itertools import combinations
import json
import math
from pathlib import Path
from statistics import mean, median

ROOT = Path(__file__).resolve().parent
DAYS_PER_YEAR = 365.25
MJD_EPOCH = date(1858, 11, 17)
SPRING = (61111, 61176)


def load(name):
    with (ROOT / name).open(newline="") as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]


def between(rows, start, stop, key="barytime"):
    return [row for row in rows if start <= row[key] < stop]


def fit(rows):
    times = [row["barytime"] for row in rows]
    frequencies = [1e6 * row["frequency"] for row in rows]
    t0, f0 = mean(times), mean(frequencies)
    return DAYS_PER_YEAR * sum((t-t0)*(f-f0) for t, f in zip(times, frequencies)) / sum((t-t0)**2 for t in times)


def delta(first, last):
    return 1e6 * (last["frequency"] - first["frequency"])


def pairs(rows, lo, hi):
    return [(b["barytime"]-a["barytime"], delta(a,b), a["barytime"], b["barytime"])
            for a, b in combinations(rows, 2) if lo <= b["barytime"]-a["barytime"] <= hi]


def gbm_analysis(rows):
    rows = sorted((r for r in rows if r["detected"] and r["frequency_err"] > 0),
                  key=lambda r: r["barytime"])
    event = between(rows, *SPRING)
    first, last = event[0], event[-1]
    days = last["barytime"] - first["barytime"]
    pre = [r for r in rows if r["barytime"] < 61088]
    ranges = [(18,90), (24,80), (30,100)]
    history = {}
    for lo, hi in ranges:
        v = pairs(pre, lo, hi)
        history[f"{lo}-{hi}"] = {"pairs":len(v), "positive":sum(u[1]>0 for u in v),
                                 "largest_change_uHz":max(u[1] for u in v)}
    annual = {}
    for year in range(2008, 2027):
        start = (date(year,1,1)-MJD_EPOCH).days
        end = (date(year+1,1,1)-MJD_EPOCH).days
        data = between(rows, start, end)
        if len(data) >= 2:
            annual[year] = {"n":len(data), "span_days":data[-1]["barytime"]-data[0]["barytime"],
                            "uHz_per_year":fit(data)}
    # Keplerian upper bound for any orbital phase, with neutron-star K=30 km/s.
    vmax_change = 30 * 2 * math.pi * days / 1160.8 * (1+0.10)**2 / (1-0.10**2)**1.5
    allpairs = pairs(rows,24,80)
    strongest = max(allpairs, key=lambda p:p[1])
    followup = {}
    for start, stop in [(61175,61297), (61296,61320), (61319,61352), (61175,61352)]:
        sample = between(rows,start,stop)
        initial, final = sample[0], sample[-1]
        followup[f"{start}-{stop}"] = {
            "epochs":len(sample),"first_mjd":initial["barytime"],"last_mjd":final["barytime"],
            "span_days":final["barytime"]-initial["barytime"], "delta_uHz":delta(initial,final),
            "sigma_endpoints_uHz":1e6*math.hypot(initial["frequency_err"],final["frequency_err"]),
            "uHz_per_year":fit(sample)}
    return {"rows":len(rows), "annual_fits":annual, "post_alert":followup,
            "2009_2022_uHz_per_year":fit(between(rows,54832,59945)),
            "pre_february_2026_pair_comparison":history,
            "all_24_80_day_pairs":len(allpairs), "largest_pair":strongest,
            "event":{"epochs":len(event), "first_mjd":first["barytime"],"last_mjd":last["barytime"],
                     "span_days":days,"delta_uHz":delta(first,last),
                     "sigma_endpoints_uHz":1e6*math.hypot(first["frequency_err"],last["frequency_err"]),
                     "uHz_per_year":fit(event),"odd_even_uHz_per_year":[fit(event[::2]),fit(event[1::2])],
                     "leave_one_out_uHz_per_year":[fit(event[:k]+event[k+1:]) for k in range(len(event))],
                     "min_search_edge_margin_uHz":min(1e6*min(r["frequency"]-r["fmin"],r["fmax"]-r["frequency"]) for r in event),
                     "min_pulsed_amplitude_snr":min(r["amplitude"]/r["ampsig"] for r in event),
                     "min_profiles":min(r["nprofiles"] for r in event),
                     "max_orbital_delta_uHz":1e6*max(r["frequency"] for r in event)*vmax_change/299792.458},
            "peak_to_final_delta_uHz":delta(last,rows[-1]),"last_mjd":rows[-1]["barytime"]}


def maxi_windows(rows):
    result = {}
    for lo,hi in [(61046,61111),SPRING,(61176,61241)]:
        window = [r for r in between(rows,lo,hi,"mjd") if r["error_2_20"]>0]
        weighted = sum(r["flux_2_20"]/r["error_2_20"]**2 for r in window) / sum(r["error_2_20"]**-2 for r in window)
        result[f"{lo}-{hi}"] = {"days":len(window),"median_ph_cm2_s":median(r["flux_2_20"] for r in window),
                                  "weighted_mean_ph_cm2_s":weighted}
    return result


if __name__ == "__main__":
    print(json.dumps({"gbm":gbm_analysis(load("gbm.csv")),"maxi":maxi_windows(load("maxi.csv"))}, indent=2))
