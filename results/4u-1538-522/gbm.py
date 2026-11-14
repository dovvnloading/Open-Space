"""Reconstruct Fermi/GBM secular timing and orbital-phase checks."""
import csv
from itertools import combinations
import json
from math import hypot
from pathlib import Path
from statistics import mean, median

ROOT = Path(__file__).resolve().parent
YEAR = 365.25
ORBIT_DAYS = 3.7283749
ORBIT_EPOCH_MJD = 2451442.489 - 2400000.5


def slope(rows):
    times = [row["barytime"] for row in rows]
    frequencies = [1e6*row["frequency"] for row in rows]
    t0, f0 = mean(times), mean(frequencies)
    return YEAR*sum((t-t0)*(f-f0) for t,f in zip(times,frequencies))/sum((t-t0)**2 for t in times)


def within(rows, lo, hi):
    return [row for row in rows if lo <= row["barytime"] < hi]


def endpoint(rows):
    first,last = rows[0],rows[-1]
    slopes = [1e6*YEAR*(b["frequency"]-a["frequency"])/(b["barytime"]-a["barytime"])
              for a,b in combinations(rows,2) if b["barytime"] != a["barytime"]]
    return {"epochs":len(rows), "first_mjd":first["barytime"], "last_mjd":last["barytime"],
            "span_days":last["barytime"]-first["barytime"],
            "period_endpoints_s":[1/first["frequency"],1/last["frequency"]],
            "delta_frequency_uHz":1e6*(last["frequency"]-first["frequency"]),
            "endpoint_sigma_uHz":1e6*hypot(first["frequency_err"],last["frequency_err"]),
            "frequency_slope_uHz_per_year":slope(rows),"theil_sen_slope_uHz_per_year":median(slopes)}


def orbital_phase(row):
    return ((row["barytime"]-ORBIT_EPOCH_MJD)/ORBIT_DAYS)%1


def analyze():
    with (ROOT/"gbm.csv").open() as stream:
        rows = [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]
    selected = [r for r in rows if r["detected"]==1 and r["frequency_err"]>0]
    # The monitor repeats the identical integration at MJD 58940.67458213849.
    valid = []
    for row in selected:
        if valid and row["barytime"]==valid[-1]["barytime"]:
            if row!=valid[-1]:
                raise ValueError("Conflicting timing measurements at one epoch")
            continue
        valid.append(row)
    epochs = [(55000,58600),(58620,60860),(60900,61354),(61131,61354)]
    fits = {f"{a}-{b}":endpoint(within(valid,a,b)) for a,b in epochs}
    sequence = within(valid,*epochs[-1])
    comparisons = {}
    for width in [180,221,300,365]:
        historical = within(valid,0,60676)
        accepted = []
        for row in historical:
            start = row["barytime"]
            window = within(historical,start,start+width)
            if len(window)>=6 and window[-1]["barytime"]-start>=0.7*width:
                accepted.append(slope(window))
        comparisons[width]={"eligible_windows":len(accepted),"steepest_uHz_per_year":min(accepted)}
    parts = [[r for r in sequence if k/3<=orbital_phase(r)<(k+1)/3] for k in range(3)]
    ranges = [(r["tstop"]-r["tstart"])/ORBIT_DAYS for r in sequence]
    runs=[]
    start=0
    for end in range(1,len(valid)+1):
        if end==len(valid) or valid[end]["frequency"]>=valid[end-1]["frequency"]:
            runs.append((end-start,valid[start]["barytime"],valid[end-1]["barytime"]))
            start=end
    return {"valid":len(valid),"history_rows":len(rows),"intervals":fits,"pre_2025_sliding_windows":comparisons,
            "two_longest_decreasing_runs":sorted(runs,reverse=True)[:2],
            "orbital_thirds":[{"n":len(part),"slope_uHz_per_year":slope(part)} for part in parts],
            "integration_lengths_in_orbits":[min(ranges),max(ranges)],
            "odd_even_slope_uHz_per_year":[slope(sequence[::2]),slope(sequence[1::2])],
            "min_edge_distance_uHz":1e6*min(min(r["frequency"]-r["fmin"],r["fmax"]-r["frequency"]) for r in sequence),
            "search_center_change_uHz":.5e6*((sequence[-1]["fmin"]+sequence[-1]["fmax"])-(sequence[0]["fmin"]+sequence[0]["fmax"]))}


if __name__ == "__main__":
    print(json.dumps(analyze(),indent=2))
