"""Recompute the Jovian energy screen and time-resolved Horizons measurements."""
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
AU = 149597870.7
GM = 126686534.0  # Jupiter alone, km^3/s^2. Horizons states use its center.
YEAR = 365.25


def screen(catalog):
    rows = [dict(zip(catalog["fields"], v)) for v in catalog["data"]]
    bound = Counter()
    for row in rows:
        energy = float(row["v_rel"])**2/2 - GM/(AU*float(row["dist"]))
        if energy < 0:
            bound[row["des"]] += 1
    return {"encounters":len(rows),"objects":len({x["des"] for x in rows}),
            "negative_encounters":bound.total(),"negative_by_object":dict(sorted(bound.items()))}


def events(states):
    t = states[:,0]
    distance = np.linalg.norm(states[:,1:4],axis=1)
    e = np.sum(states[:,4:]**2,axis=1)/2 - GM/distance
    angle = np.unwrap(np.arctan2(states[:,2],states[:,1]))
    results = []
    inside = np.where(e < 0)[0]
    for run in np.split(inside, np.where(np.diff(inside)>1)[0]+1):
        if not len(run):
            continue
        lo, hi = run[0], run[-1]
        if lo == 0 or hi == len(t)-1:
            raise ValueError("Energy event reaches the ephemeris window edge")
        start = np.interp(0,[e[lo],e[lo-1]],[t[lo],t[lo-1]])
        end = np.interp(0,[e[hi],e[hi+1]],[t[hi],t[hi+1]])
        theta_start,theta_end = np.interp([start,end],t,angle)
        dates = [str((datetime(2000,1,1,12)+timedelta(days=float(v-2451545))).date())
                 for v in (start,end)]
        results.append({"energy_interval_tdb":dates,
                        "duration_years":round((end-start)/YEAR,4),
                        "min_distance_au":round(float(distance[run].min())/AU,6),
                        "min_energy_km2_s2":round(float(min(e[run])),4),
                        "signed_turns_j2000_xy":round((theta_end-theta_start)/(2*math.pi),4)})
    return {"events":results, "window_min_energy_km2_s2":round(float(min(e)),4)}


def load_tracks():
    groups = defaultdict(list)
    with gzip.open(ROOT/"data"/"horizons.csv.gz","rt") as source:
        rows = csv.reader(source)
        next(rows)
        for name,*coordinates in rows:
            groups[name].append([float(v) for v in coordinates])
    return {k:np.array(v) for k,v in groups.items()}


if __name__ == "__main__":
    data = ROOT/"data"
    result = {"screen":screen(json.loads((data/"encounters.json").read_text())),
              "horizons":{name:events(states) for name,states in load_tracks().items()}}
    print(json.dumps(result,indent=2))
