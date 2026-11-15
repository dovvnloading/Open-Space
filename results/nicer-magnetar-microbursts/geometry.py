"""Test the event-time gradient on the physical detector array. Requires NumPy."""
import csv
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data" / "events.csv"


def plane_test(rows, obs):
    signal = [r for r in rows if r["obsid"]==obs and abs(float(r["delta_us"]))<=10
              and r["cleaned"]=="1" and 50<=int(r["pi"])<=800]
    t = np.array([float(r["delta_us"]) for r in signal])
    xy = np.array([[float(r["rawx"]),float(r["rawy"])] for r in signal])
    mpu = np.array([int(r["det_id"])//10 for r in signal])
    xy -= xy.mean(axis=0)
    t -= t.mean()
    projection = xy @ np.linalg.pinv(xy)
    total = t @ t

    def score(times):
        return ((times @ projection)*times).sum(axis=1)/total

    observed = score(t[None,:])[0]
    rng = np.random.default_rng(20261114)
    exceed = 0
    groups = [np.flatnonzero(mpu==value) for value in np.unique(mpu)]
    for _ in range(100):
        draws = np.tile(t,(1000,1))
        for group in groups:
            draws[:,group] = rng.permuted(draws[:,group],axis=1)
        exceed += int(np.sum(score(draws)>=observed-1e-10))
    print(f"{obs}  R²={observed:.5f}, {exceed}/100000 within-MPU permutations matched or exceeded it")
    return observed, exceed


def footprint_test(rows):
    """Compare the 2018 bounding rectangle with weighted draws of distinct modules."""
    flash = [r for r in rows if r["obsid"]=="1020240141" and abs(float(r["delta_us"]))<=10]
    ids = np.array([int(r["det_id"]) for r in flash])
    xy = np.array([[int(r["rawx"]),int(r["rawy"])] for r in flash])
    area = int(np.prod(np.ptp(xy,axis=0)+1))
    with (DATA.parent / "detector_counts_2018.csv").open() as f:
        reference = list(csv.DictReader(f))
    refid = np.array([int(r["det_id"]) for r in reference])
    counts = np.array([int(r["count"]) for r in reference])
    xy = np.array([[int(r["rawx"]),int(r["rawy"])] for r in reference])
    rng = np.random.default_rng(20261115)
    xs, ys = [], []
    for unit in np.unique(ids//10):
        n = int(np.sum(ids//10==unit))
        take = refid//10==unit
        weight = counts[take]/counts[take].sum()
        # Exponential races give weighted sampling without replacement.
        keys = -np.log(rng.random((100000,len(weight)))) / weight
        draw = np.argpartition(keys,n-1,axis=1)[:,:n]
        xs.append(xy[take,0][draw])
        ys.append(xy[take,1][draw])
    areas = (np.ptp(np.hstack(xs),axis=1)+1)*(np.ptp(np.hstack(ys),axis=1)+1)
    print(f"2018 footprint {area} cells, {np.sum(areas<=area)}/100000 matched draws as compact")


if __name__=="__main__":
    with DATA.open() as file:
        events = list(csv.DictReader(file))
    for obsid in ("1020240141","5505010502"):
        plane_test(events, obsid)
    footprint_test(events)
