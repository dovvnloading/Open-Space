"""Regenerate the four frozen blind-prediction tables from the pre-QC training fit."""
import csv
from pathlib import Path

import numpy as np

import fit

BASE = Path(__file__).resolve().parent
SETS = [
    ("blind_predictions.csv", "decam_eph.txt", [60]*4),
    ("deep_predictions.csv", "deep_eph.txt", [380,420,175]),
    ("second_predictions.csv", "second_eph.txt", [90]*4),
    ("last_predictions.csv", "last_eph.txt", [175,185]),
]


def reproduce(name, ephemeris, exposure):
    values = (BASE / ephemeris).read_text().split("$$SOE")[1].split("$$EOE")[0]
    data = np.array([[float(row[i]) for i in (0,3,4,5)] for row in csv.reader(values.splitlines()) if len(row)>6])
    jd, ra, dec, dist = data.T
    east, north = fit.tangent_axes(ra, dec)
    observed = list(csv.DictReader((BASE/"prediction_training.csv").open()))
    ref = next(r for r in observed if r["id"] == "HST")
    t0 = float(ref["mjd_utc"])-float(ref["distance_AU"])*499.004783836/86400
    t = jd-2400000.5 + np.array(exposure)/172800 - dist*499.004783836/86400 - t0
    out = []
    for source in list(csv.DictReader((BASE/"prediction_fits.csv").open()))[:5]:
        x = [float(source[k]) for k in fit.PARAMETERS]
        prediction = fit.predict(x, t, east, north, 40/dist)
        for moment, (e,n) in zip(jd, prediction):
            out.append(dict(mjd=round(moment-2400000.5,6), P=round(1/x[0],5),
                            east=round(e,6),north=round(n,6),separation=round(np.hypot(e,n),6),
                            pa_mod180=round(np.degrees(np.arctan2(e,n))%180,5)))
    original=list(csv.DictReader((BASE/name).open()))
    for line, recorded in zip(out, original, strict=True):
        if any(float(recorded[k]) != float(v) for k,v in line.items()):
            raise ValueError(f"{name}: prediction mismatch at {line['mjd']}")
    print(f"Verified {len(out)} forecasts in {name}")


if __name__ == "__main__":
    for schema in SETS:
        reproduce(*schema)
