"""Fit independently measured Gaia clocks on early ZTF r photometry, predict late r photometry."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

HARMONICS = 4


def design(time, freq):
    phase = 2 * np.pi * np.outer(time, np.arange(1, HARMONICS + 1)) * freq
    return np.column_stack([np.ones(len(time)), np.sin(phase), np.cos(phase)])


def fit(time, mag, frequency):
    x = design(time, frequency)
    coef = np.linalg.lstsq(x, mag, rcond=None)[0]
    return np.sum((mag - x @ coef) ** 2), coef


def prediction(time, mag, own_period, other_period):
    t = time - time[0]
    # Equal time-span split, so both train and test are multi-season when sampling allows.
    train = t <= (t[-1] / 2)
    valid = ~train
    if min(sum(train),sum(valid)) < 40 or t[-1] < 365:
        return {"status": "insufficient coverage"}
    result = {"status":"ok","n_train":int(sum(train)), "n_test":int(sum(valid)), "span_days":float(t[-1])}
    null_sse = np.sum((mag[valid] - np.mean(mag[train]))**2)
    for key, period in [("own", own_period), ("neighbor", other_period)]:
        frequency = 1 / period
        # +/- two cycles over the full time span. Only the first half sets the frequency and waveform.
        grid = frequency + np.linspace(-2, 2, 121) / t[-1]
        fits = [fit(t[train], mag[train], f) for f in grid]
        at = int(np.argmin([s for s,c in fits]))
        best_freq = grid[at]
        coefs = fits[at][1]
        resid = mag[valid] - design(t[valid], best_freq) @ coefs
        model = design(np.arange(500)/500,1) @ coefs
        result[key + "_period"] = 1 / best_freq
        result[key + "_test_r2"] = 1 - np.sum(resid**2) / null_sse
        result[key + "_amplitude_mag"] = float(np.ptp(model))
    return result


def measure(data, band, own, neighbor):
    raw = data[(data["catflags"]==0) & (data["filtercode"]==band)].dropna(subset=["hjd","mag","magerr"])
    raw = raw[(raw.magerr > 0) & (raw.magerr < .1)]
    if not len(raw):
        return {"status":"no clean photometry"}
    oid = raw.oid.value_counts().index[0]
    d = raw[raw.oid==oid].sort_values("hjd").drop_duplicates("hjd")
    # Symmetric broad clip of isolated bad measurements, preserving eclipse amplitudes.
    mag = d.mag.to_numpy()
    med = np.median(mag)
    keep = np.abs(mag - med) < max(1.5, 8 * np.median(abs(mag-med)))
    d = d[keep]
    r = prediction(d.hjd.to_numpy(), d.mag.to_numpy(),own,neighbor)
    r.update(oid=str(oid), ra_median=float(d.ra.median()),dec_median=float(d.dec.median()))
    return r


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("gaia_join", type=Path)
    p.add_argument("epochs", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    cohort = [r for r in csv.DictReader(a.gaia_join.open()) if
              float(r["sep_arcsec"]) >= 5 and float(r["dec1"]) >= -27
              and float(r["R_chance_align"]) <= .01]
    epochs = pd.read_csv(a.epochs, dtype={"source_id": str})
    out = {}
    for pair in cohort:
        for i in (1, 2):
            sid = pair[f"source_id{i}"]
            data = epochs[epochs.source_id == sid]
            out[sid] = {b: measure(data,b,float(pair[f"period{i}"]),float(pair[f"period{3-i}"]))
                        for b in ["zr","zg"]}
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
