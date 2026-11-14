"""Recreate the complete discovery screen from pinned HARPS and host positions."""
import argparse
import hashlib
import io
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.spatial import cKDTree

from verify_sources import ARCHIVE_SHA, HARPS_URL

ROOT = Path(__file__).resolve().parent


def vectors(table):
    ra, dec = np.deg2rad(table[["ra", "dec"]].to_numpy()).T
    return np.column_stack([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])


def screen(archive):
    if hashlib.sha256(archive).hexdigest() != ARCHIVE_SHA:
        raise ValueError("HARPS archive digest differs from the discovery input")
    columns = ["target", "BJD", "f_RV", "RV_mlc_nzp", "e_RV_mlc_nzp", "ra", "dec"]
    h = pd.read_csv(io.BytesIO(archive), compression="zip", usecols=columns)
    good = (h.f_RV == 0) & (h.RV_mlc_nzp.abs() < 50000) & h.e_RV_mlc_nzp.between(.01, 30)
    h = h[good]
    stars = h.groupby("target").first().reset_index()
    host_snapshot = (ROOT/"host_positions.csv.gz").read_bytes()
    expected_host_sha = "cee95a7d5fdff5941841b4f092f9b0a28bb2df26c59242657f7cb5b49ec27f91"
    if hashlib.sha256(host_snapshot).hexdigest() != expected_host_sha:
        raise ValueError("Confirmed-host position snapshot changed")
    hosts = pd.read_csv(io.BytesIO(host_snapshot), compression="gzip")
    distance, _ = cKDTree(vectors(hosts)).query(vectors(stars))
    excluded = set(stars.target[distance < 15/206265])
    rows = []
    for name, data in h.groupby("target"):
        if name in excluded or len(data) < 28:
            continue
        nightly = data.groupby(data.BJD.round().rename("night")).median(numeric_only=True).sort_values("BJD")
        if len(nightly) < 28:
            continue
        t = nightly.BJD.to_numpy()
        y = nightly.RV_mlc_nzp.to_numpy(copy=True)
        post = t > 2457165
        for era in np.unique(post):
            y[post == era] -= np.median(y[post == era])
        median = np.median(y)
        deviation = np.abs(y - median)
        keep = deviation < max(6*np.median(deviation), 10)
        t, y = t[keep], y[keep]
        if len(t) < 24 or np.ptp(t) < 25:
            continue
        freq, power = LombScargle(t, y).autopower(
            minimum_frequency=max(.0007, 1.7/np.ptp(t)), maximum_frequency=.8,
            samples_per_peak=5)
        peak = np.nanargmax(power)
        period = 1/freq[peak]
        score = -len(t)*np.log(1-power[peak]) - 3*np.log(len(t))
        phase = 2*np.pi*t/period
        model = np.column_stack([np.ones(len(t)), np.sin(phase), np.cos(phase)])
        coefficients = np.linalg.lstsq(model, y, rcond=None)[0]
        rows.append({"target":name, "nights":len(t), "span_days":round(float(np.ptp(t))),
                     "period_days":round(float(period),5), "power":round(float(power[peak]),5),
                     "screen_pass":int(power[peak] >= .20), "score":round(float(score),3),
                     "K_m_s":round(float(np.hypot(*coefficients[1:])),3),
                     "residual_rms_m_s":round(float(np.std(y-model@coefficients)),3)})
    return pd.DataFrame(rows).sort_values("score", ascending=False, kind="stable")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Local copy of HARPS_RVBank_ver02.csv.zip")
    parser.add_argument("--output", type=Path, required=True, help="Destination for the full rank table")
    options = parser.parse_args()
    if options.archive:
        data = options.archive.read_bytes()
    else:
        with urlopen(HARPS_URL, timeout=90) as response:
            data = response.read()
    result = screen(data)
    result.to_csv(options.output, index=False)
    print(f"{len(result)} targets, {result.screen_pass.sum()} above the power threshold")
