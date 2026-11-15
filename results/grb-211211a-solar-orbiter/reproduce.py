"""Identify the hard-X-ray source in STIX flare 2112111310 using public FITS."""
import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from numpy.lib.stride_tricks import sliding_window_view

HERE = Path(__file__).resolve().parent


def data_paths(entries, cache):
    cache.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, entry in entries.items():
        path = cache / entry["filename"]
        if not path.exists():
            with urllib.request.urlopen(entry["url"], timeout=90) as response:
                payload = response.read()
            if hashlib.sha256(payload).hexdigest() != entry["sha256"]:
                raise ValueError(f"Download checksum mismatch: {path.name}")
            path.write_bytes(payload)
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Cached checksum mismatch: {path.name}")
        paths[key] = path
    return paths


def subtract_mean_and_normalize(rows):
    centered = rows - rows.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(centered, axis=1, keepdims=True)
    return np.divide(centered, norm, out=np.zeros_like(centered), where=norm > 0)


def gbm_history(table, trigger, grid):
    # TRIGDAT sometimes provides overlapping resolutions. Prefer the finest row.
    start = table["TIME"] - trigger
    end = table["ENDTIME"] - trigger
    rate = table["RATE"][:, 2:5, :12].sum(axis=(1, 2))
    flux = np.full(len(grid), np.nan)
    for index in np.argsort(-(end - start)):
        flux[(grid >= start[index]) & (grid < end[index])] = rate[index]
    return flux


def compare(event, paths):
    trigger = Time(event["mjd"], format="mjd", scale="utc")
    with fits.open(paths[event["fermi"]]) as product:
        data = product["EVNTRATE"].data
        trigger_met = product[0].header["TRIGTIME"]
        nearest = np.argmin(abs(data["TIME"] - trigger_met))
        delta_r = np.asarray(event["solar_orbiter_km"]) - data["EIC"][nearest]
        grid = np.arange(-40, 92, 0.02)
        flux = gbm_history(data, trigger_met, grid)
    ra, dec = np.radians([event["ra"], event["dec"]])
    direction = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
    delay = -delta_r @ direction / 299792.458
    with fits.open(paths[event["stix"]]) as product:
        start = Time(product[0].header["DATE-OBS"], scale="utc")
        data = product["DATA"].data
        time = data["time"].astype(float) / 100 + (start - trigger).sec - delay
        widths = np.asarray(data["timedel"], dtype=float) / 100
        counts = data["counts"].astype(float)
    signal = (time >= -10) & (time < 70)
    t, width, hard = time[signal], widths[signal], counts[signal, 4]
    # Test the same lag grid on signal and on each control block.
    lags = np.arange(-8, 8.01, 0.25)
    templates = []
    for lag in lags:
        model = [np.nanmean(flux[(grid >= c-lag-w/2) & (grid < c-lag+w/2)])
                 for c, w in zip(t, width)]
        if not np.isfinite(model).all():
            raise ValueError("Missing GBM coverage")
        templates.append(model)
    model_vectors = subtract_mean_and_normalize(np.array(templates))
    scores = model_vectors @ (hard - hard.mean())
    best = int(np.argmax(scores))
    zero = list(lags).index(0.0)
    valid = abs(time[:len(time)-len(t)+1]) > 200
    valid &= np.arange(len(valid)) % len(t) == 0
    blocks = sliding_window_view(counts[:, 4], len(t))[valid]
    blocks = blocks[np.isfinite(blocks).all(axis=1)]
    control_scores = (blocks - blocks.mean(axis=1, keepdims=True)) @ model_vectors.T
    control_r = subtract_mean_and_normalize(blocks) @ model_vectors.T
    r = scores[best] / np.linalg.norm(hard - hard.mean())
    bg = (abs(time) >= 80) & (abs(time) < 240)
    summary = dict(id=event["id"], delay_s=round(float(delay), 4),
                   baseline_light_s=round(float(np.linalg.norm(delta_r)/299792.458), 3),
                   lag_offset_s=float(lags[best]), r_best=round(float(r), 5),
                   r_zero_offset=round(float(scores[zero]/np.linalg.norm(hard-hard.mean())), 5),
                   same_day_controls=len(blocks),
                   controls_ge_amplitude=int((control_scores.max(axis=1) >= scores[best]).sum()),
                   controls_ge_r=int((control_r.max(axis=1) >= r).sum()),
                   stix_hard_peak=int(hard.max()), hard_background_median=float(np.median(counts[bg,4])),
                   stix_bin_s=float(np.median(width)),
                   fixed_r_by_stix_band=[round(float(np.corrcoef(templates[zero], values)[0,1]),5)
                                         for values in counts[signal].T])
    rows = [[event["id"], round(float(c), 3), *values.astype(int).tolist(),
             round(float(g), 3)] for c, values, g in zip(t, counts[signal], templates[zero])]
    extra = dict(time=time, counts=counts, trigger_unix=trigger.unix, delay=delay)
    return summary, rows, extra


def check_flare(paths, event, flare):
    t = event["time"] + event["trigger_unix"] + event["delay"]
    counts = event["counts"]
    peak = Time(flare["peak_utc"], scale="utc").unix
    on = (t >= flare["start_unix"]) & (t < flare["end_unix"])
    off = (abs(t - peak) >= 80) & (abs(t - peak) < 240)
    background = counts[off].mean(axis=0)
    excess = counts[on].sum(axis=0) - np.count_nonzero(on) * background
    with fits.open(paths["science"]) as product:
        data = product["DATA"].data
        energy = product["ENERGIES"].data
        center = Time(product[0].header["DATE-OBS"], scale="utc").unix + data["time"] / 100
        science_on = (center >= flare["start_unix"]) & (center < flare["end_unix"])
        science_off = abs(center - peak) > 80
        band = (energy["e_low"] >= 50) & (energy["e_high"] <= 150)
        hard = data["counts"][:, :, :, band].astype(float).sum(axis=(2, 3))
        # All science rows here cover 20 s with identical exposure.
        if not np.all(data["timedel"] == 2000):
            raise ValueError("Unexpected science bin width")
        detector_excess = hard[science_on].sum(axis=0) - sum(science_on)*hard[science_off].mean(axis=0)
    with fits.open(paths["flareflag"]) as product:
        data = product["DATA"].data
        time = Time(product[0].header["DATE-OBS"], scale="utc").unix + data["time"] / 100
        on_flag = (time >= flare["start_unix"]) & (time < flare["end_unix"])
        status = {k: np.unique(data[k][on_flag]).tolist() for k in
                  ["thermal_index","non_thermal_index","location_status","flare_progress"]}
    return dict(flare_id=flare["flare_id"], n_quicklook_centers=int(on.sum()),
                quicklook_background_mean=background.round(2).tolist(),
                quicklook_excess_counts=excess.round(1).tolist(),
                n_science_on=int(science_on.sum()), n_science_off=int(science_off.sum()),
                science_50_150_excess_counts=round(float(detector_excess.sum()),1),
                science_detector_positive=int((detector_excess>0).sum()), flareflag=status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((HERE/"inputs.json").read_text())
    paths = data_paths(manifest["files"], args.cache)
    result, bins = {}, []
    for event in manifest["events"]:
        result[event["id"]], rows, extra = compare(event, paths)
        if event["id"] == "211211549":
            result["flare"] = check_flare(paths, extra, manifest["flare"])
        bins.extend(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out/"measurements.json").write_text(json.dumps(result, indent=2)+"\n")
    with (args.out/"profiles.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["gbm_id","seconds_from_predicted_arrival","stix_4_10","stix_10_15",
                         "stix_15_25","stix_25_50","stix_50_84","gbm_channels_2_4_rate"])
        writer.writerows(bins)


if __name__ == "__main__":
    main()
