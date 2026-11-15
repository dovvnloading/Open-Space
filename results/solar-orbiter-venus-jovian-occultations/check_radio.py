"""Test both Jupiter contacts with local spectral contrasts and time-shifted controls."""
import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path

import cdflib
import numpy as np

from reproduce import ROOT, write_csv

BANDS = [(4, 8), (8, 16)]


def products(cdf, day):
    midnight = np.datetime64(day)
    def seconds(epoch):
        return (cdflib.cdfepoch.to_datetime(epoch) - midnight) / np.timedelta64(1, "s")

    if day == "2023-08-21":
        for antenna in ("ch0_V1V2", "ch1_V3V4"):
            tag = "hfr_auto_averages_" + antenna
            var = "psp_fld_l2_rfs_" + tag
            frequency = cdf.varget("frequency_" + tag)
            if np.any(frequency != frequency[0]):
                raise ValueError("PSP frequency grid varied")
            valid = (cdf.varget(var + "_gain") == 1) & (cdf.varget(var + "_hl") == 0)
            valid &= cdf.varget(var + "_nsum") == 40
            time = seconds(cdf.varget("epoch_" + tag))[valid]
            spectrum = cdf.varget(var)[valid].astype(float)
            yield antenna[-4:], frequency[0], np.broadcast_to(time[:, None], spectrum.shape), spectrum
    else:
        epoch = seconds(cdf.varget("Epoch"))
        elapsed = cdf.varget("ELAPSED_TIME_SINCE_EPOCH")
        time = epoch[:, None] + elapsed
        spectrum = cdf.varget("PSD_SFU").astype(float)
        config = 10 if day == "2021-08-09" else 9
        good = (cdf.varget("QUALITY_FLAG") == 3) & (cdf.varget("SENSOR_CONFIG") == config)
        keep = good[:, None] & (cdf.varget("INTERPOL_FLAG") == 1)
        spectrum[~keep] = np.nan
        yield "RPW-" + str(config), cdf.varget("FREQUENCY"), time, spectrum


def channel_samples(time, flux, freq, low, high):
    result = []
    for j in np.flatnonzero((freq >= low * 1e6) & (freq < high * 1e6)):
        mask = np.isfinite(flux[:, j]) & (flux[:, j] > 0) & (time[:, j] >= 0) & (time[:, j] < 86400)
        if np.any(mask):
            x = time[mask, j]
            if np.any(np.diff(x) < 0):
                raise ValueError("Sample times are out of order")
            result.append((x, 10 * np.log10(flux[mask, j])))
    return result


def estimate(samples, ingress, egress, width, minimum):
    intervals = [(ingress-width, ingress), (ingress, ingress+width),
                 (egress-width, egress), (egress, egress+width)]
    ingress_contrast, egress_contrast = [], []
    for time, db in samples:
        indices = np.searchsorted(time, intervals, side="left")
        if np.any(np.diff(indices, axis=1) < minimum):
            continue
        pre, inside1, inside2, post = [np.median(db[a:b]) for a, b in indices]
        ingress_contrast.append(pre-inside1)
        egress_contrast.append(post-inside2)
    if len(ingress_contrast) < 3:
        return None
    enter, leave = np.median(ingress_contrast), np.median(egress_contrast)
    return float(enter), float(leave), float(min(enter, leave)), len(ingress_contrast)


def local_jump(samples, event, width, sign):
    # Independent linear fits on the two sides give the fitted step at the contact.
    jumps = []
    for time, db in samples:
        before, edge, after = np.searchsorted(time, [event - width, event, event + width])
        if min(edge-before, after-edge) < 3:
            continue
        x = (time[before:after] - event) / width
        right = (x >= 0).astype(float)
        matrix = np.column_stack((np.ones(len(x)), x, right, x * right))
        fit = np.linalg.lstsq(matrix, db[before:after], rcond=None)[0]
        jumps.append(sign * fit[2])
    return float(np.median(jumps)) if len(jumps) >= 3 else np.nan


def utc_offset(day, value):
    return float((np.datetime64(value.removesuffix("Z")) - np.datetime64(day)) / np.timedelta64(1, "s"))


def run(cache, refresh):
    with (ROOT / "contacts.csv").open() as stream:
        contacts = {r["date"]: r for r in csv.DictReader(stream) if r["direction"] == "jupiter"
                    and r["venus_height_km"] == "0" and r["source_radius_in_target_radii"] == "0"}
    manifest = json.loads((ROOT / "radio_inputs.json").read_text())
    reports, controls, sensitivity = [], [], []
    for day, item in manifest.items():
        path = cache / item["url"].split("/")[-1]
        if refresh and not path.exists():
            cache.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(item["url"], timeout=120) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise ValueError(f"Hash mismatch for {item['url']}")
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Hash mismatch for {path}")
        width, minimum = item["window_s"], item["minimum_samples_per_channel_window"]
        ingress, egress = [utc_offset(day, contacts[day][k]) for k in ("ingress_utc", "egress_utc")]
        for antenna, frequency, time, flux in products(cdflib.CDF(path), day):
            for low, high in BANDS:
                samples = channel_samples(time, flux, frequency, low, high)
                actual = estimate(samples, ingress, egress, width, minimum)
                if actual is None:
                    reports.append([day, antenna, low, high, width, 0, "", "", "", 0, ""])
                    continue
                if day == "2025-02-18":
                    for halfwidth in (5, 10, 20, 30):
                        sensitivity.append([day, low, high, "linear_step", 0, halfwidth,
                            round(local_jump(samples, ingress, halfwidth, -1), 5),
                            round(local_jump(samples, egress, halfwidth, 1), 5)])
                    for shift in range(-60, 61, 10):
                        shifted = estimate(samples, ingress+shift, egress+shift, width, minimum)
                        sensitivity.append([day, low, high, "shifted_contrast", shift, width,
                                            round(shifted[0], 5), round(shifted[1], 5)])
                scores = []
                for shift in range(-86400, 86401, width*4):
                    first, last = ingress+shift, egress+shift
                    if first-width < 0 or last+width > 86400 or not (last+width < ingress-width or first-width > egress+width):
                        continue
                    ctl = estimate(samples, first, last, width, minimum)
                    if ctl is None or ctl[3] != actual[3]:
                        continue
                    controls.append([day, antenna, low, high, width, shift, ctl[3], *[round(n, 5) for n in ctl[:3]]])
                    scores.append(ctl[2])
                rank = (1+sum(s >= actual[2] for s in scores)) / (1+len(scores)) if scores else ""
                reports.append([day, antenna, low, high, width, actual[3],
                                *[round(n, 5) for n in actual[:3]], len(scores), round(rank, 5) if scores else ""])
    write_csv("radio_contacts.csv", ["date","antenna","band_lo_MHz","band_hi_MHz","side_window_s","n_channels",
              "ingress_db","egress_db","smaller_limb_db","n_control_windows","empirical_upper_fraction"], reports)
    write_csv("radio_sensitivity.csv", ["date","band_lo_MHz","band_hi_MHz","estimator","shift_s",
              "side_window_s","ingress_db","egress_db"], sensitivity)
    write_csv("radio_controls.csv", ["date","antenna","band_lo_MHz","band_hi_MHz","side_window_s","shift_s",
              "n_channels","ingress_db","egress_db","smaller_limb_db"], controls)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--download", action="store_true", help="Retrieve any missing pinned CDF")
    args = parser.parse_args()
    run(args.cache, args.download)
