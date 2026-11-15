"""Reproduce the archived GRB subset, Mars timing scans and sky-band intersections."""
import argparse
import csv
import hashlib
from datetime import date, datetime, timedelta
from pathlib import Path
import struct
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import correlate

ROOT = Path(__file__).resolve().parent
EVENTS = {"bn250902062": "GRB 250902A", "bn260307478": "GRB 260307A"}
STEP = 0.25
C_KM_S = 299792.458


def read_csv(name):
    with (ROOT / name).open() as file:
        return list(csv.DictReader(file))


def save_csv(path, rows):
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def obtain(cache, raw_manifest, fname):
    row = raw_manifest[fname]
    path = cache / fname
    if not path.exists():
        with urlopen(row["url"], timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Archive hash changed for {fname}")
        path.write_bytes(data)
    elif hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise ValueError(f"Unexpected cached bytes for {fname}")
    return path


def mars_series(path):
    """Use type-3, good-flag frames. T1 is 1/1024 s, T2 is modulo 65536 quarter seconds."""
    data = path.read_bytes()
    if len(data) % 852:
        raise ValueError(f"Incomplete HEND record in {path}")
    t, counts = [], []
    for start in range(0, len(data), 852):
        record = data[start:start + 852]
        if record[283] != 3 or any(record[273:277]):
            continue
        utc = datetime.fromisoformat(record[24:47].decode())
        seconds = (utc - utc.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
        sync, trigger = struct.unpack(">2I", record[306:314])
        phase = ((256 * trigger - sync + 2**23) % 2**24 - 2**23) / 1024
        t.extend(seconds - 30 + (np.arange(240) + 0.5) * STEP + phase - 39.75)
        counts.extend(record[576:816])
    i = np.rint(np.asarray(t) / STEP).astype(int)
    keep = (i >= 0) & (i < 345600)
    n = np.bincount(i[keep], minlength=345600)
    sums = np.bincount(i[keep], weights=np.asarray(counts)[keep], minlength=345600)
    rate = sums / np.maximum(n, 1)
    good = n > 0
    short = np.nanmedian(np.where(good, rate, np.nan).reshape(-1, 40), axis=1)
    valid = np.isfinite(short)
    short = np.interp(np.arange(len(short)), np.flatnonzero(valid), short[valid])
    background = np.interp(np.arange(len(n)), 40 * np.arange(len(short)) + 19.5,
                           median_filter(short, size=21))
    series = (rate - background) * np.sqrt(n / np.maximum(background, 1))
    series[~good] = 0
    return series, good


def fermi_profile(row, paths, channels):
    duration = float(row["t90_s"])
    grid = np.arange(-20, max(20, duration) + 25, STEP)
    profile = np.zeros_like(grid)
    for path in paths:
        with fits.open(path) as product:
            spectrum = product["SPECTRUM"].data
            t = (spectrum["TIME"] + spectrum["ENDTIME"]) / 2 - product[0].header["TRIGTIME"]
            exposure = spectrum["EXPOSURE"]
            good = (spectrum["QUALITY"] == 0) & (exposure > .001)
            bg = good & (((t > -140) & (t < -30)) | ((t > duration + 50) & (t < duration + 150)))
            counts = spectrum["COUNTS"][:, channels].sum(1)
            baseline = counts[bg].sum() / exposure[bg].sum()
            inside = good & (t >= grid[0] - STEP / 2) & (t < grid[-1] + STEP / 2)
            i = np.floor((t[inside] - grid[0]) / STEP + .5).astype(int)
            profile += np.bincount(i, weights=counts[inside] - exposure[inside] * baseline,
                                   minlength=len(grid))
    profile = gaussian_filter1d(profile, 1.5)
    tail = profile[(grid < -6) | (grid > duration + 12)]
    noise = 1.4826 * np.median(np.abs(tail - np.median(tail)))
    profile[profile < max(2 * noise, profile.max() * .02)] = 0
    return grid, profile


def direction(ra, dec):
    lon, lat = np.deg2rad(ra), np.deg2rad(dec)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def scores(series, good, profile):
    score = correlate(series, profile, "valid", method="fft") / np.linalg.norm(profile)
    covered = correlate(good.astype(float), profile**2, "valid", method="fft") / sum(profile**2)
    score[covered <= .90] = np.nan
    return score


def controls(score, band, include):
    width = int(np.ceil(2 * band / STEP))
    size = len(score) // width * width
    blocks = score[:size].reshape(-1, width)
    valid = (np.isfinite(blocks).sum(1) >= .90 * width) & include[:size].reshape(-1, width).all(1)
    return np.nanmax(blocks[valid], axis=1)


def measure(score, lags, predicted, halfwidth):
    ids = np.flatnonzero((np.abs(lags - predicted) < halfwidth) & np.isfinite(score))
    k = ids[np.argmax(score[ids])]
    return float(lags[k]), float(score[k])


def sky_intersection(path, baseline, delay, widths):
    with fits.open(path) as file:
        tab = file[1]
        if tab.header["ORDERING"] != "NESTED":
            raise ValueError("Expected nested sky map")
        prob = np.asarray(tab.data["PROBABILITY"], dtype=float).ravel()
        ns = tab.header["NSIDE"]
    native = HEALPix(nside=ns, order="nested")
    subgrid = HEALPix(nside=ns * 8, order="nested")
    ranked = np.argsort(-prob)
    size = np.searchsorted(np.cumsum(prob[ranked]), .9 * prob.sum()) + 1
    children = (ranked[:size, None] * 64 + np.arange(64)).ravel()
    lon, lat = subgrid.healpix_to_lonlat(children)
    predicted = -(baseline @ direction(lon.deg, lat.deg))
    native_area = size * native.pixel_area.to_value(u.deg**2)
    return [(width, native_area, np.sum(np.abs(predicted - delay) <= width)
             * subgrid.pixel_area.to_value(u.deg**2)) for width in widths]


def run(args):
    manifest = {r["filename"]: r for r in read_csv("manifest.csv")}
    selected = read_csv("cohort.csv")
    rows, nulls, bands, annuli, sensitivity, curves, checks = [], [], [], [], [], [], []
    references = {r["trigger"]: r for r in read_csv("reference_positions.csv")}
    args.cache.mkdir(exist_ok=True, parents=True)
    args.out.mkdir(exist_ok=True, parents=True)

    for row in selected:
        name = row["trigger"]
        if args.targets_only and name not in EVENTS:
            continue
        day = date.fromisoformat(f"20{name[2:4]}-{name[4:6]}-{name[6:8]}")
        files = [obtain(args.cache, manifest, k) for k in manifest if k.startswith("glg_ctime_") and name in k]
        series, good = mars_series(obtain(args.cache, manifest, day.strftime("%Y%m%d") + ".dat"))
        grid, profile = fermi_profile(row, files, [3, 4, 5, 6])
        score = scores(series, good, profile)
        t0 = float(row["mjd"]) % 1 * 86400
        lags = np.arange(len(score)) * STEP - grid[0] - t0
        baseline = np.array([float(row[f"{axis}_km"]) for axis in "xyz"]) / C_KM_S
        predicted = -baseline @ direction(float(row["ra"]), float(row["dec"]))
        length = np.linalg.norm(baseline)
        band = 80 if row["group"] == "calibration" else 5 + 3 * np.deg2rad(np.hypot(3.7, float(row["err_deg"]))) * length
        delay, height = measure(score, lags, predicted, band)
        ctrl = controls(score, band, np.abs(lags - predicted) > 2 * band)
        rows.append(dict(trigger=name, group=row["group"], catalog_lag_s=predicted, peak_lag_s=delay,
                         residual_s=delay-predicted, score=height, half_search_s=band,
                         off_blocks=len(ctrl), off_exceed=int(sum(ctrl >= height)), off_max=max(ctrl)))
        if name in references:
            reference = references[name]
            refdelay = -baseline @ direction(float(reference["ra_deg"]), float(reference["dec_deg"]))
            cross = abs(np.arccos(-refdelay / length) - np.arccos(-delay / length))
            checks.append(dict(trigger=name, reference=reference["reference"], predicted_s=refdelay,
                               fitted_minus_ref_s=delay-refdelay, cross_track_deg=np.rad2deg(cross)))
        if name not in EVENTS:
            continue
        # This sign gives angular radius about the Earth-to-Mars vector (not its antipode).
        center = baseline / length
        annuli.append(dict(trigger=name, center_ra_deg=np.rad2deg(np.arctan2(center[1], center[0])) % 360,
                           center_dec_deg=np.rad2deg(np.arcsin(center[2])), baseline_light_s=length,
                           delay_mars_minus_earth_s=delay, radius_deg=np.rad2deg(np.arccos(-delay/length))))
        for channel in ([2, 3], [3, 4], [4, 5, 6], [3, 4, 5, 6]):
            _, alternate = fermi_profile(row, files, channel)
            altlag, altpeak = measure(scores(series, good, alternate), lags, predicted, band)
            sensitivity.append(dict(trigger=name, channels=",".join(map(str,channel)), lag_s=altlag, score=altpeak))
        sky = next(k for k in manifest if k.startswith("glg_healpix") and name in k)
        for width, base, sliced in sky_intersection(obtain(args.cache, manifest, sky), baseline, delay, (1, 2, 5, 24)):
            bands.append(dict(trigger=name, half_width_s=width, gbm_90_deg2=base, intersection_deg2=sliced))
        indices = np.rint((t0 + delay + grid) / STEP).astype(int)
        for t, f, m, ok in zip(grid, profile, series[indices], good[indices]):
            curves.append(dict(trigger=name, seconds_from_gbm=t, fermi_template=f, mars_standardized=m, covered=int(ok)))
        for offset in range(-8, 9):
            if not offset:
                continue
            d = (day + timedelta(days=offset)).strftime("%Y%m%d")
            y, g = mars_series(obtain(args.cache, manifest, d + ".dat"))
            s = scores(y, g, profile)
            blocks = controls(s, band, np.ones(len(s), dtype=bool))
            nulls.append(dict(trigger=name, day=d, n=len(blocks), exceed=int(sum(blocks >= height)),
                              max_score=max(blocks)))
    for fname, table in [("measurements",rows), ("controls",nulls), ("areas",bands), ("annuli",annuli),
                         ("channels",sensitivity), ("profiles",curves), ("reference_checks",checks)]:
        save_csv(args.out / (fname + ".csv"), table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--targets-only", action="store_true")
    run(parser.parse_args())
