"""Decode the frozen NOAA and NASA products into common UTC-minute inputs."""

import argparse
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path

import astropy.units as units
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
import cdflib
import erfa
import netCDF4
import numpy as np
from scipy.interpolate import BarycentricInterpolator

START = np.datetime64("2026-11-05")
MINUTES = 1440
SECONDS = np.arange(3 * MINUTES) * 60.0 + 30.0
ROOT = Path(__file__).resolve().parent


def verify_sources(data_root):
    sources = json.loads((ROOT / "sources.json").read_text(encoding="utf-8"))
    for source in sources:
        path = data_root / source["mission"] / source["file"]
        content = path.read_bytes()
        if len(content) != source["bytes"] or hashlib.sha256(content).hexdigest() != source["sha256"]:
            raise ValueError(f"Unexpected size or checksum for {path}")
    return sources


def utc_seconds(epochs):
    utc = cdflib.cdfepoch.to_datetime(epochs)
    return (utc - START) / np.timedelta64(1, "s")


def minute_mean(times, vectors, good, day, minimum):
    seconds = times - day * 86400
    in_day = (seconds >= 0) & (seconds < 86400)
    if not np.all(in_day) or np.any(np.diff(times) <= 0):
        raise ValueError("Expected strictly increasing time labels within a UTC day")
    bins = (seconds // 60).astype(int)
    selected = good & np.all(np.isfinite(vectors) & (np.abs(vectors) < 1000), axis=1)
    count = np.bincount(bins[selected], minlength=MINUTES)
    sums = np.array([
        np.bincount(bins[selected], weights=component[selected], minlength=MINUTES)
        for component in vectors.T
    ]).T
    result = np.full((MINUTES, 3), np.nan)
    np.divide(sums, count[:, None], out=result, where=count[:, None] >= minimum)
    return result, count


def as_float(variable):
    return np.ma.asarray(variable[:], dtype=float).filled(np.nan)


def solar_magnetic(file, day):
    epoch = np.datetime64("1958-01-01")
    start_seconds = (START - epoch) / np.timedelta64(1, "s")
    with netCDF4.Dataset(file) as dataset:
        times = as_float(dataset["time_sec"]) / 1e6 - start_seconds
        vectors = as_float(dataset["b_gse_sec"])
        flags = as_float(dataset["flags_summary"])
    return minute_mean(times, vectors, flags == 0, day, minimum=45)


def wind_magnetic(file, day):
    cdf = cdflib.CDF(file)
    return minute_mean(
        utc_seconds(cdf.varget("Epoch3")), cdf.varget("B3GSE"), cdf.varget("NUM3_PTS") > 0,
        day, minimum=15,
    )


def wind_plasma(file):
    cdf = cdflib.CDF(file)
    times = utc_seconds(cdf.varget("Epoch"))
    vectors = cdf.varget("V_GSE")
    good = (
        (cdf.varget("DQF") == 0)
        & np.all(np.isfinite(vectors) & (np.abs(vectors) <= 1500), axis=1)
        & (vectors[:, 0] >= -1000)
        & (vectors[:, 0] <= -200)
    )
    limited = (cdf.varget("QF_V").astype(int) & 8192) != 0
    return times[good], vectors[good], limited[good]


def local_velocity(times, vectors):
    left = np.searchsorted(times, SECONDS - 900, side="left")
    right = np.searchsorted(times, SECONDS + 900, side="right")
    velocity = np.full((len(SECONDS), 3), np.nan)
    for index in np.flatnonzero(right - left >= 8):
        velocity[index] = np.median(vectors[left[index]:right[index]], axis=0)
    return velocity, right - left


def gse_basis(time):
    sun = (get_body_barycentric("sun", time) - get_body_barycentric("earth", time))
    x = sun.xyz.to_value(units.km).T
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    pole = erfa.ecm06(time.tt.jd1, time.tt.jd2)[..., 2, :]
    y = np.cross(pole, x)
    y /= np.linalg.norm(y, axis=1, keepdims=True)
    return np.stack((x, y, np.cross(x, y)), axis=1)


def solar_orbit(file):
    rows = [line.split() for line in file.read_text(encoding="ascii").splitlines()
            if line.startswith("2026-")]
    seconds = (np.array([row[0] for row in rows], dtype="datetime64[ms]") - START)
    seconds = seconds / np.timedelta64(1, "s")
    positions = np.array([row[1:4] for row in rows], dtype=float)
    start = np.clip(np.searchsorted(seconds, SECONDS) - 4, 0, len(seconds) - 8)
    if np.any(SECONDS < seconds[0]) or np.any(SECONDS > seconds[-1]):
        raise ValueError("SOLAR-1 orbit cannot be extrapolated")
    interpolated = np.empty((len(SECONDS), 3))
    for index in np.unique(start):
        subset = slice(index, index + 8)
        target = start == index
        interpolated[target] = BarycentricInterpolator(
            seconds[subset], positions[subset], axis=0
        )(SECONDS[target])
    utc = START + SECONDS.astype("timedelta64[s]")
    return np.einsum("nij,nj->ni", gse_basis(Time(utc)), interpolated)


def wind_orbit(files):
    seconds, gse, gci = [], [], []
    for file in files:
        cdf = cdflib.CDF(file)
        seconds.extend(utc_seconds(cdf.varget("Epoch")))
        gse.extend(cdf.varget("GSE_POS"))
        gci.extend(cdf.varget("GCI_POS"))
    seconds, gse, gci = np.asarray(seconds), np.asarray(gse), np.asarray(gci)
    utc = START + seconds.astype("timedelta64[s]")
    fit = np.einsum("nij,nj->ni", gse_basis(Time(utc)), gci)
    difference = np.max(np.abs(fit - gse), axis=0)
    if np.any(difference >= 2):
        raise ValueError("Calculated axes fail the Wind GSE check")
    right = np.clip(np.searchsorted(seconds, SECONDS), 1, len(seconds) - 1)
    covered = (
        (seconds[right - 1] <= SECONDS) & (SECONDS <= seconds[right])
        & (seconds[right] - seconds[right - 1] <= 600)
    )
    interpolated = np.stack(
        [np.interp(SECONDS, seconds, component) for component in gse.T], axis=1
    )
    interpolated[~covered] = np.nan
    return interpolated, difference


def run(data_root, output):
    sources = verify_sources(data_root)
    solar_files, wind_files = data_root / "solar1", data_root / "wind"
    solar_b, solar_n, wind_b, wind_n = [], [], [], []
    times, vectors, limited, orbit_files = [], [], [], []
    for day in range(3):
        date = (datetime(2026, 11, 5) + timedelta(days=day)).strftime("%Y%m%d")
        solar_file = next(solar_files / entry["file"] for entry in sources
                          if entry.get("date") == date and entry["mission"] == "solar1")
        field, count = solar_magnetic(solar_file, day)
        solar_b.append(field)
        solar_n.append(count)
        field, count = wind_magnetic(wind_files / f"wi_h0_mfi_{date}_v03.cdf", day)
        wind_b.append(field)
        wind_n.append(count)
        time, vector, limited_flag = wind_plasma(wind_files / f"wi_k0_swe_{date}_v01.cdf")
        times.extend(time)
        vectors.extend(vector)
        limited.append([len(time), int(np.count_nonzero(limited_flag))])
        orbit_files.append(wind_files / f"wi_or_pre_{date}_v10.cdf")
    velocity, velocity_count = local_velocity(np.asarray(times), np.asarray(vectors))
    wind_position, frame_check = wind_orbit(orbit_files)
    solar_position = solar_orbit(solar_files / "SWFO_Ephem_28day_2026308_20263081704.oem")
    np.savez_compressed(
        output, solar=np.vstack(solar_b), wind=np.vstack(wind_b), solar_count=np.hstack(solar_n),
        wind_count=np.hstack(wind_n), solar_position=solar_position, wind_position=wind_position,
        flow=velocity, flow_count=velocity_count, plasma_counts=np.asarray(limited),
        wind_gse_check_km=frame_check,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.data_root, arguments.output)
