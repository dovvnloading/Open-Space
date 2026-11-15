"""Replay EXO 2030+375 using primary mission products. Python >=3.10, NumPy, Astropy, SciPy."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
from astropy.constants import c
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time
from scipy.interpolate import CubicSpline

ORBITS = (150, 152, 154, 156, 158, 160, 162, 166)
EPOCH, PERIOD = 52756.17, 46.0213
GBM = "https://gammaray.msfc.nasa.gov/gbm/science/pulsars/lightcurves/exo2030.fits.gz"
BAT = "https://swift.gsfc.nasa.gov/results/transients/EXO2030p375.lc.txt"
IDS = {"EXO": "J2032+376", "Cyg X-3": "J2032+409", "Cyg X-1": "J1958+352"}
TARGET, CYG = (308.0633, 37.6375), (308.108, 40.958)
HERE = Path(__file__).resolve().parent


def fetch(cache, url, manifest, optional=False):
    path = cache / url.rsplit("/", 1)[-1]
    if url in manifest:
        return path if manifest[url] else None
    if not path.exists():
        try:
            with urlopen(url, timeout=90) as response:
                body = response.read()
        except HTTPError as exc:
            if optional and exc.code == 404:
                manifest[url] = ""
                return None
            raise
        temporary = path.with_suffix(path.suffix + ".part")
        temporary.write_bytes(body)
        temporary.replace(path)
    manifest[url] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def contrast(x, center, width=2.5, offset=6.0):
    windows = []
    counts = []
    for delta in (-offset, 0, offset):
        cut = (abs(x[:, 0] - center - delta) <= width) & (x[:, 2] > 0)
        counts.append(int(cut.sum()))
        if cut.any():
            t, y, err = x[cut].T
            w = err**-2
            windows.append((np.average(t, weights=w), np.average(y, weights=w), w.sum()**-.5))
    if min(counts) == 0:
        return np.nan, np.nan, counts
    t, y, err = np.asarray(windows).T
    fraction = (t[1] - t[0]) / (t[2] - t[0])
    weights = np.array([fraction-1, 1, -fraction])
    return weights @ y, np.linalg.norm(weights*err), counts


def monitors(cache, manifest):
    hist = fits.getdata(fetch(cache, GBM, manifest), "PSR_HIST")
    hist = hist[hist["DETECTED"] == 1]
    _, indices = np.unique(hist["BARYTIME"], return_index=True)
    hist = hist[indices]
    bat = np.loadtxt(fetch(cache, BAT, manifest))
    bat = bat[(bat[:, 7] == 0) & (bat[:, 2] > 0), :3].copy()
    bat[:, 0] += .5
    data = {"BAT": bat, "GBM": np.column_stack([hist[k] for k in ("BARYTIME", "AMPLITUDE", "AMPSIG")])}
    for name, key in IDS.items():
        url = f"https://maxi.riken.jp/star_data/{key}/{key}_g_lc_1day_all.dat"
        x = np.loadtxt(fetch(cache, url, manifest))
        for column, band in zip((1, 3, 5, 7), ("2-20", "2-4", "4-10", "10-20")):
            data[f"{name} {band}"] = x[:, [0, column, column+1]]
    return hist, data


def angle(ra, dec, pos):
    r, d = np.deg2rad(pos)
    half = np.sin((dec-d)/2)**2 + np.cos(dec)*np.cos(d)*np.sin((ra-r)/2)**2
    return np.rad2deg(2*np.arcsin(np.sqrt(np.clip(half, 0, 1))))


def light_time(center):
    grid = np.arange(center-7, center+8, .5)
    r, d = np.deg2rad(TARGET)
    vector = np.array([np.cos(r)*np.cos(d), np.sin(r)*np.cos(d), np.sin(d)])
    earth = get_body_barycentric("earth", Time(grid, format="mjd", scale="tt"))
    return CubicSpline(grid, earth.xyz.to_value("m").T @ vector / c.value)


def emission_time(local, correction):
    bary = local + correction(local)/86400
    mean = 2*np.pi*(bary-EPOCH)/PERIOD
    e = .41
    anomaly = mean.copy()
    for _ in range(6):
        anomaly -= (anomaly-e*np.sin(anomaly)-mean)/(1-e*np.cos(anomaly))
    omega = np.deg2rad(211.9 + .18*(bary-EPOCH)/365.25)
    delay = 246*((np.cos(anomaly)-e)*np.sin(omega) + np.sqrt(1-e*e)*np.sin(anomaly)*np.cos(omega))
    return bary - delay/86400


def photons(cache, manifest, center):
    accumulated = [[], []]
    correction = light_time(center)
    for day in range(int(center)-4, int(center)+5):
        for tile in (139, 171):
            name = f"mx_mjd{day}_gsc_low_{tile}.evt.gz"
            url = (f"https://heasarc.gsfc.nasa.gov/FTP/maxi/data/obs/MJD{day//1000*1000}/"
                   f"MJD{day}/events/gsc_low/{name}")
            path = fetch(cache, url, manifest, optional=True)
            if path is None:
                continue
            with fits.open(path) as hdus:
                a, h = hdus["EVENTS"].data, hdus["EVENTS"].header
                ok = ((a["PI"] >= 60) & (a["PI"] < 300) & (a["QUAL_POS"] == 0)
                      & (a["QUAL_PI"] == 0) & (a["DFLAG"] == 0))
                ra, dec = (np.deg2rad(a[k][ok].astype(float)) for k in ("RA", "DEC"))
                local = h["MJDREFI"] + h["MJDREFF"] + (a["TIME"][ok]+h.get("TIMEZERO",0))/86400
                t = emission_time(local, correction)
                for output, pos in zip(accumulated, (TARGET, CYG)):
                    output.extend(t[(angle(ra,dec,pos) < 1.25) & (abs(t-center) < 3)])
    return [np.sort(values) for values in accumulated]


def fold(times, center, frequency):
    if not len(times):
        return [0, 0, np.nan, np.nan, np.nan]
    scans = np.split(times, np.flatnonzero(np.diff(times)*86400 > 400)+1)
    h1, h2 = [], []
    for scan in scans:
        v = np.exp(2j*np.pi*(scan-center)*86400*frequency)
        h1.append(v.sum())
        h2.append((v*v).sum())
    h1, h2 = np.array(h1), np.array(h2)
    norm = 2/len(times)
    first, second = norm*abs(sum(h1))**2, norm*abs(sum(h2))**2
    # Rotate both harmonics consistently; preserve the complete pattern within each scan.
    rng = np.random.default_rng(646)
    exceed = 0
    for _ in range(10):
        rot = np.exp(2j*np.pi*rng.random((1000,len(scans))))
        simulated = norm*(abs(rot@h1)**2 + abs((rot*rot)@h2)**2)
        exceed += np.count_nonzero(simulated >= first+second)
    return len(times),len(scans),first,second,(1+exceed)/10001


def write(path, header, rows):
    with path.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Writable cache and output directory")
    parser.add_argument("--monitors-only", action="store_true")
    args = parser.parse_args()
    cache = args.directory
    cache.mkdir(parents=True, exist_ok=True)
    manifest = {}
    hist, data = monitors(cache, manifest)
    monitor_rows, pulse_rows = [], []
    for n in ORBITS:
        epoch = EPOCH + PERIOD*n
        time = hist["BARYTIME"]
        possible = hist[(time > epoch-1) & (time < epoch+12)]
        peak = possible[np.argmax(possible["AMPLITUDE"])]
        center, frequency = peak["BARYTIME"], peak["FREQUENCY"]
        for name, x in data.items():
            excess, error, counts = contrast(x, center)
            monitor_rows.append([n, center, 1/frequency, name, *counts, excess, error, excess/error])
        if not args.monitors_only:
            for region, times in zip(("EXO", "Cyg X-3"), photons(cache, manifest, center)):
                pulse_rows.append([n, region, *fold(times,center,frequency)])
    write(cache/"monitor-replay.csv", ["orbit","center_mjd","period_s","series","left_n","center_n",
          "right_n","excess","formal_error","ratio"], monitor_rows)
    if pulse_rows:
        write(cache/"photon-replay.csv", ["orbit","region","photons","scans","power_1","power_2",
              "scan_rotation_p_two_harmonics"], pulse_rows)
    write(cache/"retrieval.csv", ["url","sha256_or_empty_if_404"], sorted(manifest.items()))


if __name__ == "__main__":
    main()
