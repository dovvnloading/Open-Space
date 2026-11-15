"""Build one-minute receiver power from public Parker Solar Probe CDFs and JPL geometry."""
import argparse
import csv
import hashlib
import gzip
import json
from datetime import date, datetime, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import cdflib
import numpy as np

HERE = Path(__file__).resolve().parent
DATES = ["2018-10-03", "2019-12-26", "2020-07-11", "2021-02-20",
         "2021-10-16", "2023-08-21", "2024-11-06"]
BANDS = [(1.3, 2), (2, 4), (4, 8), (8, 16)]
ANTENNAS = ["ch0_V1V2", "ch1_V3V4"]


def horizons(day):
    query = dict(format="json", COMMAND="'299'", CENTER="'500@-96'",
                 EPHEM_TYPE="'OBSERVER'", START_TIME=f"'{day}'",
                 STOP_TIME=f"'{date.fromisoformat(day) + timedelta(days=1)}'",
                 STEP_SIZE="'1 m'", QUANTITIES="'13,20'", CSV_FORMAT="'YES'")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(query)
    with urlopen(url, timeout=90) as response:
        receipt = json.load(response)["result"]
    data = receipt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    times, diameters = [], []
    for line in csv.reader(data):
        times.append(datetime.strptime(line[0].strip(), "%Y-%b-%d %H:%M"))
        diameters.append(float(line[3]))
    return url, times, np.array(diameters) / 7200


def median_spectra(times, values, channels):
    # Equal weight for each one-minute bin, irrespective of survey cadence.
    order = np.argsort(times)
    times, values = times[order], values[order]
    positions = np.searchsorted(times, np.arange(1441))
    means = np.full((1440, len(channels)), np.nan)
    counts = np.diff(positions)
    for i in np.flatnonzero(counts >= 3):
        spec = np.nanmedian(values[positions[i]:positions[i+1]], axis=0)
        for j, use in enumerate(channels):
            if np.isfinite(spec[use]).sum() >= use.sum() / 2:
                means[i, j] = np.nanmedian(10 * np.log10(spec[use]))
    return means, counts


def extract(path, day, antenna):
    cdf = cdflib.CDF(path)
    key = "hfr_auto_averages_" + antenna
    var = "psp_fld_l2_rfs_" + key
    flux = cdf.varget(var).astype(float)
    freq = cdf.varget("frequency_" + key)
    if not np.all(freq == freq[0]):
        raise ValueError("Frequency grid changed within the day")
    times = cdflib.cdfepoch.to_datetime(cdf.varget("epoch_" + key))
    minutes = (times - np.datetime64(day)) / np.timedelta64(1, "m")
    quality = (cdf.varget(var + "_gain") == 1) & (cdf.varget(var + "_hl") == 0)
    nsum = cdf.varget(var + "_nsum")
    flux[(flux <= 0) | (flux > float(cdf.varattsget(var)["VALIDMAX"]))] = np.nan
    channels = [(freq[0] >= lo*1e6) & (freq[0] < hi*1e6) for lo, hi in BANDS]
    full = {}
    for setting in [40, 80]:
        good = quality & (nsum == setting)
        if good.any():
            full[setting] = median_spectra(minutes[good], flux[good], channels)
    return full


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path, help="Directory for the seven original daily CDFs")
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    manifest_path = HERE/"manifest.json"
    expected = {r["date"]:r["sha256"] for r in json.loads(manifest_path.read_text())} if manifest_path.exists() else {}
    metadata = []
    columns = ["date", "utc_minute", "venus_radius_deg", "antenna", "nsum",
               "spectra", *[f"db_{lo}_{hi}" for lo, hi in BANDS]]
    with gzip.open(HERE/"minutes.csv.gz", "wt", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(columns)
        for day in DATES:
            geo_url, times, radius = horizons(day)
            minute = int(np.argmax(radius))
            filename = "psp_fld_l2_rfs_hfr_" + day.replace("-", "") + "_v03.cdf"
            url = f"https://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/rfs_hfr/{day[:4]}/{filename}"
            path = args.cache/filename
            if not path.exists():
                with urlopen(url, timeout=120) as response:
                    path.write_bytes(response.read())
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if day in expected and digest != expected[day]:
                raise ValueError(f"Changed CDF bytes for {day}")
            metadata.append(dict(date=day, maximum_radius_utc_minute=minute,
                                 maximum_radius_deg=round(float(radius[minute]), 6),
                                 cdf=url, sha256=digest,
                                 geometry=geo_url))
            for ant in ANTENNAS:
                for nsum, (power, count) in extract(path, day, ant).items():
                    for i in np.flatnonzero(count):
                        writer.writerow([day, i, round(radius[i], 6), ant, nsum, count[i],
                                         *[round(v, 5) if np.isfinite(v) else "" for v in power[i]]])
    (HERE/"manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
