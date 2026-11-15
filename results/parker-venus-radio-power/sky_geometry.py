"""Rebuild apparent Venus and source positions at 15-second cadence near each flyby."""
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent


def sample(args):
    day, center, body = args
    date = datetime.fromisoformat(day) + timedelta(minutes=center)
    start, stop = (date + timedelta(minutes=x) for x in [-60,60])
    query = dict(format="json", COMMAND=f"'{body}'", CENTER="'500@-96'", EPHEM_TYPE="'OBSERVER'",
                 START_TIME=f"'{start:%Y-%m-%d %H:%M}'", STOP_TIME=f"'{stop:%Y-%m-%d %H:%M}'",
                 STEP_SIZE="'480'", QUANTITIES="'2,13'", CSV_FORMAT="'YES'",
                 ANG_FORMAT="'DEG'", EXTRA_PREC="'YES'", TIME_DIGITS="'SECONDS'")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(query)
    with urlopen(url, timeout=90) as response:
        report = json.load(response)["result"]
    raw = report.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    stamp, xyz, radius = [], [], []
    for line in csv.reader(raw):
        ra, dec = np.radians([float(line[3]), float(line[4])])
        stamp.append(datetime.strptime(line[0].strip().split(".")[0], "%Y-%b-%d %H:%M:%S"))
        xyz.append([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
        radius.append(float(line[5]) / 7200)
    if len(stamp) != 481 or stamp[1] - stamp[0] != timedelta(seconds=15):
        raise ValueError("Unexpected Horizons grid")
    return (day, body), (stamp, np.array(xyz), np.array(radius))


def main():
    manifest = json.loads((HERE/"manifest.json").read_text())
    tasks = [(r["date"], r["maximum_radius_utc_minute"], b) for r in manifest for b in [299,10,399,599]]
    with ThreadPoolExecutor(4) as workers:
        data = dict(workers.map(sample,tasks))
    with (HERE/"sightlines.csv").open("w",newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["utc","venus_radius_deg","solar_center_clearance_deg","earth_center_clearance_deg","jupiter_center_clearance_deg"])
        for r in manifest:
            day = r["date"]
            times, venus, rad = data[(day,299)]
            cols = []
            for body in [10,399,599]:
                stamps, vec, _ = data[(day,body)]
                if stamps != times:
                    raise ValueError("Horizons epochs differ")
                sep = np.degrees(np.arctan2(np.linalg.norm(np.cross(venus,vec),axis=1),(venus*vec).sum(axis=1)))
                cols.append(sep-rad)
            for i,t in enumerate(times):
                writer.writerow([t.isoformat(),round(rad[i],5),*[round(v[i],5) for v in cols]])


if __name__ == "__main__":
    main()
