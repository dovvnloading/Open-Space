"""Compare reception times from Horizons observer tables with the light-ray replay."""
import argparse
import csv
import gzip
import json
from datetime import datetime, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
START = datetime(2050, 3, 12, 3, 30)
DATA = ROOT / "inputs" / "observer.csv.gz"


def update():
    with gzip.open(DATA, "wt", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["body", "time_tt", "astrometric_ra_deg", "astrometric_dec_deg", "range_au"])
        for body in (701, 702):
            params = dict(format="json", COMMAND=f"'{body}'", CENTER="'500@399'",
                          MAKE_EPHEM="'YES'", EPHEM_TYPE="'OBSERVER'", START_TIME="'2050-03-12 03:30'",
                          STOP_TIME="'2050-03-13 02:40'", STEP_SIZE="'1 m'", QUANTITIES="'1,20'",
                          CSV_FORMAT="'YES'", ANG_FORMAT="'DEG'", EXTRA_PREC="'YES'",
                          TIME_TYPE="'TT'", OBJ_DATA="'NO'")
            url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
            with urlopen(url, timeout=90) as reply:
                raw = json.load(reply)["result"].split("$$SOE")[1].split("$$EOE")[0]
            for row in csv.reader(raw.strip().splitlines()):
                writer.writerow([body, row[0].strip(), *[item.strip() for item in row[3:6]]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    if parser.parse_args().refresh:
        update()
    with gzip.open(DATA, "rt") as f:
        rows = list(csv.DictReader(f))
    directions, ranges = {}, {}
    for body in ("701", "702"):
        selected = [r for r in rows if r["body"] == body]
        ra, dec, distance = np.array([[float(row[k]) for k in
                                     ("astrometric_ra_deg", "astrometric_dec_deg", "range_au")]
                                    for row in selected]).T
        minutes = np.array([(datetime.strptime(r["time_tt"], "%Y-%b-%d %H:%M") - START
                             ).total_seconds()/60 for r in selected])
        ra, dec = np.radians(ra), np.radians(dec)
        directions[body] = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], axis=1)
        ranges[body] = distance * 149597870.7
    separation = np.linalg.norm(np.cross(directions["701"], directions["702"]), axis=1) * ranges["701"]
    curve = CubicSpline(minutes, separation)
    predictions = [r for r in json.loads((ROOT / "events.json").read_text())["ura184"]
                   if r["kind"] == "occultation"]
    results = []
    for event in predictions:
        expected = datetime.fromisoformat(event["minimum_received_tdb"])
        center = (expected - START).total_seconds()/60
        result = minimize_scalar(curve, bounds=(center-30, center+30), method="bounded")
        time = START + timedelta(minutes=result.x)
        results.append(dict(observer_minimum_tt=time.isoformat(timespec="milliseconds"),
                            observer_impact_km=round(float(result.fun), 2),
                            observer_tt_minus_replay_tdb_s=round((time-expected).total_seconds(), 2)))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
