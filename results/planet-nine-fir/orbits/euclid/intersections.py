"""Cross each archive footprint at its UTC midpoint with the conditioned Phan family."""
import argparse
import csv
import hashlib
import io
import json
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.time import Time
from scipy.interpolate import CubicSpline
from matplotlib.path import Path as Polygon

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import propagate as p


def fetch(url, path):
    data = urlopen(url, timeout=90).read()
    path.write_bytes(data)
    return dict(url=url, bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def spacecraft():
    qs = dict(format="json", COMMAND="'-680'", CENTER="'500@0'", MAKE_EPHEM="'YES'",
              EPHEM_TYPE="'VECTORS'", START_TIME="'2024-06-25'", STOP_TIME="'2025-01-05'",
              STEP_SIZE="'1 d'", OUT_UNITS="'AU-D'", REF_PLANE="'FRAME'",
              VEC_TABLE="'1'", CSV_FORMAT="'YES'", VEC_CORR="'NONE'")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(qs)
    dest = HERE/"horizons.json"
    if not dest.exists():
        receipt = fetch(url, dest)
        (HERE/"horizons_receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    data = json.loads(dest.read_bytes())["result"].split("$$SOE")[1].split("$$EOE")[0].strip()
    arr = np.array([[float(row[k]) for k in (0,2,3,4)] for row in csv.reader(io.StringIO(data))])
    return CubicSpline(arr[:,0]-2451545, arr[:,1:], extrapolate=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata", type=Path, help="euclid_dated.response from the archive queries")
    args = parser.parse_args()
    sc = spacecraft()
    models = p.family(json.loads((HERE.parent/"orbit_input.json").read_text()))
    table = []
    covered = np.zeros(len(models[0]), dtype=bool)
    cache = {}
    raw=args.metadata.read_bytes()
    if args.metadata.suffix=='.gz':
        import gzip
        raw=gzip.decompress(raw)
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    for row in rows:
        t = Time((float(row["t_min"])+float(row["t_max"]))/2, format="mjd", scale="utc")
        et = t.tt.jd-2451545.
        if et not in cache:
            cache[et] = np.array(p.sky(models, et, sc(t.tdb.jd-2451545))).T
        points = cache[et]
        corners = np.array([float(x) for x in row["s_region"].split()[1:]]).reshape(-1,2)
        # Project to the tangent plane so polygon edges follow great circles.
        n = p.unit(float(row["s_ra"]),float(row["s_dec"]))
        ra0 = np.radians(float(row["s_ra"]))
        east = np.array([-np.sin(ra0),np.cos(ra0),0])
        north = np.cross(n,east)
        def tangent(positions):
            u = p.unit(positions[:,0], positions[:,1]).T
            return np.array([(u@east)/(u@n),(u@north)/(u@n)]).T
        inside = Polygon(tangent(corners)).contains_points(tangent(points))
        covered |= inside
        table.append([row["obs_id"],t.isot,row["access_url"].split("file_name=")[-1],
                      inside.sum(),f'{inside.mean():.6f}'])
    with (HERE/"coverage.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["obs_id","utc_mid","product","tracks_in_polygon","grid_fraction"])
        w.writerows(table)
    info = dict(archive_rows=len(rows),positive_rows=int(sum(x[3]>0 for x in table)),
                geomet_matched=int(covered.sum()),grid_size=int(covered.size),
                metadata_sha256=hashlib.sha256(raw).hexdigest())
    print(info)
    (HERE/"summary.json").write_text(json.dumps(info,indent=2)+"\n")


if __name__ == "__main__":
    main()
