"""Geocentric candidate generator. Downloads the fixed RFC release and Horizons grid."""
import csv
import hashlib
import io
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.spatial import cKDTree
from reproduce import EVIDENCE, J2000, direction, sexagesimal
import json
from datetime import timedelta


def main():
    src = json.loads((EVIDENCE/"catalog.json").read_text())["rfc_2026c_cat.txt"]
    with urlopen(src["url"], timeout=120) as response:
        catalog = response.read()
    if hashlib.sha256(catalog).hexdigest() != src["full_sha256"]:
        raise ValueError("Unexpected RFC release content")
    sources = [line.split() for line in catalog.decode().splitlines() if line.startswith("RFC ")]
    ra = [15*sexagesimal(":".join(row[3:6])) for row in sources]
    dec = [sexagesimal(":".join(row[6:9])) for row in sources]
    positions = direction(ra, dec)
    tree = cKDTree(positions)
    params = {"COMMAND":"699","OBJ_DATA":"NO","MAKE_EPHEM":"YES","EPHEM_TYPE":"OBSERVER",
              "CENTER":"500@399","START_TIME":"2026-01-01","STOP_TIME":"2036-01-01",
              "STEP_SIZE":"12h","QUANTITIES":"1,13,20","ANG_FORMAT":"DEG",
              "CSV_FORMAT":"YES","CAL_FORMAT":"JD","EXTRA_PREC":"YES"}
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text&"+urlencode(
        {k:"'"+v+"'" for k,v in params.items()})
    with urlopen(url, timeout=120) as response:
        ephem = response.read().decode()
    rows = csv.reader(io.StringIO(ephem.split("$$SOE")[1].split("$$EOE")[0].strip()))
    jd, ra, dec, diameter = np.array([[float(r[i]) for i in (0,3,4,5)] for r in rows]).T
    saturn = direction(ra, dec)
    hits = {}
    for i in range(1, len(jd)-1):
        velocity = saturn[i+1] - saturn[i-1]   # direction per day
        curvature = 2*(saturn[i+1]+saturn[i-1]-2*saturn[i])
        for j in tree.query_ball_point(saturn[i],0.0027):
            offset = positions[j]-saturn[i]
            dt = (offset @ velocity)/(velocity @ velocity)
            if abs(dt)>.26:
                continue
            miss = np.linalg.norm(offset-dt*velocity-dt*dt*curvature)*206264.806
            if miss < 1.35*diameter[i]+3:
                day = jd[i]+dt
                key = (sources[j][1], int(day+0.5))
                if key not in hits or hits[key][1]>miss:
                    hits[key] = (day,miss)
    for (name,_), (day,miss) in sorted(hits.items(),key=lambda h:h[1]):
        epoch = J2000+timedelta(days=day-2451545)
        print(name, epoch.strftime("%Y-%m-%d %H:%M"), f"{miss:.2f} arcsec")


if __name__ == "__main__":
    main()
