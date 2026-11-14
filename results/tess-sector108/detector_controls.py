"""Same-detector controls for the TIC 257433294 photometric reduction discrepancy."""

import csv
import json
import math
import re
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from urllib.parse import urlencode

import numpy as np
from astropy.io import fits

from fetch import ROOT, get
from measure import measure_events, write_rows


def angular_distance(ra, dec, ra0, dec0):
    a, d, a0, d0 = map(math.radians, [ra, dec, ra0, dec0])
    h = (
        math.sin((d - d0) / 2) ** 2
        + math.cos(d) * math.cos(d0) * math.sin((a - a0) / 2) ** 2
    )
    return math.degrees(2 * math.asin(math.sqrt(h)))


def main():
    params = {
        "-source": "IV/39/tic82",
        "-c": "57.228327 -27.405763",
        "-c.rd": "4",
        "-out": "TIC,RAJ2000,DEJ2000,Tmag",
        "Tmag": "10..15",
        "-out.max": "unlimited",
    }
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(params)
    path = ROOT / "neighbour_field.tsv"
    get(url, path)
    lines = [
        line
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    nearby = list(
        csv.DictReader(StringIO("\n".join(lines[:1] + lines[3:])), delimiter="\t")
    )
    with (ROOT / "ephemerides.csv").open() as inp:
        c = next(r for r in csv.DictReader(inp) if int(r["TIC"]) == 257433294)
    names = re.findall(
        r"-o (tess\d+-s0108-(\d{16})-\d+-s_lc\.fits) ",
        (ROOT / "targets.sh.txt").read_text(),
    )
    by_id = {int(t): n for n, t in names}
    candidates = []
    for r in nearby:
        tic = int(r["TIC"])
        if tic == 257433294 or tic not in by_id:
            continue
        sep = angular_distance(
            float(r["RAJ2000"]), float(r["DEJ2000"]), float(c["ra"]), float(c["dec"])
        )
        candidates.append((sep, tic, float(r["Tmag"]), by_id[tic]))
    # Rank without reference to flux, test the 12 nearest, then restrict to camera 2 CCD 1.
    candidates = sorted(candidates)[:12]

    def fetch_control(row):
        sep, tic, mag, name = row
        padded = f"{tic:016}"
        directory = "/".join(padded[i : i + 4] for i in range(0, 16, 4))
        url = f"https://stpubdata.s3.amazonaws.com/tess/public/tid/s0108/{directory}/{name}"
        path = ROOT / "fits" / name
        digest = get(url, path)
        with fits.open(path) as h:
            camera, ccd = h[0].header["CAMERA"], h[0].header["CCD"]
            tab = h[1].data
            receipt = {
                "tic": tic,
                "Tmag": mag,
                "separation_deg": sep,
                "camera": camera,
                "ccd": ccd,
                "url": url,
                "sha256": digest,
            }
            if (camera, ccd) != (2, 1):
                return receipt, []
            t = np.array(tab["TIME"], float)
            a = np.array(tab["SAP_FLUX"], float)
            b = np.array(tab["PDCSAP_FLUX"], float)
            m = (
                (tab["QUALITY"] == 0)
                & np.isfinite(t)
                & np.isfinite(a)
                & np.isfinite(b)
                & (a > 0)
                & (b > 0)
            )
            events = []
            for label, values in [("SAP", a), ("PDCSAP", b)]:
                ev = measure_events(
                    t[m],
                    values[m],
                    float(c["T0_BJD_TDB"]) - 2457000,
                    float(c["period_d"]),
                    float(c["duration_h"]) / 24,
                    h[1].header["TIMEDEL"],
                    2.0,
                )
                events.extend(
                    {"tic": tic, "extraction": label, "cycle": n, "depth_ppm": v}
                    for n, _, v in ev
                )
        return receipt, events

    with ThreadPoolExecutor(4) as pool:
        results = list(pool.map(fetch_control, candidates))
    (ROOT / "detector_provenance.json").write_text(
        json.dumps({"selection_url": url, "fits": [r for r, _ in results]}, indent=2)
        + "\n"
    )
    write_rows("detector_events.csv", [e for _, ev in results for e in ev])


if __name__ == "__main__":
    main()
