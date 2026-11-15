"""Refresh input extracts into a new directory. Live orbital solutions and records can change."""

import argparse
import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


def fetch(url, payload=None):
    options = {"headers": {"Content-Type": "application/json"}, "method": "GET"}
    if payload is not None:
        options["data"] = json.dumps(payload).encode()
    with urlopen(Request(url, **options), timeout=90) as response:
        return response.read()


def nordic_rows(page, table_id):
    table = re.search(fr'<figure[^>]*id="{table_id}".*?</figure>', page, re.S).group()
    rows = []
    for row in re.findall(r"<tr.*?</tr>", table, re.S)[1:]:
        cells = re.findall(r"<td.*?>(.*?)</td>", row, re.S)
        rows.append([html.unescape(re.sub("<[^>]+>", "", cell)).strip() for cell in cells])
    return rows


def acquire(dest):
    dest.mkdir(parents=True, exist_ok=False)
    payload = {"desigs": ["240P"], "output_format": ["XML", "OBS80"]}
    raw = json.loads(fetch("https://data.minorplanetcenter.net/api/get-obs", payload))
    # MPC served both components inside the parent query. Column 12 identifies B.
    lines = [s for s in raw[0]["OBS80"].splitlines()
             if s[15:19] in ("2025", "2026") and s[14] in "CBK"]
    (dest / "astrometry.mpc").write_text("\n".join(lines) + "\n")
    headers = ["jd", "x", "y", "z", "vx", "vy", "vz"]
    states = []
    for body in ("90001213", "90001214", "399", "5", "6"):
        params = dict(format="json", COMMAND=f"'{body}'", CENTER="'500@10'", MAKE_EPHEM="'YES'",
                      EPHEM_TYPE="'VECTORS'", START_TIME="'2014-01-01'", STOP_TIME="'2026-07-01'",
                      STEP_SIZE="'4 d'", REF_PLANE="'FRAME'", VEC_TABLE="'2'", CSV_FORMAT="'YES'",
                      OUT_UNITS="'AU-D'", VEC_CORR="'NONE'")
        response = json.loads(fetch("https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)))
        lines = response["result"].split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
        for line in lines:
            cells = line.split(",")
            jd = float(cells[0])
            if 2460895 <= jd <= 2461170:
                states.append([jd, *map(float, cells[2:8]), body])
    pd.DataFrame(states, columns=headers + ["body"]).to_csv(dest / "states.csv", index=False, float_format="%.16g")
    page = fetch("https://arxiv.org/html/2606.21609v1").decode()
    observations = []
    for epoch, offset in zip(nordic_rows(page, "S5.T1"), nordic_rows(page, "S5.T2"), strict=True):
        if epoch[0] == "Apr 7":  # B had no measured photometric detection
            continue
        year = "2025" if int(epoch[2]) <= 365 else "2026"
        times = [datetime.strptime(f"{year} {epoch[0]} {hour}", "%Y %b %d %H:%M")
                 .replace(tzinfo=timezone.utc).timestamp() for hour in epoch[1].split(" - ")]
        observations.append([epoch[0], 2440587.5 + sum(times) / (2 * 86400), *map(float, offset[2:4])])
    pd.DataFrame(observations, columns=["date", "jd", "E", "N"]).to_csv(dest / "not.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_directory", type=Path)
    acquire(parser.parse_args().new_directory)
