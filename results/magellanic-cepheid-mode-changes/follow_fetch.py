"""Fetch extended OGLE-IV curves and EROS cone queries for selected Cepheids."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from urllib.parse import urlencode
from hashlib import sha256
import csv
from fetch import ROOT

TARGETS = ["2332", "1661", "2170", "3278", "3610", "0076"]


def degrees(text):
    a, b, c = map(float, text.split(":"))
    return (-1 if text.startswith("-") else 1) * (abs(a) + b / 60 + c / 3600)


def sources():
    ids = (ROOT / "inputs/smc4_ident.dat").read_text().splitlines()
    for target in TARGETS:
        star = f"OGLE-SMC-CEP-{target}"
        item = next(l.split() for l in ids if l.startswith(star + " "))
        ra = 15 * degrees(item[2])
        dec = degrees(item[3])
        for band in ["I", "V"]:
            yield f"{target}_extended_{band}.dat", f"https://ogledb.astrouw.edu.pl/~ogle/OCVS/data/{band}/{target[-2:]}/{star}.dat"
        for table in ["xcat", "xtime"]:
            pars = {
                "-source": f"II/390/{table}",
                "-c": f"{ra:.7f} {dec:.7f}",
                "-c.rs": "6",
                "-out.max": "30000",
            }
            if table == "xcat":
                pars["-out.all"] = "1"
            else:
                pars["-out"] = "EROS,HJD,Rmag,e_Rmag,Bmag,e_Bmag"
            yield f"{target}_eros_{table}.tsv", "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(
                pars
            )


def get(item):
    name, url = item
    path = ROOT / "inputs" / name
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=100) as response:
                    raw = response.read()
                if b"<html" in raw[:500].lower():
                    raise ValueError(f"HTML at {url}")
                path.write_bytes(raw)
                break
            except OSError:
                if attempt == 2:
                    raise
    raw = path.read_bytes()
    return name, sha256(raw).hexdigest(), url, len(raw)


if __name__ == "__main__":
    with ThreadPoolExecutor(5) as pool:
        rows = list(pool.map(get, sources()))
    with (ROOT / "follow_manifest.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["file", "sha256", "url", "bytes"])
        w.writerows(rows)
