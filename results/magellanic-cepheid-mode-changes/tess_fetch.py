"""Fetch the seven public TESS visits covering OGLE-SMC-CEP-2332."""

import csv, hashlib, io
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile
from screen import ROOT, INPUT

RA, DEC = 13.762583333333, -73.233083333333
BASE = "https://mast.stsci.edu/tesscut/api/v0.1/"


def acquire(s):
    url = (
        BASE
        + "astrocut?"
        + urlencode(dict(ra=RA, dec=DEC, sector=s, x=11, y=11, units="px"))
    )
    p = INPUT / "tess" / f"s{s}.fits"
    if not p.exists():
        for i in range(3):
            try:
                payload = urlopen(url, timeout=150).read()
                with ZipFile(io.BytesIO(payload)) as archive:
                    (member,) = [n for n in archive.namelist() if n.endswith(".fits")]
                    p.write_bytes(archive.read(member))
                break
            except OSError:
                if i == 2:
                    raise
    return s, url, p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()


if __name__ == "__main__":
    (INPUT / "tess").mkdir(parents=True, exist_ok=True)
    # Frozen visits in the sector inventory at time of analysis.
    with ThreadPoolExecutor(3) as pool:
        rows = list(pool.map(acquire, [1, 27, 28, 67, 68, 94, 95]))
    with (ROOT / "tess_manifest.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["sector", "url", "bytes", "fits_sha256"])
        w.writerows(rows)
