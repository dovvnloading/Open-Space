"""Freeze archive fields used by the period analysis; all raw tables can be re-downloaded."""

import csv
import gzip
import io
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
POSITION = (25.68910759, 13.26281275)


def query(url, params):
    full = url + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(full, timeout=120) as response:
        return response.read().decode()


def reduced(source, output, columns):
    with gzip.GzipFile(filename=str(HERE / output), mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, newline="") as stream:
            writer = csv.DictWriter(
                stream, columns, extrasaction="ignore", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(csv.DictReader(io.StringIO(source)))


if __name__ == "__main__":
    # Discovery selection used the J2000 WDJ position and a 3-arcsec radius.
    ra = 15 * (1 + 42 / 60 + 45.37 / 3600)
    dec = 13 + 15 / 60 + 46.38 / 3600
    ztf = query(
        "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves",
        {"POS": f"CIRCLE {ra:.8f} {dec:.8f} 0.000833", "FORMAT": "CSV"},
    )
    reduced(
        ztf,
        "ztf.csv.gz",
        [
            "oid",
            "expid",
            "hjd",
            "mag",
            "magerr",
            "catflags",
            "filtercode",
            "ra",
            "dec",
            "limitmag",
        ],
    )
    columns = [
        "objID",
        "obsTime",
        "filterID",
        "psfFlux",
        "psfFluxErr",
        "psfQfPerfect",
        "infoFlag3",
    ]
    ps1 = query(
        "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv",
        {"objID": "123910256890265907", "columns": "[" + ",".join(columns) + "]"},
    )
    reduced(ps1, "panstarrs.csv.gz", columns)
