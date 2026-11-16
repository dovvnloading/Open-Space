"""Retrieve cumulative eROSITA catalogue cone and DR1/DR2 aperture upper limits."""

import argparse
import json
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree

from astropy.io.votable import parse_single_table

ROOT = "https://erosita.mpe.mpg.de/erodat"
SURVEYS = ("DR1_eRASS1", "DR2_eRASSc3")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ra", type=float)
    parser.add_argument("dec", type=float)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = {}
    for release in SURVEYS:
        params = dict(ra=args.ra, dec=args.dec, band="024", dr_survey=release)
        url = ROOT + "/upperlimit/service?" + urlencode(params)
        with urlopen(url, timeout=60) as response:
            limit = json.load(response)
        if not limit.get("de_sky"):
            raise ValueError(f"Outside footprint {limit}")
        output[release] = dict(url=url, result=limit)
    params = dict(CAT="DR2_Main", RA=args.ra, DEC=args.dec, SR=30 / 3600.0, VERB=1)
    url = ROOT + "/catalogue/SCS?" + urlencode(params)
    with urlopen(url, timeout=60) as response:
        data = response.read()
    for elem in ElementTree.fromstring(data).findall(".//{*}INFO"):
        if elem.get("name") == "Error" or elem.get("value") in ("ERROR", "OVERFLOW"):
            raise ValueError(elem.text or str(elem.attrib))
    table = parse_single_table(BytesIO(data)).to_table()
    output["DR2_Main"] = dict(
        url=url,
        columns=table.colnames,
        rows=[{name: str(row[name]) for name in table.colnames} for row in table],
    )
    args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
