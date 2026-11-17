"""Requery the public DR11 v2 LRS combined catalogue around the prioritized Gaia sources."""
import argparse
import csv
import io
import re
from http.cookiejar import CookieJar
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import HTTPCookieProcessor, Request, build_opener

HERE = Path(__file__).resolve().parent
SEARCH = "https://www.lamost.org/dr11/search"
QUERY = "https://www.lamost.org/dr11/v2.0/table/combined/q"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--medium-resolution", action="store_true")
    args = parser.parse_args()
    search_url = SEARCH
    query_url = QUERY
    table = "combined"
    dest = "lamost_dr11.csv"
    if args.medium_resolution:
        search_url = "https://www.lamost.org/dr11/medcas/search"
        query_url = QUERY.replace("/combined/", "/med_combined/")
        table = "med_combined"
        dest = "lamost_mrs_dr11.csv"
    with (HERE.parent/"radial_followup_gaia.csv").open() as infile:
        positions = [f"{row['ra']},{row['dec']},8" for row in csv.DictReader(infile)]
    session = build_opener(HTTPCookieProcessor(CookieJar()))
    with session.open(search_url, timeout=60) as response:
        form = response.read().decode()
    token = re.search(r'name="csrf_token" value="([^"]+)"', form)
    if not token:
        raise ValueError("The public query form changed")
    fields = ["obsid","designation","gaia_source_id","ra","dec","z","rv","rv_err",
              "class","subclass","teff","logg","feh","snrg"]
    if args.medium_resolution:
        fields = ["obsid","designation","gaia_source_id","ra","dec","rv_b0","rv_b0_err","rv_r0","rv_r0_err"]
    params = {
        "csrf_token": token[1],
        "output.fmt": "csv",
        "pos.type": "proximity",
        "pos.radecTextarea": "\n".join(positions),
        "pos.defaultRadius": "8",
        "showcol": [table + "." + name for name in fields],
    }
    request = Request(query_url, data=urlencode(params,doseq=True).encode(),headers={"Referer":search_url})
    with session.open(request, timeout=90) as response:
        data = response.read()
    header = next(csv.reader(io.StringIO(data.decode("utf-8-sig"))), [])
    if not {"obsid","ra","dec","gaia_source_id"}.issubset(header):
        raise ValueError("The returned file has an unexpected schema")
    (HERE/dest).write_bytes(data)


if __name__ == "__main__":
    main()
