"""Requery the public DR11 v2 LRS combined catalogue around the prioritized Gaia sources."""
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
    with (HERE.parent/"radial_followup_gaia.csv").open() as infile:
        positions = [f"{row['ra']},{row['dec']},8" for row in csv.DictReader(infile)]
    session = build_opener(HTTPCookieProcessor(CookieJar()))
    with session.open(SEARCH, timeout=60) as response:
        form = response.read().decode()
    token = re.search(r'name="csrf_token" value="([^"]+)"', form)
    if not token:
        raise ValueError("The public query form changed")
    fields = ["obsid","designation","gaia_source_id","ra","dec","z","rv","rv_err",
              "class","subclass","teff","logg","feh","snrg"]
    params = {
        "csrf_token": token[1],
        "output.fmt": "csv",
        "pos.type": "proximity",
        "pos.radecTextarea": "\n".join(positions),
        "pos.defaultRadius": "8",
        "showcol": ["combined." + name for name in fields],
    }
    request = Request(QUERY, data=urlencode(params,doseq=True).encode(),headers={"Referer":SEARCH})
    with session.open(request, timeout=90) as response:
        data = response.read()
    header = next(csv.reader(io.StringIO(data.decode("utf-8-sig"))), [])
    if not {"obsid","ra","dec","rv","rv_err"}.issubset(header):
        raise ValueError("The returned file has an unexpected schema")
    (HERE/"lamost_dr11.csv").write_bytes(data)


if __name__ == "__main__":
    main()
