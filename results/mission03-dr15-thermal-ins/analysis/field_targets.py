"""Read the observation target from the PPS source-list header for selected crowded fields."""

import csv
import hashlib
import gzip
import io
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from astropy.io import fits

OBSID = (
    "0111090301",
    "0503220101",
    "0112870401",
    "0113050601",
    "0077340201",
    "0077340101",
)


def fetch(obs):
    url = f"https://heasarc.gsfc.nasa.gov/FTP/xmm/data/rev0/{obs}/PPS/P{obs}EPX000OBSMLI0000.FTZ"
    with urlopen(url, timeout=45) as r:
        content = r.read()
    with fits.open(
        io.BytesIO(
            gzip.decompress(content)
            if content.startswith(bytes([31, 139]))
            else content
        )
    ) as f:
        target = f[0].header["OBJECT"]
    return dict(
        obsid=obs, target=target, url=url, sha256=hashlib.sha256(content).hexdigest()
    )


if __name__ == "__main__":
    with ThreadPoolExecutor(3) as ex:
        print(json.dumps(list(ex.map(fetch, OBSID)), indent=2))
