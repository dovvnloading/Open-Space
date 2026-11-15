"""Fetch the native images and catalog rows used by the comparison, with byte checks."""

import csv
import hashlib
from urllib.parse import urlencode
from urllib.request import urlopen

from screen import CACHE, HERE, URLS

ESO = "https://dataportal.eso.org/dataPortal/soda/sync?"
IDS = [
    "2016-08-09T11:09:52.408",
    "2016-08-09T11:09:52.411",
    "2016-08-09T11:09:52.413",
    "2016-08-09T11:09:53.840",
    "2016-08-09T11:09:53.845",
    "2019-10-07T14:31:55.021",
    "2019-10-07T14:31:55.024",
]
MASH = "https://vizier.cfa.harvard.edu/ftp/cats/V/127A/"
INPUTS = {
    **URLS,
    "radio.fits": "https://archive-gw-1.kat.ac.za/public/repository/10.48479/0n8c-5q84/data/SMGPS_fits/SMGPS_Candidate_G276.229%2B00.419.fits",
    "ha.fits": MASH + "fits/1102_ha.fits",
    "sr.fits": MASH + "fits/1102_sr.fits",
    "mash-object.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=V/127A&-out.all&-c=145.505%20-52.345&-c.rm=2",
    "spectrum-index.dat": MASH + "sp.dat",
}
for obs in IDS:
    INPUTS["ADP." + obs + ".fits"] = ESO + urlencode(
        dict(ID="ivo://eso.org/ID?ADP." + obs, CIRCLE="145.505 -52.345 0.06")
    )


def digest(data, tsv=False):
    # ASU includes request timestamp and server metadata in its comment header.
    if tsv:
        data = b"\n".join(l for l in data.splitlines() if l and not l.startswith(b"#"))
    return hashlib.sha256(data).hexdigest()


def main():
    CACHE.mkdir(exist_ok=True)
    with (HERE / "image-inputs.csv").open() as f:
        frozen = {r["file"]: r for r in csv.DictReader(f)}
    for name, url in INPUTS.items():
        path = CACHE / name
        if not path.exists():
            with urlopen(url, timeout=180) as resp:
                data = resp.read()
            # Verify before committing the cache file.
        else:
            data = path.read_bytes()
        if digest(data, name.endswith(".tsv")) != frozen[name]["sha256"]:
            raise ValueError(f"Input differs from frozen release: {name}")
        if not path.exists():
            path.write_bytes(data)


if __name__ == "__main__":
    main()
