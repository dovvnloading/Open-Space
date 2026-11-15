"""Download the two public archive inventories and the pinned UltracoolSheet."""

import argparse
import hashlib
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

BASE = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync"
UC = "https://zenodo.org/api/records/15802304/files/UltracoolSheet%20-%20Main.csv/content"
TILES = """SELECT obs_id,s_ra,s_dec,s_region FROM ivoa.ObsCore
WHERE obs_collection='VLASS' AND obs_publisher_did LIKE '%quicklook' AND obs_id LIKE 'VLASS4%'"""
FILES = """SELECT o.observationID,a.uri FROM caom2.Observation o
JOIN caom2.Plane p ON o.obsID=p.obsID JOIN caom2.Artifact a ON p.planeID=a.planeID
WHERE o.collection='VLASS' AND o.observationID LIKE 'VLASS4%'
AND p.productID LIKE '%.quicklook' AND a.uri LIKE '%.fits'"""


def tap_url(query):
    return (
        BASE
        + "?"
        + urlencode(
            dict(
                QUERY=" ".join(query.split()),
                LANG="ADQL",
                REQUEST="doQuery",
                FORMAT="csv",
                MAXREC=150000,
            )
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    for name, url in [
        ("ucd.csv", UC),
        ("tiles.csv", tap_url(TILES)),
        ("files.csv", tap_url(FILES)),
    ]:
        with urlopen(url, timeout=180) as response:
            data = response.read()
        sha = hashlib.sha256(data).hexdigest()
        if (
            name == "ucd.csv"
            and sha
            != "2307303072ebf729be6e7af6ec67c0a9feed6c2ab6f0adb35444bfab5724a184"
        ):
            raise ValueError("Unexpected UltracoolSheet release digest")
        (args.cache / name).write_bytes(data)
        print(name, sha)


if __name__ == "__main__":
    main()
