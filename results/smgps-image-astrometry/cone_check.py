"""Query the live VizieR image index around native and published image centers."""

import csv
import time
from urllib.parse import urlencode

import astropy.units as u
from astropy.coordinates import SkyCoord
from harvest import HERE, get

BASE = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
EXAMPLES = [
    "SMGPS_Candidate_G043.942+01.442.fits",
    "SMGPS_Candidate_G051.746+01.384.fits",
    "SMGPS_Candidate_G057.868+00.221.fits",
]


def main():
    results = []
    for row in csv.DictReader((HERE / "coordinates.csv").open()):
        if row["filename"] not in EXAMPLES:
            continue
        for origin in ["native", "listed"]:
            c = SkyCoord(
                l=float(row[origin + "_l"]) * u.deg,
                b=float(row[origin + "_b"]) * u.deg,
                frame="galactic",
            ).icrs
            url = BASE + urlencode(
                {
                    "-source": "J/A+A/693/A247/list",
                    "-c": f"{c.ra.deg:.7f} {c.dec.deg:+.7f}",
                    "-c.rm": "5",
                    "-out": "FileName",
                    "-out.max": "unlimited",
                }
            )
            for attempt in range(5):
                response = get(url).decode()
                # An HTML outage page is not a null cone search.
                if (
                    "VizieR Astronomical Server" in response
                    and "#RESOURCE=" in response
                ):
                    break
                time.sleep(2**attempt)
            else:
                raise RuntimeError(f"Unverified VizieR response for {url}")
            filenames = [
                line.strip()
                for line in response.splitlines()
                if line.strip().endswith(".fits")
            ]
            results.append(
                [
                    row["filename"],
                    origin,
                    int(row["filename"] in filenames),
                    len(filenames),
                    url,
                ]
            )
    with (HERE / "cone-check.csv").open("w", newline="") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(
            ["filename", "query_center", "self_retrieved", "number_of_rows", "url"]
        )
        w.writerows(results)


if __name__ == "__main__":
    main()
