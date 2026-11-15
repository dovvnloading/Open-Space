"""Compare the published radial velocity track to the coordinates of the measured stars."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

from reproduce import get_sources

URL = ("https://raw.githubusercontent.com/cmateu/galstreams/"
       "6a304c05bea7f774390a11082ea127719e99b70a/galstreams/tracks/"
       "track.st.NGC1261b.ibata2024.ecsv")
SHA256 = "03a116fa782fc7c9a59ec31daff0c0bcc8758e5a65461755a89c27ffdf7f3a31"
C_KMS = 299792.458


def audit(cache, output):
    get_sources(cache)
    path = cache/"track.ecsv"
    if path.exists():
        content = path.read_bytes()
    else:
        with urlopen(URL, timeout=90) as response:
            content = response.read()
    if hashlib.sha256(content).hexdigest() != SHA256:
        raise ValueError("Unexpected Galstreams ECSV bytes")
    if not path.exists():
        path.write_bytes(content)
    lines = content.decode().splitlines()
    if not any("name: radial_velocity, unit: km / s," in row for row in lines):
        raise ValueError("Expected signed radial velocity in km/s")
    rows = list(csv.DictReader(row for row in lines if row and not row.startswith("#")))
    samples = [(float(row["ra"]), float(row["radial_velocity"])) for row in rows]
    measured = [row.split() for row in (cache/"ibata.txt").read_text().splitlines()
                if row[:19].strip().isdigit() and row[95:97].strip() == "14"
                and int(row[92:94]) > 0 and float(row[85:91]) < 300]
    ra_lo, ra_hi = min(float(row[1]) for row in measured), max(float(row[1]) for row in measured)
    if any(samples[i][0] >= samples[i+1][0] for i in range(len(samples)-1)):
        raise ValueError("RA must increase along this track to define on-track support")
    inside = [v for ra, v in samples if ra_lo <= ra <= ra_hi]
    speeds = [v for _, v in samples]
    result = dict(track_url=URL, track_sha256=SHA256, measured_heliocentric_rvs=len(measured),
                  measured_ra_min_deg=ra_lo, measured_ra_max_deg=ra_hi,
                  track_rows=len(samples), rows_in_measured_ra_span=len(inside),
                  min_velocity_kms=min(speeds), max_velocity_kms=max(speeds),
                  c_kms=C_KMS, rows_abs_velocity_over_c=sum(abs(v)>C_KMS for v in speeds),
                  rows_abs_velocity_over_c_in_measured_ra_span=sum(abs(v)>C_KMS for v in inside))
    output.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.cache, args.output)
