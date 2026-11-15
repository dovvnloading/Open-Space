"""List ZTF alert objects at each flare position using the ALeRCE broker.

Usage: python ztf_alerts.py SUMMARY.csv OUT.json [--radius ARCSEC]

Only positions north of declination -28 deg are queried. For every alert
object within the radius the number of detections and the first and last
detection dates are stored. Many detections spread over years indicate a
persistently variable nucleus; a short cluster indicates a transient.
"""
import argparse
import json

import pandas as pd
import requests

URL = "https://api.alerce.online/ztf/v1/objects"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("summary")
    parser.add_argument("out")
    parser.add_argument("--radius", type=float, default=5.0)
    args = parser.parse_args()
    table = pd.read_csv(args.summary, dtype={"id": str})
    out = {}
    for _, row in table[table.dec > -28].iterrows():
        response = requests.get(URL, params={"ra": row.ra, "dec": row.dec, "radius": args.radius, "page_size": 50},
                                timeout=120)
        response.raise_for_status()
        out[row.id] = [{"oid": item["oid"], "detections": item.get("ndet"), "first_mjd": item.get("firstmjd"),
                        "last_mjd": item.get("lastmjd")} for item in response.json().get("items", [])]
    json.dump(out, open(args.out, "w"), indent=1)
    print(len(out), "positions;", sum(bool(v) for v in out.values()), "with alert objects")


if __name__ == "__main__":
    main()
