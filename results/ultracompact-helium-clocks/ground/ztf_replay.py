"""Localize positive difference light from the 2023 December MASTER OT J0558+3915 eruption."""
import csv
import gzip
import hashlib
import json
import math
from datetime import datetime, timedelta
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
FINGERPRINT = "df9f8011b4cddd811fe8901a5d301e749225245718f3c55e14ff3e143bd41ed3"
MJD0 = datetime(1858, 11, 17)


def offset(ra, dec, ra0, dec0):
    return math.hypot((ra-ra0)*math.cos(math.radians(dec0)), dec-dec0)*3600


def main():
    raw = gzip.decompress((ROOT / "ztf-alerts.json.gz").read_bytes())
    if hashlib.sha256(raw).hexdigest() != FINGERPRINT:
        raise ValueError("Incorrect ALeRCE snapshot")
    rows = [r for r in json.loads(raw) if 60288 <= r["mjd"] < 60309]
    # Some exposures have two alert candidate IDs but identical subtraction photometry.
    seen = set()
    selected = []
    alerts = []
    for r in sorted(rows, key=lambda row: (row["mjd"], row["candid"])):
        key = r["pid"]
        used = (r["mjd"] >= 60290 and r["isdiffpos"] == 1 and r["magpsf"] < 16 and
                key not in seen)
        if used:
            selected.append(r)
            seen.add(key)
        alerts.append(dict(mjd=r["mjd"], candid=r["candid"], pid=key, fid=r["fid"],
                           isdiffpos=r["isdiffpos"], difference_mag=r["magpsf"],
                           mag_uncertainty=r["sigmapsf"], ra=r["ra"], dec=r["dec"],
                           rb=r["rb"], used_for_centroid=used))
    ra = median(r["ra"] for r in selected)
    dec = median(r["dec"] for r in selected)
    mjd = median(r["mjd"] for r in selected)
    years = ((MJD0 + timedelta(days=mjd) - datetime(2016, 1, 1, 12)).total_seconds() /
             (365.25 * 86400))
    matches = []
    with (ROOT / "gaia-cone.csv").open() as f:
        for r in csv.DictReader(f):
            a, d = float(r["ra"]), float(r["dec"])
            # When proper motion is unavailable, use the published position.
            a += float(r["pmra"] or 0)*years/(3.6e6*math.cos(math.radians(d)))
            d += float(r["pmdec"] or 0)*years/3.6e6
            matches.append(dict(source_id=r["source_id"], epoch_ra=a, epoch_dec=d,
                                separation_arcsec=round(offset(a, d, ra, dec), 4)))
    matches.sort(key=lambda row: row["separation_arcsec"])
    scatter = [offset(r["ra"], r["dec"], ra, dec) for r in selected]
    filters = {}
    for band in (1, 2):
        subset = [r for r in selected if r["fid"] == band]
        filters[str(band)] = {"exposures": len(subset), "ra": median(r["ra"] for r in subset),
                              "dec": median(r["dec"] for r in subset)}
    summary = dict(n=len(selected), ra_deg=ra, dec_deg=dec, median_mjd=mjd,
                   median_scatter_arcsec=round(median(scatter), 4),
                   max_scatter_arcsec=round(max(scatter), 4),
                   first_mjd=selected[0]["mjd"], last_mjd=selected[-1]["mjd"],
                   distinct_mjd_nights=len({int(r["mjd"]) for r in selected}), by_filter=filters,
                   nearest_gaia=matches[:5])
    print(json.dumps(dict(summary=summary, alerts=alerts), indent=2))


if __name__ == "__main__":
    main()
