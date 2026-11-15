"""Fold three spatially resolved ZTF light curves on the Gaia binary ephemeris."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

BASE = Path(__file__).resolve().parent
STARS = {
    "source": (296.18363924, 43.97936911),
    "target": (296.18488951, 43.98079315),
    "south": (296.18520307, 43.97862836),
}
FIELDS = ["oid", "mjd", "filtercode", "mag", "magerr", "catflags"]
GAIA = dict(reference_time=2234.0540609516934, frequency=1.9302825156283367480,
            primary_phase=0.7466577, secondary_phase=0.2377278)
PALOMAR = EarthLocation.from_geodetic(-116.863*u.deg,33.357*u.deg,1706*u.m)
iers.conf.auto_download = False


def retrieve():
    """Save the measured columns. Keep the row flags so filtering can be repeated."""
    receipts = {}
    for name, (ra, dec) in STARS.items():
        params = {"POS": f"CIRCLE {ra} {dec} 0.0005", "BANDNAME": "g,r,i",
                  "FORMAT": "csv", "BAD_CATFLAGS_MASK": "32768"}
        url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(params)
        raw = urlopen(url, timeout=120).read()
        rows = list(csv.DictReader(io.StringIO(raw.decode())))
        receipts[name] = {"url": url, "rows": len(rows), "raw_sha256": hashlib.sha256(raw).hexdigest()}
        with gzip.open(BASE/f"ztf_{name}.csv.gz", "wt") as f:
            writer = csv.DictWriter(f, FIELDS, extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (BASE/"ztf_receipts.json").write_text(json.dumps(receipts, indent=2)+"\n")


def read_star(name):
    with gzip.open(BASE/f"ztf_{name}.csv.gz", "rt") as f:
        rows = [r for r in csv.DictReader(f) if int(r["catflags"]) == 0]
    mjd = np.array([float(r["mjd"]) for r in rows])
    coord = SkyCoord(*STARS[name], unit=u.deg)
    time = Time(mjd, format="mjd", scale="utc", location=PALOMAR)
    bjd = (time.tdb + time.light_travel_time(coord)).jd
    return bjd, np.array([r['filtercode'] for r in rows]), np.array([r['oid'] for r in rows]), np.array([float(r['mag']) for r in rows])


def wrap(phi):
    return (phi + 0.5) % 1 - 0.5


def contrast(mag, inside, outside):
    return np.median(mag[inside])-np.median(mag[outside])


def main():
    reference = Time(2455197.5+GAIA["reference_time"],format="jd",scale="tcb")
    # Scale frequency as well as epoch when converting a TCB day to TDB.
    period = (Time(reference.jd+1/GAIA["frequency"],format="jd",scale="tcb").tdb-reference.tdb).to_value(u.day)
    output=[]
    for name in STARS:
        bjd, bands, oids, mag = read_star(name)
        # A distinct reference zero point for every field/filter object ID.
        phase = ((bjd-reference.tdb.jd)/period) % 1
        first = wrap(phase-GAIA["primary_phase"]); second = wrap(phase-GAIA["secondary_phase"])
        outside = (np.abs(first)>.165) & (np.abs(second)>.15)
        for oid in np.unique(oids):
            group = oids==oid
            if np.any(group & outside):mag[group] -= np.median(mag[group & outside])
        for band in np.unique(bands):
            use=bands==band
            for label, dist in [('primary',first),('secondary',second)]:
                inside=use & (np.abs(dist)<.04)
                out=use & outside
                days=np.floor(bjd[inside | out]-.5)
                unique=np.unique(days)
                ii=np.flatnonzero(inside | out)
                by_night=[ii[days==d] for d in unique]
                rng=np.random.default_rng(8110754)
                sample=[]
                for _ in range(1000):
                    ids=np.concatenate([by_night[j] for j in rng.integers(len(unique),size=len(unique))])
                    sample.append(contrast(mag[ids],inside[ids],out[ids]))
                lo,hi=np.percentile(sample,[2.5,97.5])
                output.append(dict(star=name,band=band,eclipse=label,n=int(use.sum()),
                    eclipse_points=int(inside.sum()),reference_points=int(out.sum()),nights=len(unique),
                    depth_mag=round(contrast(mag,inside,out),6),ci95_low=round(lo,6),ci95_high=round(hi,6)))
    with (BASE/"ztf_contrasts.csv").open('w') as f:
        writer=csv.DictWriter(f,output[0].keys(),lineterminator="\n");writer.writeheader();writer.writerows(output)
    print('Wrote contrasts for',len(output),'star, band and eclipse combinations.')


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download",action="store_true",help="refresh ZTF snapshots from IRSA")
    if parser.parse_args().download:retrieve()
    main()
