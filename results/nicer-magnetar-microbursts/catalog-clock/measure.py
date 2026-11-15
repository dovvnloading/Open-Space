"""Compare the time scale on each published start to NICER event and orbit data."""

import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time, TimeDelta

ROOT = Path(__file__).resolve().parent
BASE = "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/"
SHIFTS = (0.0, 69.184, -69.184, 138.368)
COUNT_FIELDS = ("literal", "plus69", "minus69", "plus138")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def acquire(path, url, expected, download):
    if not path.exists() and download:
        tmp = path.with_suffix(path.suffix + ".partial")
        with urlopen(url, timeout=180) as resp, tmp.open("wb") as out:
            shutil.copyfileobj(resp, out)
        if digest(tmp) != expected:
            raise ValueError(f"Unexpected downloaded bytes at {url}")
        tmp.replace(path)
    if digest(path) != expected:
        raise ValueError(f"Unexpected input at {path}")


def epoch(header):
    if header["TIMESYS"].strip() != "TT":
        raise ValueError("Expected TT mission epoch")
    return Time(header["MJDREFI"], header["MJDREFF"], format="mjd", scale="tt")


def predicted_time(labels, ra, dec, eph, orbit):
    # Iteratively invert SSB arrival time = spacecraft time + geometric delay.
    # Geocentric spacecraft XYZ are in metres; RA/Dec are ICRS degrees.
    u, v = np.radians([ra, dec])
    direction = np.array([np.cos(v)*np.cos(u), np.cos(v)*np.sin(u), np.sin(v)])
    target = Time(labels, format="isot", scale="tdb")
    guess = (target.tt - eph).sec
    ot, xyz = orbit
    for _ in range(3):
        when = eph + TimeDelta(guess, format="sec")
        earth = get_body_barycentric("earth", when).xyz.to_value("m")
        spacecraft = np.stack([np.interp(guess, ot, axis) for axis in xyz])
        delay = direction @ (earth + spacecraft) / 299792458.0
        guess += (target - (when.tdb + TimeDelta(delay, format="sec"))).sec
    return guess


def measure(row, events, cached):
    obs = row["obsid"]
    with fits.open(cached / f"ni{obs}.orb.gz") as orbit, fits.open(cached / f"ni{obs}_0mpu7_cl.evt.gz") as fits_events:
        oh = orbit[1].header
        traj = orbit[1].data
        ev = fits_events["EVENTS"]
        eh = ev.header
        ref = epoch(eh)
        ot = traj["TIME"] + oh["TIMEZERO"] + (epoch(oh) - ref).sec
        xyz = [traj[k] for k in ("X", "Y", "Z")]
        ra, dec = eh["RA_OBJ"], eh["DEC_OBJ"]
        if events[0]["source"] == "CXOU J010043.1-721134":
            ra, dec = 15.17975, -72.19272222
        predicted = predicted_time([r["start"] for r in events], ra, dec, ref, (ot, xyz))
        pi = ev.data["PI"]
        t = np.sort(ev.data["TIME"][(pi >= 50) & (pi <= 800)] + eh["TIMEZERO"])
        gti = fits_events["GTI"]
        gh = gti.header
        offset = gh.get("TIMEZERO", 0.) + (epoch(gh) - ref).sec
        gs, ge = gti.data["START"] + offset, gti.data["STOP"] + offset
        out = []
        for src, center in zip(events, predicted):
            result = dict(src, ra=round(ra, 7), dec=round(dec, 7), predicted_met_tt_s=round(center, 6))
            # A millisecond offset matters for an isolated microburst; this
            # population comparison deliberately uses wider windows. Orbit
            # interpolation outside its endpoints is used only to test coverage.
            for margin in (0.25, 1.0):
                tag = "narrow" if margin == 0.25 else "wide"
                for key, shift in zip(COUNT_FIELDS, SHIFTS):
                    lo = center + shift - margin
                    hi = center + shift + float(src["t90_ms"]) / 1000. + margin
                    covered = (ot[0] <= lo and hi <= ot[-1] and np.any((gs <= lo) & (hi <= ge)))
                    result[f"{tag}_{key}"] = int(np.searchsorted(t, hi) - np.searchsorted(t, lo)) if covered else ""
            out.append(result)
        return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--download", action="store_true")
    args = p.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    master = json.loads((ROOT.parent / "inputs.json").read_text())
    catalog = args.cache / "pubtable.txt"
    acquire(catalog, master[catalog.name]["url"], master[catalog.name]["sha256"], args.download)
    group = defaultdict(list)
    for line in catalog.read_text().splitlines():
        obs = line[26:36]
        if re.fullmatch(r"[0-9]{10}", obs):
            group[obs].append(dict(obsid=obs, source=line[:21].strip(), number=line[22:25].strip(),
                                   start=line[37:60], t90_ms=line[61:69].strip()))
    with (ROOT / "inputs.csv").open() as f:
        files = list(csv.DictReader(f))
    if set(group) != {r["obsid"] for r in files}:
        raise ValueError("Archive inputs do not cover all catalog observations")
    with args.output.open("w", newline="") as f:
        writer = None
        for row in files:
            obs, month = row["obsid"], row["month"]
            for folder, name, key in [("xti/event_cl", f"ni{obs}_0mpu7_cl.evt.gz", "events_sha256"),
                                      ("auxil", f"ni{obs}.orb.gz", "orbit_sha256")]:
                acquire(args.cache/name, f"{BASE}{month}/{obs}/{folder}/{name}", row[key], args.download)
            measured = measure(row, group[obs], args.cache)
            if writer is None:
                writer = csv.DictWriter(f, measured[0], lineterminator="\n")
                writer.writeheader()
            writer.writerows(measured)


if __name__ == "__main__":
    main()
