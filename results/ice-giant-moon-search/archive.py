"""Inventory public planet fields, enumerate calibrated exposures and hash downloaded bytes.

Examples (from this directory):
    python archive.py inventory manifests --end 2026-11-17
    python archive.py products manifests
    python archive.py fetch manifests/products.csv /data/moons --pattern 'jw08975005*_cal.fits'
End is exclusive UTC. Inventories use MAST metadata at the time of retrieval.
"""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import fnmatch
import hashlib
import json
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://mast.stsci.edu/api/v0/invoke"
DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file"
PLANETS = {"Uranus": "URAN", "Neptune": "NEPTUN"}
COLUMNS = ("obsid obs_id obs_collection instrument_name target_name t_min t_max "
           "t_exptime proposal_id calib_level filters dataRights dataproduct_type "
           "s_ra s_dec mtFlag intentType").split()
OBS_FIELDS = ["planet", "selected_by", "sep_arcmin_mid", *COLUMNS]
PRODUCT_FIELDS = ["planet", "parent_obsids", "obs_collection", "proposal_id", "obs_id",
                  "filters", "productFilename", "productSubGroupDescription", "size", "dataURI"]
FILE_FIELDS = ["dataURI", "productFilename", "bytes", "sha256"]


def csvrows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, fields, rows):
    temp = path.with_suffix(".tmp")
    with temp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def values(name, items):
    return {"paramName": name, "values": items}


def between(name, lower, upper):
    return values(name, [{"min": float(lower), "max": float(upper)}])


def invoke(service, params, page_size=5000):
    page = 1
    while True:
        payload = dict(service=service, params=params, format="json", pagesize=page_size, page=page)
        body = urlencode({"request": json.dumps(payload)}).encode()
        for attempt in range(6):
            try:
                with urlopen(Request(API, data=body), timeout=120) as response:
                    result = json.load(response)
                if result["status"] == "COMPLETE":
                    break
                raise RuntimeError(result.get("msg") or result["status"])
            except (OSError, ValueError, RuntimeError):
                if attempt == 5:
                    raise
                time.sleep(min(2 ** attempt, 16))
        yield from result["data"]
        if page >= result.get("paging", {}).get("pagesFiltered", 1):
            break
        page += 1


def image_mode(row):
    instr = row["instrument_name"].upper().split("/")[0]
    return (instr in {"WFC3", "MIRI", "NIRCAM", "NIRISS"} and
            row["dataproduct_type"] == "image")


def inventory(args):
    import numpy as np
    from astropy.coordinates import get_body
    from astropy.time import Time

    start, end = Time(args.start).mjd, Time(args.end).mjd
    if start >= end:
        raise ValueError("end must follow start")
    days = np.arange(np.floor(start), np.ceil(end) + 1.)
    base = [values("obs_collection", ["HST", "JWST"]), values("dataRights", ["PUBLIC"]),
            between("t_min", start, end)]
    jobs = []
    ephem = {}
    for planet, stem in PLANETS.items():
        # Geocentric daily apparent positions. The 10 arcmin radius more than covers
        # Earth/JWST parallax, field size and low-precision built-in ephemeris errors.
        c = get_body(planet.lower(), Time(days, format="mjd"))
        ra = np.rad2deg(np.unwrap(c.ra.radian))
        dec = c.dec.degree
        ephem[planet] = (ra, dec)
        jobs.append((planet, "name", [dict(paramName="target_name", values=[], freeText=f"%{stem}%")]))
        pad = 0.35
        low, high = ra.min() - pad, ra.max() + pad
        bounds = [(low % 360, high % 360)]
        if int(low // 360) != int(high // 360):
            bounds = [(0, high % 360), (low % 360, 360)]
        for lo, hi in bounds:
            jobs.append((planet, "sky", [between("s_ra", lo, hi), between("s_dec", dec.min()-pad, dec.max()+pad)]))
    selected = {}
    def query(job):
        planet, mode, extra = job
        rows = list(invoke("Mast.Caom.Filtered", dict(columns=",".join(COLUMNS), filters=base+extra)))
        return planet, mode, rows
    with ThreadPoolExecutor(3) as pool:
        for planet, mode, rows in pool.map(query, jobs):
            ra, dec = ephem[planet]
            for row in rows:
                if not image_mode(row) or row["t_min"] >= end:
                    continue
                mid = (row["t_min"] + (row["t_max"] or row["t_min"])) / 2
                sep = None
                if row["s_ra"] is not None and row["s_dec"] is not None:
                    phi = np.interp(mid, days, dec)
                    delta_ra = (row["s_ra"] - np.interp(mid, days, ra) + 180) % 360 - 180
                    sep = 60 * np.hypot(delta_ra * np.cos(np.deg2rad(phi)), row["s_dec"] - phi)
                if mode == "sky" and (sep is None or sep > 10):
                    continue
                key = (planet, row["obsid"])
                if key in selected:
                    selected[key]["selected_by"] = "name+sky"
                else:
                    selected[key] = dict(row, planet=planet, selected_by=mode,
                                         sep_arcmin_mid="" if sep is None else round(sep, 4))
    args.out.mkdir(parents=True, exist_ok=True)
    ordered = sorted(selected.values(), key=lambda r: (r["planet"], r["t_min"], r["obs_id"]))
    write_csv(args.out / "observations.csv", OBS_FIELDS, ordered)
    info = dict(start_utc=args.start, end_exclusive_utc=args.end, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                endpoint=API, service="Mast.Caom.Filtered", columns=COLUMNS, filters=base,
                query_union=jobs, sky_radius_arcmin=10, unique_observations=len(ordered))
    (args.out / "inventory_query.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"{len(ordered)} candidate observation groups")


def products(args):
    obs = csvrows(args.out / "observations.csv")
    parents = {r["obsid"]: r for r in obs}
    ids = sorted(parents)
    batches = [ids[i:i+20] for i in range(0, len(ids), 20)]
    def query(batch):
        rows = invoke("Mast.Caom.Products", {"obsid": ",".join(batch)}, page_size=5000)
        for r in rows:
            if (r.get("dataRights") == "PUBLIC" and r.get("productType") == "SCIENCE" and
                    r["productSubGroupDescription"] in {"CAL","CALINTS","RATEINTS","FLT","FLC"} and
                    r["productFilename"].lower().endswith(".fits")):
                yield r
    def collect(batch):
        return list(query(batch))
    selected = {}
    with ThreadPoolExecutor(4) as pool:
        for rows in pool.map(collect, batches):
            for row in rows:
                parent = str(row["parent_obsid"])
                if parent not in parents:
                    raise ValueError(f"Unknown product parent {parent}")
                key = (parents[parent]["planet"], row["dataURI"])
                if key not in selected:
                    selected[key] = dict(row, planet=key[0], parent_obsids=set())
                selected[key]["parent_obsids"].add(parent)
    ordered = sorted(selected.values(), key=lambda r: (r["planet"], r["productFilename"]))
    for row in ordered:
        row["parent_obsids"] = "|".join(sorted(row["parent_obsids"]))
    write_csv(args.out / "products.csv", PRODUCT_FIELDS, ordered)
    print(f"{len(ordered)} distinct planet/product rows")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fetch(args):
    rows = [r for r in csvrows(args.products) if fnmatch.fnmatchcase(r["productFilename"], args.pattern)]
    if not rows:
        raise ValueError(f"No matching product filenames for {args.pattern}")
    args.dest.mkdir(parents=True, exist_ok=True)
    ledger = args.dest / "downloaded.csv"
    previous = {r["productFilename"]: r for r in csvrows(ledger)} if ledger.exists() else {}
    unique = {r["productFilename"]:r for r in rows}
    for filename, row in sorted(unique.items()):
        if Path(filename).name != filename:
            raise ValueError(filename)
        path = args.dest / filename
        size = int(row["size"])
        if path.exists() and path.stat().st_size == size:
            sha = digest(path)
            if filename in previous and previous[filename]["sha256"] != sha:
                raise ValueError(f"Downloaded file changed: {path}")
        else:
            tmp = path.with_suffix(".part")
            url = DOWNLOAD + "?" + urlencode({"uri": row["dataURI"]})
            with urlopen(url, timeout=180) as src, tmp.open("wb") as dest:
                while data := src.read(2**20):
                    dest.write(data)
            if tmp.stat().st_size != size:
                raise OSError(f"Incomplete transfer: {path.name}")
            with tmp.open("rb") as handle:
                if not handle.read(80).startswith(b"SIMPLE  =                    T"):
                    raise ValueError(f"Invalid FITS signature: {path.name}")
            tmp.replace(path)
            sha = digest(path)
        previous[filename] = dict(dataURI=row["dataURI"], productFilename=filename, bytes=size, sha256=sha)
        write_csv(ledger, FILE_FIELDS, sorted(previous.values(), key=lambda r:r["productFilename"]))
    print(f"{len(unique)} files hashed")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="mode", required=True)
    for name, func in (("inventory", inventory), ("products", products)):
        sub = subs.add_parser(name)
        sub.add_argument("out", type=Path)
        sub.set_defaults(run=func)
        if name == "inventory":
            sub.add_argument("--start", default="2022-01-01")
            sub.add_argument("--end", default="2026-11-17")
    sub = subs.add_parser("fetch")
    sub.add_argument("products", type=Path)
    sub.add_argument("dest", type=Path)
    sub.add_argument("--pattern", required=True)
    sub.set_defaults(run=fetch)
    args = p.parse_args()
    args.run(args)
