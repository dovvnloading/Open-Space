"""Identify a PPS source by sky position, then fetch its OGIP spectral products."""

import argparse
import hashlib
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError
from astropy.io import fits
from astropy.coordinates import SkyCoord

BASE = "https://heasarc.gsfc.nasa.gov/FTP/xmm/data"
RESPONSE = BASE + "/responses/pipeline/"


def get(url, path):
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(url, timeout=90) as response:
                    payload = response.read()
                break
            except HTTPError as exc:
                if exc.code == 404 or attempt == 2:
                    raise
                time.sleep(2**attempt)
        temp = path.with_suffix(path.suffix + ".part")
        temp.write_bytes(payload)
        temp.replace(path)
    return path


def retrieve(obsid, ra, dec, dest, max_sep):
    base = f"{BASE}/rev0/{obsid}/PPS/"
    manifest = {"obsid": obsid, "input_ra": ra, "input_dec": dec, "products": []}

    def product(name, prefix=base):
        path = get(prefix + name, dest / name)
        manifest["products"].append(
            {
                "file": name,
                "url": prefix + name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
        return path

    index = get(base, dest / (obsid + ".html")).read_text()
    files = re.findall(r'href="([^"/]+)"', index)
    srclist = product(f"P{obsid}EPX000OBSMLI0000.FTZ")
    with fits.open(srclist) as f:
        cat = f[1].data
        angles = (
            SkyCoord(cat["RA"], cat["DEC"], unit="deg")
            .separation(SkyCoord(ra, dec, unit="deg"))
            .arcsec
        )
        best = int(angles.argmin())
        manifest.update(
            pps_source=int(cat["SRC_NUM"][best]),
            offset_arcsec=float(angles[best]),
            pps_ra=float(cat["RA"][best]),
            pps_dec=float(cat["DEC"][best]),
            target=f[0].header.get("OBJECT"),
        )
    if manifest["offset_arcsec"] > max_sep:
        manifest["status"] = "no PPS source inside match radius"
        return manifest
    tail = f"SRSPEC0{manifest['pps_source']:03X}.FTZ"
    spectra = [n for n in files if n[11:13] in ("M1", "M2", "PN") and n.endswith(tail)]
    for spectrum in spectra:
        with fits.open(product(spectrum)) as f:
            h = f["SPECTRUM"].header
            arf = h["ANCRFILE"].strip()
            rmf = h["RESPFILE"].strip()
            bg = h["BACKFILE"].strip()
        product(bg)
        product(arf)
        try:
            product(rmf, RESPONSE)
        except HTTPError as exc:
            if exc.code != 404:
                raise
            product(rmf + ".gz", RESPONSE)
    manifest["status"] = "spectra fetched" if spectra else "matched, no PPS spectrum"
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ra", type=float)
    p.add_argument("dec", type=float)
    p.add_argument("dest", type=Path)
    p.add_argument("obsids", nargs="+")
    p.add_argument(
        "--max-sep",
        type=float,
        default=4,
        help="Maximum PPS identification distance (arcsec)",
    )
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    for obs in a.obsids:
        if not re.fullmatch(r"\d{10}", obs):
            raise ValueError("Expected ten-digit XMM observation ID")
        m = retrieve(obs, a.ra, a.dec, a.dest, a.max_sep)
        (a.dest / (obs + ".json")).write_text(json.dumps(m, indent=2) + "\n")
        print(
            obs,
            m["target"],
            m["status"],
            m["pps_source"],
            round(m["offset_arcsec"], 2),
            flush=True,
        )


if __name__ == "__main__":
    main()
