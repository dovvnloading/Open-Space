"""Recreate the image sample and verify the per-image SHA-256 digests."""
import argparse
import csv
import gzip
import hashlib
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

def rows(name):
    with (DATA / name).open(newline="") as stream:
        return list(csv.DictReader(stream))

def choose_controls(targets, peers):
    """Prefer up to three same-CCD peers, then fill to eight by S/N and PSF size."""
    pairs = []
    for target in targets:
        selection = [p for p in peers if p["visit"] == target["visit"] and p["oid"] != target["oid"]]
        def distance(peer):
            tsize = float(target["ixxpsf"]) + float(target["iyypsf"])
            psize = float(peer["ixxpsf"]) + float(peer["iyypsf"])
            return abs(math.log(float(peer["snr"]) / float(target["snr"]))) + 2 * abs(math.log(psize / tsize))
        selection = sorted(selection, key=lambda p: (distance(p), int(p["measurement_id"])))
        same = [p for p in selection if p["detector"] == target["detector"]][:3]
        matched = same + [p for p in selection if p not in same][:8-len(same)]
        for peer in matched:
            pairs.append(dict(target_id=target["measurement_id"], control_id=peer["measurement_id"],
                              visit=peer["visit"], detector=peer["detector"], snr=peer["snr"],
                              same_detector=(peer["detector"] == target["detector"])))
    return pairs

def url(mid, kind):
    return "https://api.lsst.fink-portal.org/api/v1/cutouts?" + urlencode(
        {"diaSourceId": mid, "kind": kind, "output-format": "FITS"})

def read_hdus(content):
    """Read these uncompressed image FITS HDUs, verifying supported data layout."""
    hdus = []
    offset = 0
    while offset < len(content):
        header, pos = {}, offset
        while pos + 80 <= len(content):
            card = content[pos:pos+80].decode("ascii")
            pos += 80
            if card[:8].strip() == "END":
                break
            if card[8:10] == "= ":
                header[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
        if header.get("NAXIS") != "2":
            raise ValueError("Expected a two-dimensional FITS image")
        nx, ny = int(header["NAXIS1"]), int(header["NAXIS2"])
        dtype = {-32: ">f4", -64: ">f8", 32: ">i4", 16: ">i2"}[int(header["BITPIX"])]
        start = offset + math.ceil((pos-offset)/2880)*2880
        size = nx*ny*np.dtype(dtype).itemsize
        plane = np.frombuffer(content[start:start+size], dtype=dtype).reshape(ny,nx).astype(float)
        hdus.append((header, plane))
        offset = start + math.ceil(size/2880)*2880
    if [h.get("EXTNAME") for h, _ in hdus] != [None, "UNCERT", "PSFIMAGE"]:
        raise ValueError("Unexpected Rubin cutout HDUs")
    return hdus

def read_cutout(mid, kind):
    return read_hdus(gzip.decompress((DATA / "stamps" / f"{mid}_{kind}.fits.gz").read_bytes()))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Download the same packet images from Fink")
    args = parser.parse_args()
    for r in rows("stamps.csv"):
        file = DATA / "stamps" / r["file"]
        if args.refresh:
            with urlopen(r["url"], timeout=60) as response:
                body = response.read()
        else:
            body = gzip.decompress(file.read_bytes())
        if hashlib.sha256(body).hexdigest() != r["sha256_fits"]:
            raise ValueError(f"{file} does not match its recorded image")
        read_hdus(body)
        if args.refresh:
            file.parent.mkdir(exist_ok=True)
            file.write_bytes(gzip.compress(body, mtime=0))
    print("Verified all FITS cutout digests and HDU layouts.")
