"""Verify archived inputs; retrieve and verify the TESS pixel arrays on demand."""

import argparse
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent


def array_digest(path):
    # TESSCut regenerates DATE and CHECKSUM header cards. Hash the actual analysis inputs.
    digest = hashlib.sha256()
    with fits.open(path) as hdus:
        for name in ("TIME", "TIMECORR", "CADENCENO", "QUALITY", "FLUX"):
            array = hdus[1].data[name]
            digest.update(
                np.ascontiguousarray(
                    array.astype(array.dtype.newbyteorder("<"))
                ).tobytes()
            )
        for key in (
            "CRPIX1",
            "CRPIX2",
            "CRVAL1",
            "CRVAL2",
            "PC1_1",
            "PC1_2",
            "PC2_1",
            "PC2_2",
            "CDELT1",
            "CDELT2",
            "CTYPE1",
            "CTYPE2",
        ):
            digest.update(str(hdus[2].header[key]).encode() + b"\n")
    return digest.hexdigest()


def verify(item):
    name, source = item
    path = ROOT / name
    if name.startswith("raw/"):
        if not path.exists():
            path.parent.mkdir(exist_ok=True)
            with (
                urlopen(source["url"], timeout=180) as r,
                ZipFile(io.BytesIO(r.read())) as archive,
            ):
                files = [f for f in archive.namelist() if f.endswith(".fits")]
                if len(files) != 1:
                    raise ValueError(files)
                tmp = path.with_suffix(".tmp")
                tmp.write_bytes(archive.read(files[0]))
                if array_digest(tmp) != source["analysis_sha256"]:
                    raise ValueError(f"Archive pixel data or WCS changed for {name}")
                tmp.replace(path)
        actual = array_digest(path)
        expected = source["analysis_sha256"]
    else:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        expected = source["sha256"]
    if actual != expected:
        raise ValueError(f"Digest mismatch for {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pixels", action="store_true", help="retrieve the six original TESS cutouts"
    )
    args = parser.parse_args()
    sources = json.loads((ROOT / "sources.json").read_text())
    items = [i for i in sources.items() if args.pixels or not i[0].startswith("raw/")]
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(verify, items))
    print(f"Verified {len(items)} inputs")
