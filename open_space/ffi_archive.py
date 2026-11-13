"""Access coordinate-defined SPOC full-frame cutouts from the public MAST TESScut service."""

from hashlib import sha256
from io import BytesIO
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from zipfile import ZipFile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .tess_archive import atomic_bytes, download


SERVICE = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?"
MAX_BYTES = 160_000_000


def request_url(field, sector):
    ra, dec, size = float(field["ra_deg"]), float(field["dec_deg"]), int(field["size_pixels"])
    if not (0 <= ra < 360 and -90 < dec < 90 and math.isfinite(ra) and math.isfinite(dec)):
        raise ValueError("Expected finite ICRS coordinates")
    if not 11 <= size <= 41 or not 1 <= sector <= 9999:
        raise ValueError("Use 11–41 pixels and a positive sector")
    return SERVICE + urlencode({"ra": ra, "dec": dec, "y": size, "x": size, "units": "px", "sector": sector})


def acquire(field, product, cache):
    """Hash scientific arrays, source image names and WCS independently of generated header timestamps."""
    sector = product["sector"]
    url = request_url(field, sector)
    path = Path(cache) / (sha256(url.encode()).hexdigest() + ".fits")
    if path.exists():
        raw = path.read_bytes()
    else:
        with ZipFile(BytesIO(download(url, MAX_BYTES))) as archive:
            images = [entry for entry in archive.infolist() if entry.filename.endswith("_astrocut.fits")]
            if len(images) != 1 or images[0].file_size > MAX_BYTES:
                raise ValueError("Expected exactly one bounded SPOC detector cutout; move fields off detector edges")
            raw = archive.read(images[0])
    with fits.open(BytesIO(raw), memmap=False) as hdus:
        head, rows = hdus[0].header, hdus["PIXELS"].data
        if (head["FFI_TYPE"] != "SPOC" or head["TELESCOP"] != "TESS" or head["SIMDATA"]
                or head["SECTOR"] != sector or abs(head["RA_OBJ"]-field["ra_deg"]) > 1e-6
                or abs(head["DEC_OBJ"]-field["dec_deg"]) > 1e-6
                or head["TIMESYS"] != "TDB" or head["BJDREFI"]+head["BJDREFF"] != 2457000):
            raise ValueError("TESScut measurement identity differs from the plan")
        wcs = WCS(hdus["APERTURE"].header).celestial
        wcs_header = wcs.to_header(relax=True)
        science = sha256()
        science.update(json.dumps(dict(sorted(wcs_header.items())), sort_keys=True).encode())
        for key in ("SECTOR", "CAMERA", "CCD", "DATA_REL", "RA_OBJ", "DEC_OBJ", "TIMESYS", "TIMEUNIT",
                    "BJDREFI", "BJDREFF", "FFI_TYPE"):
            science.update(str(head[key]).encode() + b"\0")
        arrays = []
        for name, dtype in (("TIME", ">f8"), ("FLUX", ">f4"), ("FLUX_ERR", ">f4"), ("QUALITY", ">i4")):
            a = np.asarray(rows[name], dtype=dtype)
            science.update(a.tobytes(order="C"))
            arrays.append(a)
        source_images = [str(v).strip() for v in rows["FFI_FILE"]]
        science.update("\n".join(source_images).encode())
        if not source_images or not all(re.fullmatch(r"tess\d{13}-s\d{4}-\d-\d-\d{4}-s_ffi(?:c\.fits)?", v)
                                        for v in source_images):
            raise ValueError("Missing MAST parent FFI identifiers")
        n = field["size_pixels"]
        if arrays[1].shape[1:] != (n, n) or arrays[2].shape != arrays[1].shape:
            raise ValueError("Unexpected dimensions in TESScut")
        digest = science.hexdigest()
        if product.get("science_sha256") and product["science_sha256"] != digest:
            raise ValueError("TESScut scientific content differs from the pinned digest")
        valid = (arrays[3] == 0) & np.isfinite(arrays[0])
        valid &= (np.isfinite(arrays[1]) & np.isfinite(arrays[2]) & (arrays[2] > 0)).mean(axis=(1,2)) >= 0.98
        if valid.sum() < 200:
            raise ValueError("Insufficient valid FFI cadences")
        t, images, err = [np.asarray(a[valid], dtype=float) for a in arrays[:3]]
        if np.any(np.diff(t) <= 0) or np.ptp(t) < 8:
            raise ValueError("Insufficient unique FFI time baseline")
        meta = {"url": url, "sector": sector, "camera": head["CAMERA"], "ccd": head["CCD"],
                "data_release": head["DATA_REL"], "science_sha256": digest,
                "retrieved_fits_sha256": sha256(raw).hexdigest(), "retrieved_fits_bytes": len(raw),
                "input_images": len(source_images), "quality_and_finite_pass": len(t),
                "first_parent_ffi": source_images[0], "last_parent_ffi": source_images[-1],
                "parent_ffi_list_sha256": sha256("\n".join(source_images).encode()).hexdigest(),
                "start_btjd": t[0], "stop_btjd": t[-1], "time_reference": "BJD_TDB - 2457000"}
    if not path.exists():
        atomic_bytes(path, raw)
    return (t, images, err, wcs), meta
