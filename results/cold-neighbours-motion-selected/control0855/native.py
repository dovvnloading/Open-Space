"""Read the header and required image rows of an uncompressed public NEOWISE L1b FITS."""
import hashlib
import math
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def fetch_range(url, start, end):
    request = Request(f"{url}?x-id=GetObject&range={start}-{end}",
                      headers={"Range": f"bytes={start}-{end}"})
    for attempt in range(4):
        try:
            with urlopen(request, timeout=75) as response:
                must_start = f"bytes {start}-{end}/"
                if not response.headers.get("Content-Range", "").startswith(must_start):
                    raise ValueError("Unexpected FITS byte range")
                result = response.read()
            if len(result) != end-start+1:
                raise ValueError("Truncated FITS byte range")
            return result
        except (HTTPError, URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2**attempt)


def cutout(url, ra, dec, half_width=43):
    initial = fetch_range(url, 0, 65535)
    end = next(i+80 for i in range(0, len(initial)-80, 80)
               if initial[i:i+80] == b"END".ljust(80))
    offset = 2880 * math.ceil(end/2880)
    raw_header = initial[:offset]
    header = fits.Header.fromstring(raw_header.decode("ascii"), sep="")
    if header["BITPIX"] != -32 or header.get("BSCALE",1) != 1 or header.get("BZERO",0) != 0:
        raise ValueError("Unsupported calibrated image format")
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    cx, cy = WCS(header).celestial.world_to_pixel_values(ra, dec)
    x0, x1 = max(0,int(cx)-half_width), min(nx,int(cx)+half_width+1)
    y0, y1 = max(0,int(cy)-half_width), min(ny,int(cy)+half_width+1)
    if x1 <= x0 or y1 <= y0:
        raise ValueError("Image misses the field")
    strip = fetch_range(url, offset+4*nx*y0, offset+4*nx*y1-1)
    pixels = np.frombuffer(strip, dtype=">f4").reshape(y1-y0,nx)[:,x0:x1].copy()
    header["CRPIX1"] -= x0
    header["CRPIX2"] -= y0
    header["INHDRSHA"] = hashlib.sha256(raw_header).hexdigest()
    header["INRAWSHA"] = hashlib.sha256(strip).hexdigest()
    return fits.PrimaryHDU(pixels, header)
