"""Extract original VLASS pixels using two HTTP byte ranges and archive WCS."""

import io
from urllib.parse import quote
from urllib.request import Request, urlopen
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ARCHIVE = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VLASS/"


def byte_range(url, start, end):
    request = Request(url, headers={"Range": f"bytes={start}-{end}"})
    with urlopen(request, timeout=35) as response:
        if response.status != 206 or not response.headers.get(
            "Content-Range", ""
        ).startswith(f"bytes {start}-{end}/"):
            raise ValueError("Archive returned a different byte range")
        data = response.read()
    if len(data) != end - start + 1:
        raise ValueError("Incomplete archive byte range")
    return data


def extract(uri, ra, dec):
    url = ARCHIVE + quote(uri.rsplit("/", 1)[-1])
    raw = byte_range(url, 0, 143999)
    while True:
        found = next(
            (j + 80 for j in range(0, len(raw), 80) if raw[j : j + 8] == b"END     "),
            None,
        )
        if found is not None:
            end = found
            break
        if len(raw) >= 2304000:
            raise ValueError("FITS primary header exceeds 2.3 MB")
        raw += byte_range(url, len(raw), 2 * len(raw) - 1)
    header = fits.Header.fromstring(raw[:end].decode("ascii"))
    offset = 2880 * ((end + 2879) // 2880)
    if (
        header["BITPIX"] != -32
        or header.get("BSCALE", 1) != 1
        or header.get("BZERO", 0) != 0
    ):
        raise ValueError("Unsupported primary FITS encoding")
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    x, y = np.rint(WCS(header).celestial.world_to_pixel_values(ra, dec)).astype(int)
    left, right = max(0, x - 45), min(nx, x + 46)
    bottom, top = max(0, y - 45), min(ny, y + 46)
    if right <= left or top <= bottom:
        raise ValueError("Target outside stored image")
    data = byte_range(url, offset + bottom * nx * 4, offset + top * nx * 4 - 1)
    tile = np.frombuffer(data, dtype=">f4").reshape(top - bottom, nx)[:, left:right]
    header["CRPIX1"] -= left
    header["CRPIX2"] -= bottom
    header["OSRANGE"] = (True, "Cut from original NRAO pixels using HTTP Range")
    out = io.BytesIO()
    fits.PrimaryHDU(data=tile[None, None, :, :], header=header).writeto(
        out, output_verify="silentfix"
    )
    return out.getvalue()
