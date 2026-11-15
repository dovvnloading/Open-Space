"""Read small pieces of the public SPHEREx DR1 HiPS spectral cubes over HTTP.

The cubes are order-5 HEALPix tiles (nested scheme, equatorial frame) of
102 spectral channels x 512 x 512 pixels, with IMAGE, NHITS, FLAGS and
SPECTRAL_CHANNELS extensions. Byte ranges are fetched so that only the
needed planes leave the archive.
"""
import math
import time
from urllib.request import Request, urlopen

ROOT = "https://nasa-irsa-spherex.s3.amazonaws.com/dr1/hips5_cube/hips_20261026/Norder5"


def healpix_nest(order, ra_deg, dec_deg):
    """Nested HEALPix index of (ra, dec) in degrees (Gorski et al. 2005)."""
    nside = 1 << order
    z = math.sin(math.radians(dec_deg))
    tt = (math.radians(ra_deg) % (2 * math.pi)) / (math.pi / 2)     # in [0, 4)
    za = abs(z)
    if za <= 2.0 / 3.0:
        t1 = nside * (0.5 + tt)
        t2 = nside * z * 0.75
        jp = int(t1 - t2)
        jm = int(t1 + t2)
        ifp, ifm = jp // nside, jm // nside
        face = (ifp | 4) if ifp == ifm else (ifp if ifp < ifm else ifm + 8)
        ix = jm & (nside - 1)
        iy = nside - (jp & (nside - 1)) - 1
    else:
        ntt = min(int(tt), 3)
        tp = tt - ntt
        tmp = nside * math.sqrt(3 * (1 - za))
        jp = min(int(tp * tmp), nside - 1)
        jm = min(int((1 - tp) * tmp), nside - 1)
        if z >= 0:
            face, ix, iy = ntt, nside - jm - 1, nside - jp - 1
        else:
            face, ix, iy = ntt + 8, jp, jm
    nest = 0
    for bit in range(order):
        nest |= ((ix >> bit) & 1) << (2 * bit) | ((iy >> bit) & 1) << (2 * bit + 1)
    return face * nside * nside + nest


def tile_url(npix):
    return f"{ROOT}/Dir{(npix // 10000) * 10000}/Npix{npix}_cube.fits"


def ranged(url, lo, hi, tries=4, allow_short=False):
    for attempt in range(tries):
        try:
            with urlopen(Request(url, headers={"Range": f"bytes={lo}-{hi}"}), timeout=300) as r:
                data = r.read()
            if len(data) != hi - lo + 1 and not (allow_short and data):
                raise IOError(f"short read {len(data)} of {hi - lo + 1}")
            return data
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(5 * (attempt + 1))


def read_header(url, offset):
    """Return (cards dict, header bytes, raw header text) for the HDU starting at offset."""
    chunk, text = b"", ""
    while True:
        block = ranged(url, offset + len(chunk), offset + len(chunk) + 28799, allow_short=True)
        chunk += block
        text = chunk.decode("ascii", "replace")
        end = next((i for i in range(0, len(text), 80) if text[i:i + 8] == "END     "), None)
        if end is not None:
            size = ((end + 80 + 2879) // 2880) * 2880
            cards = {}
            for i in range(0, end, 80):
                card = text[i:i + 80]
                if card[8:10] == "= ":
                    value = card[10:].split(" /")[0].strip().strip("'").strip()
                    cards[card[:8].strip()] = value
            return cards, size, text[:end + 80]


def hdu_layout(url):
    """Offsets, sizes and shapes of every HDU in a tile, read from headers only."""
    out, offset = [], 0
    for _ in range(5):
        try:
            cards, size, text = read_header(url, offset)
        except Exception:
            break
        naxis = int(cards.get("NAXIS", 0))
        shape = [int(cards[f"NAXIS{i}"]) for i in range(1, naxis + 1)]
        bitpix = int(cards["BITPIX"])
        nbytes = (abs(bitpix) // 8 * math.prod(shape) if shape else 0) + int(cards.get("PCOUNT", 0) or 0)
        out.append({"name": cards.get("EXTNAME", "PRIMARY"), "header_offset": offset, "data_offset": offset + size,
                    "data_bytes": nbytes, "shape": shape, "bitpix": bitpix, "header_text": text})
        offset += size + ((nbytes + 2879) // 2880) * 2880
        if cards.get("EXTNAME") == "SPECTRAL_CHANNELS":
            break
    return out
