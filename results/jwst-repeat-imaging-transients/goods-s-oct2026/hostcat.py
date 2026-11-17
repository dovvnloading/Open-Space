from urllib.request import urlopen, Request
from astropy.io import fits
from io import BytesIO
from math import prod, ceil
import time, os

url = "https://archive.stsci.edu/hlsps/jades/dr2/goods-s/catalogs/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits"
if os.getenv("AWS"):
    url = "https://stpubdata.s3.amazonaws.com/hlsp/jades/dr2/goods-s/catalogs/hlsp_jades_jwst_nircam_goods-s-deep_photometry_v2.0_catalog.fits"


def read_once(start, num):
    r = urlopen(
        Request(url, headers={"Range": f"bytes={start}-{start+num-1}"}), timeout=35
    )
    assert r.status == 206 and r.getheader("Content-Range").startswith(
        f"bytes {start}-"
    ), (r.status, r.getheader("Content-Range"))
    return r.read(), int(r.getheader("Content-Range").split("/")[-1])


def read_range(start, num):
    for i in range(5):
        try:
            return read_once(start, num)
        except Exception:
            if i == 4:
                raise
            time.sleep(i + 1)


def header(start):
    length = 14400
    while True:
        buf, end = read_range(start, length)
        if any(buf[i : i + 8] == b"END     " for i in range(0, len(buf), 80)):
            break
        length *= 2
    b = BytesIO(buf)
    h = fits.Header.fromfile(b)
    return h, b.tell(), end


if __name__ == "__main__":
    import json

    pos = 0
    exts = []
    while True:
        h, sz, total = header(pos)
        data = (
            prod(h[f"NAXIS{x}"] for x in range(1, h["NAXIS"] + 1))
            * abs(h["BITPIX"])
            // 8
            if h["NAXIS"]
            else 0
        )
        off = pos + sz
        name = h.get("EXTNAME", "PRIMARY")
        print(
            name,
            h.get("NAXIS2", 0),
            round(data * 1e-6),
            [h[k] for k in h if k.startswith("TTYPE")][:12],
            flush=True,
        )
        exts.append(
            dict(
                name=name,
                offset=off,
                size=data,
                header=h.tostring(sep="\n", padding=False),
            )
        )
        pos = off + 2880 * ceil((data + h.get("PCOUNT", 0)) / 2880)
        json.dump(exts, open("host_exts.json", "w"))
        if pos >= total:
            break
    json.dump(exts, open("host_exts.json", "w"))
