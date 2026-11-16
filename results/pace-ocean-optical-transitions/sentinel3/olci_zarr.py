"""Small read-only reader for the two Copernicus Marine Zarr v2 arrays used here."""

import ctypes
import ctypes.util
import datetime as dt
import hashlib
import json
import os
import tempfile
from pathlib import Path
from time import sleep
from urllib.error import HTTPError
from urllib.request import urlopen

import numpy as np

DATASET = "cmems_obs-oc_glo_bgc-reflectance_my_l3-olci-4km_P1D_202207"
STAC = (
    "https://stac.marine.copernicus.eu/metadata/"
    f"OCEANCOLOUR_GLO_BGC_L3_MY_009_103/{DATASET}/dataset.stac.json"
)
ROOT = Path(__file__).resolve().parent
CACHE = Path(
    os.environ.get("OLCI_CACHE", Path(tempfile.gettempdir()) / "mission06-olci")
)
BLOSC = ctypes.CDLL(ctypes.util.find_library("blosc") or "libblosc.so.1")
BLOSC.blosc_decompress.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t]
RECEIPTS = {}


def read(url):
    """Cache responses by URL and record the content digest for offline comparison."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / hashlib.sha256(url.encode()).hexdigest()
    if path.exists():
        data = path.read_bytes()
    else:
        for attempt in range(3):
            try:
                with urlopen(url, timeout=45) as response:
                    data = response.read()
                with tempfile.NamedTemporaryFile(dir=CACHE, delete=False) as tmp:
                    tmp.write(data)
                os.replace(tmp.name, path)
                break
            except OSError as exc:
                if attempt == 2 or isinstance(exc, HTTPError) and exc.code == 404:
                    raise
                sleep(1 + attempt)
    RECEIPTS[url] = (len(data), hashlib.sha256(data).hexdigest())
    return data


class Grid:
    def __init__(self, layout):
        self.base = json.loads(read(STAC))["assets"][layout]["href"]
        self.meta = json.loads(read(self.base + "/.zmetadata"))["metadata"]
        self.lat = self.chunk("latitude", [0])
        self.lon = self.chunk("longitude", [0])
        meta = self.meta["time/.zarray"]
        (ntime,) = meta["shape"]
        (block,) = meta["chunks"]
        days = np.concatenate(
            [self.chunk("time", [i]) for i in range((ntime - 1) // block + 1)]
        )[:ntime]
        if self.meta["time/.zattrs"]["units"] != "days since 1900-01-01":
            raise ValueError("Unrecognized time units")
        self.dates = {
            dt.date(1900, 1, 1) + dt.timedelta(days=int(d)): i
            for i, d in enumerate(days)
        }

    def chunk(self, name, indices):
        meta = self.meta[name + "/.zarray"]
        key = meta.get("dimension_separator", ".").join(map(str, indices))
        url = f"{self.base}/{name}/{key}"
        try:
            blob = read(url)
        except HTTPError as exc:
            if exc.code != 404:
                raise
            RECEIPTS[url] = (0, "HTTP_404_unstored_fill_chunk")
            return np.full(meta["chunks"], meta["fill_value"], dtype=meta["dtype"])
        data = np.empty(meta["chunks"], dtype=meta["dtype"])
        if (
            meta["compressor"]["id"] != "blosc"
            or meta["order"] != "C"
            or meta.get("filters")
        ):
            raise ValueError(f"Unsupported data encoding at {url}")
        if (
            len(blob) < 16
            or int.from_bytes(blob[4:8], "little") != data.nbytes
            or int.from_bytes(blob[12:16], "little") != len(blob)
        ):
            raise ValueError(f"Truncated or incompatible compressed block at {url}")
        size = BLOSC.blosc_decompress(blob, data.ctypes.data, data.nbytes)
        if size != data.nbytes:
            raise ValueError(f"Decompression failed at {url}")
        return data

    def bounds(self, south, north, west, east):
        return (
            np.flatnonzero((self.lat >= south) & (self.lat <= north)),
            np.flatnonzero((self.lon >= west) & (self.lon <= east)),
        )

    def extract(self, name, date, y, x):
        ti = self.dates[date]
        meta = self.meta[name + "/.zarray"]
        nt, ny, nx = meta["chunks"]
        result = np.full((len(y), len(x)), np.nan)
        for a in np.unique(y // ny):
            rows = y // ny == a
            for b in np.unique(x // nx):
                cols = x // nx == b
                chunk = self.chunk(name, [ti // nt, a, b])[ti % nt]
                values = chunk[np.ix_(y[rows] % ny, x[cols] % nx)].astype(float)
                values[values == meta["fill_value"]] = np.nan
                result[np.ix_(rows, cols)] = values
        return result
