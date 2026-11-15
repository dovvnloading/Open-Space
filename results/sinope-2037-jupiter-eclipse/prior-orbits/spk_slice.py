"""Read the selected type-2 NAIF SPK records, without downloading the full kernel."""

import struct
from pathlib import Path
from urllib.request import Request, urlopen
import numpy as np

BASE = (
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/"
)
J0 = 2451545.0


class Kernel:
    def __init__(self, key, cache):
        self.key = key
        self.used = set()
        self.url = BASE + key + ".bsp"
        self.cache = Path(cache)
        self.cache.mkdir(exist_ok=True, parents=True)
        record = self.read(0, 1024)
        assert record[:8] == b"DAF/SPK "
        self.order = {b"LTL-IEEE": "<", b"BIG-IEEE": ">"}[record[88:96]]
        self.first = struct.unpack_from(self.order + "I", record, 76)[0]

    def read(self, start, n):
        self.used.add((start, n))
        path = self.cache / f"{self.key}-{start}-{n}.bin"
        if not path.exists():
            end = start + n - 1
            with urlopen(
                Request(self.url, headers={"Range": f"bytes={start}-{end}"}), timeout=90
            ) as response:
                if response.status != 206 or not response.headers[
                    "Content-Range"
                ].startswith(f"bytes {start}-{end}/"):
                    raise ValueError(
                        f'Invalid range response for {start}-{end}: {response.status} {response.headers["Content-Range"]}'
                    )
                raw = response.read()
            if len(raw) != n:
                raise ValueError("Incomplete kernel range")
            path.write_bytes(raw)
        return path.read_bytes()

    def segment(self, target, epoch):
        rec = self.first
        while rec:
            raw = self.read((rec - 1) * 1024, 1024)
            nxt, _, n = struct.unpack_from(self.order + "3d", raw)
            for j in range(int(n)):
                lo, hi, body, center, frame, kind, first, last = struct.unpack_from(
                    self.order + "2d6i", raw, 24 + j * 40
                )
                if body == target and lo <= epoch <= hi:
                    return center, frame, kind, first, last
            rec = int(nxt)
        raise ValueError(f"No segment for {target}")

    def coefficients(self, target, epoch):
        center, frame, kind, first, last = self.segment(target, epoch)
        epoch0, span, record_len, count = np.frombuffer(
            self.read((last - 4) * 8, 32), dtype=self.order + "f8"
        )
        ind = min(int((epoch - epoch0) // span), int(count) - 1)
        n = int(record_len)
        addr = (first - 1 + ind * n) * 8
        coefficients = np.frombuffer(self.read(addr, n * 8), dtype=self.order + "f8")
        return center, frame, kind, addr, coefficients


def evaluate(record, seconds):
    midpoint, radius = record[:2]
    x = (seconds - midpoint) / radius
    if abs(x) > 1.0000001:
        raise ValueError("Outside record")
    return np.array(
        [np.polynomial.chebyshev.chebval(x, co) for co in record[2:].reshape(3, -1)]
    )
