"""Streaming I/O for the fixed-width FITS binary tables in the two cited releases."""

import csv
import io
import re
import time
import urllib.request
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import hashlib
import http.client
import numpy as np

BH_CONTROLS = {4373465352415301632, 5870569352746779008, 4318465066420528000}


def target_ids(path):
    with open(path) as stream:
        ids = {int(row["source_id"]) for row in csv.DictReader(stream)}
    return np.array(sorted(ids | BH_CONTROLS), dtype=np.int64)


def read_header(stream):
    header = {}
    while True:
        block = stream.read(2880)
        if not block and not header:
            return None
        if len(block) != 2880:
            raise ValueError("Incomplete FITS header")
        for start in range(0, 2880, 80):
            card = block[start : start + 80]
            key = card[:8].decode("ascii").strip()
            if key == "END":
                return header
            if card[8:10] == b"= ":
                header[key] = (
                    card[10:].decode("ascii").split("/")[0].strip().strip("'").strip()
                )


def table_dtype(header):
    types = {
        "D": ">f8",
        "E": ">f4",
        "K": ">i8",
        "J": ">i4",
        "I": ">i2",
        "B": "u1",
        "L": "S1",
    }
    fields = []
    for i in range(1, int(header["TFIELDS"]) + 1):
        repeat, code = re.fullmatch(r"(\d*)([DEKJIBLA])", header[f"TFORM{i}"]).groups()
        count = int(repeat or 1)
        name = header[f"TTYPE{i}"]
        if code == "A":
            fields.append((name, f"S{count}"))
        elif count == 1:
            fields.append((name, types[code]))
        else:
            fields.append((name, types[code], (count,)))
    dtype = np.dtype(fields)
    if dtype.itemsize != int(header["NAXIS1"]):
        raise ValueError("Table width disagrees with FITS header")
    return dtype


def byte_range(url, first, last):
    """Require the requested interval, including the archive's reported full size."""
    for retry in range(4):
        try:
            request = urllib.request.Request(
                url, headers={"Range": f"bytes={first}-{last}"}
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                match = re.fullmatch(
                    r"bytes (\d+)-(\d+)/(\d+)", response.getheader("Content-Range", "")
                )
                if response.status != 206 or match is None:
                    raise OSError("Archive did not honor the requested HTTP Range")
                start, end, total = map(int, match.groups())
                data = response.read()
                if (
                    start != first
                    or end != min(last, total - 1)
                    or len(data) != end - start + 1
                ):
                    raise OSError("Wrong archive interval or incomplete table")
                return data, total
        except (OSError, http.client.IncompleteRead):
            if retry == 3:
                raise
            time.sleep(2**retry)


class RangeStream(io.RawIOBase):
    """Ordered hash and read stream with a bounded number of parallel HTTP reads."""

    def __init__(self, url, length, chunk_size=4 * 1024**2, workers=7):
        self.url, self.length, self.chunk_size = url, length, chunk_size
        self.digest = hashlib.sha256()
        self.seen = 0
        self.starts = iter(range(0, length, chunk_size))
        self.pool = ThreadPoolExecutor(workers)
        self.pending = deque()
        self.buffer = memoryview(b"")
        for _ in range(workers):
            self.enqueue()

    def enqueue(self):
        start = next(self.starts, None)
        if start is not None:
            self.pending.append(
                self.pool.submit(
                    byte_range,
                    self.url,
                    start,
                    min(start + self.chunk_size, self.length) - 1,
                )
            )

    def readable(self):
        return True

    def readinto(self, dest):
        if not self.buffer:
            if not self.pending:
                return 0
            chunk, total = self.pending.popleft().result()
            if total != self.length:
                raise OSError("Product size differs from the pinned release")
            self.enqueue()
            self.digest.update(chunk)
            self.seen += len(chunk)
            self.buffer = memoryview(chunk)
        count = min(len(dest), len(self.buffer))
        dest[:count] = self.buffer[:count]
        self.buffer = self.buffer[count:]
        return count

    def close(self):
        self.pool.shutdown(wait=True, cancel_futures=True)
        super().close()
