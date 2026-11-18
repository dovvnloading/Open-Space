"""Seekable PDS HTTP byte ranges with an exact record of bytes read."""

import hashlib
import os
from pathlib import Path
from urllib.request import Request, urlopen


class RangeSource:
    chunk_size = 8 * 1024 * 1024

    def __init__(self, url, size, cache):
        self.url = url
        self.size = size
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.position = 0
        self.chunks = {}
        self.current = None

    def seek(self, offset, whence=os.SEEK_SET):
        origins = {os.SEEK_SET: 0, os.SEEK_CUR: self.position, os.SEEK_END: self.size}
        position = origins[whence] + offset
        if position < 0:
            raise ValueError("Negative video byte position")
        self.position = position
        return position

    def tell(self):
        return self.position

    def seekable(self):
        return True

    def readable(self):
        return True

    def block(self, start):
        if self.current is not None and self.current[0] == start:
            return self.current[1]
        stop = min(self.size, start + self.chunk_size) - 1
        length = stop - start + 1
        path = self.cache / f"{start:012d}-{stop:012d}.bin"
        if path.is_file():
            data = path.read_bytes()
        else:
            request = Request(self.url, headers={"Range": f"bytes={start}-{stop}",
                                                 "User-Agent": "Open-Space-research/1.0"})
            with urlopen(request, timeout=40) as response:
                if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{stop}/{self.size}":
                    raise IOError(f"PDS response did not match requested byte range: {start}-{stop}")
                data = response.read(length + 1)
            if len(data) != length:
                raise IOError(f"Incomplete PDS byte range: {start}-{stop}")
            path.write_bytes(data)
        if len(data) != length:
            raise IOError(f"Invalid cached byte range: {path.name}")
        self.chunks[start] = {"start": start, "stop_inclusive": stop,
                              "sha256": hashlib.sha256(data).hexdigest()}
        self.current = (start, data)
        return data

    def read(self, size):
        end = min(self.position + size, self.size)
        pieces = []
        while self.position < end:
            start = self.position // self.chunk_size * self.chunk_size
            block = self.block(start)
            start_inside = self.position - start
            part = block[start_inside : min(end-start, len(block))]
            pieces.append(part)
            self.position += len(part)
        return b"".join(pieces)

    def receipt(self):
        return {"url": self.url, "source_bytes": self.size,
                "ranges": [self.chunks[key] for key in sorted(self.chunks)]}
