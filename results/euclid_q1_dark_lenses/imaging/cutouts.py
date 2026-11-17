"""Read selected products from the immutable Q1 public ZIP using HTTP ranges."""
import io
import math
from pathlib import Path
import struct
import time
import urllib.request
import zipfile
import zlib

import numpy as np

RECORD = "https://zenodo.org/api/records/15025832/files"
# Archive sizes and whole-file MD5 are published with the record. Partial downloads
# verify Content-Range, and each decompressed member's length and ZIP CRC32.
ARCHIVES = {"lens": (3051362622, "fc09586703c560c2af55dd50d6490c8b"),
            "group": (255134676, "22e5bad271c411bf27eef58887cd612b"),
            "unsuccess": (315674226, "681047ab3238eae7b3eae93417f0bccd"),
            "recenter": (56358693, "111faa4ee1ff1d16a5531539a8768986")}


def read_range(name, first, last):
    size = ARCHIVES[name][0]
    url = f"{RECORD}/{name}.zip/content"
    req = urllib.request.Request(url, headers={"Range": f"bytes={first}-{last}", "Accept-Encoding": "identity"})
    for retry in range(5):
        try:
            with urllib.request.urlopen(req, timeout=55) as response:
                expected = f"bytes {first}-{last}/{size}"
                if response.status != 206 or response.headers.get("Content-Range") != expected:
                    raise ValueError(f"Unexpected Content-Range, expected {expected}")
                b = response.read(last-first+2)
                if len(b) != last-first+1:
                    raise IOError("Truncated byte range")
                return b
        except (OSError, ValueError):
            if retry == 4:
                raise
            time.sleep(2**retry)


class CutoutArchive:
    def __init__(self, name, cache=None):
        self.name = name
        size = ARCHIVES[name][0]
        index = Path(cache, name+".index") if cache else None
        if index and index.exists():
            central = index.read_bytes()
        else:
            end = read_range(name, size-65557, size-1)
            pos = end.rfind(b"PK\x05\x06")
            _, _, _, _, _, cdsize, start, comment = struct.unpack_from("<4s4H2IH", end, pos)
            if size-start != cdsize+22+comment:
                raise ValueError("Unrecognized ZIP directory")
            central = read_range(name, start, size-1)
            if index:
                index.parent.mkdir(parents=True, exist_ok=True)
                index.write_bytes(central)
        self.offset = size-len(central)
        self.entries = {e.filename: e for e in zipfile.ZipFile(io.BytesIO(central)).infolist()}

    def members(self, obj, products, cached_chunk=None):
        names = [f"{self.name}/{obj}/{name}" for name in products]
        entries = [self.entries[n] for n in names if n in self.entries]
        if not entries:
            raise KeyError(obj)
        start = self.offset+min(e.header_offset for e in entries)
        last = max(entries, key=lambda e:e.header_offset)
        stop = min(ARCHIVES[self.name][0], self.offset+last.header_offset+last.compress_size+4096)
        buf = read_range(self.name, start, stop-1) if cached_chunk is None else cached_chunk[start:stop]
        output = {}
        for e in entries:
            pos = self.offset+e.header_offset-start
            if buf[pos:pos+4] != b"PK\x03\x04":
                raise ValueError("Invalid member signature")
            n, extra = struct.unpack_from("<HH", buf, pos+26)
            if buf[pos+30:pos+30+n].decode() != e.filename:
                raise ValueError("Mismatched ZIP member")
            pos += 30+n+extra
            compressed = buf[pos:pos+e.compress_size]
            if e.compress_type not in (0, 8):
                raise ValueError("Unsupported compression")
            data = zlib.decompress(compressed,-15) if e.compress_type==8 else compressed
            if zlib.crc32(data)!=e.CRC or len(data)!=e.file_size:
                raise ValueError(f"Member failed CRC or byte count: {e.filename}")
            output[e.filename.split("/",2)[2]] = data
        return output


def image_extensions(raw):
    """Read the simple, unscaled, uncompressed image HDUs in this specific release."""
    pos, images = 0, {}
    while pos < len(raw):
        h, off = {}, pos
        while off+80 <= len(raw):
            c = raw[off:off+80].decode("ascii")
            off += 80
            if c[:8].strip() == "END":
                break
            if c[8:10] == "= ":
                val = c[10:].strip()
                h[c[:8].strip()] = val.split("'")[1].strip() if val.startswith("'") else val.split("/")[0].strip()
        else:
            raise ValueError("Missing FITS END")
        off = pos + ((off-pos+2879)//2880)*2880
        naxis = int(h["NAXIS"])
        if naxis==0:
            pos=off
            continue
        if naxis!=2 or h.get("BSCALE","1")!="1" or h.get("BZERO","0")!="0" or h.get("XTENSION","IMAGE").strip()!="IMAGE":
            raise ValueError("Unexpected FITS layout")
        dtype = {"-64":">f8","-32":">f4","8":"u1","16":">i2","32":">i4"}[h["BITPIX"]]
        shape=(int(h["NAXIS2"]),int(h["NAXIS1"]))
        data=np.frombuffer(raw,dtype=dtype,count=math.prod(shape),offset=off).reshape(shape)
        images[h.get("EXTNAME","PRIMARY")] = (h,data)
        pos=off + ((data.nbytes+2879)//2880)*2880
    if pos!=len(raw):
        raise ValueError("FITS byte count mismatch")
    return images
