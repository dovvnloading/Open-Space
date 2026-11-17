"""Extract dated S/2025 U1 and Uranus states from the JPL URA184 SPK."""
import argparse
import hashlib
import re
from pathlib import Path
from urllib.request import Request, urlopen
import numpy as np
from jplephem.daf import DAF
from jplephem.spk import SPK
from jplephem.excerpter import write_excerpt
from geometry import AU, C, direction

URL = "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/ura184.bsp"


class RangeReader:
    """Bounded, verified, cached HTTP ranges of the 4.4 GB kernel."""
    def __init__(self, cache):
        self.pos = 0
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)

    def seek(self, pos, whence=0):
        if whence:
            raise ValueError("expected absolute seek")
        self.pos = pos

    def read(self, amount):
        lo = self.pos//65536*65536
        hi = ((self.pos+amount+65535)//65536)*65536 - 1
        name = self.cache / f"ura184-{lo}-{hi}.bin"
        if not name.exists():
            url = URL + f"?bytes={lo}-{hi}"
            with urlopen(Request(url, headers={"Range": f"bytes={lo}-{hi}"}), timeout=90) as f:
                hdr = f.headers.get("Content-Range", "")
                match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", hdr)
                if (not match or f.status != 206 or int(match[1]) != lo
                        or int(match[2]) != min(hi, int(match[3])-1)) :
                    raise IOError(f"invalid kernel range {hdr}")
                data = f.read(hi-lo+2)
                if len(data)!=int(match[2])-lo+1:
                    raise IOError("short kernel range")
            tmp = name.with_suffix(".tmp")
            tmp.write_bytes(data)
            tmp.replace(name)
        with name.open("rb") as f:
            f.seek(self.pos-lo)
            value = f.read(amount)
        self.pos += amount
        return value


def excerpt(start_jd, end_jd, output, cache, check=None):
    src = SPK(DAF(RangeReader(cache)))
    # Excerpt only the active segment for each body. Limit the requested span to
    # one side of a segment boundary; handle larger intervals in separate runs.
    requested = {799,75052} | ({check} if check is not None else set())
    summaries = []
    for name, row in src.daf.summaries():
        start, end, target, center, frame, kind = row[:6]
        if target in requested and center==7 and frame==1 and kind==2:
            if start <= (start_jd-2451545)*86400 and end >= (end_jd-2451545)*86400:
                summaries.append((name, row))
    if sorted(s[1][2] for s in summaries) != sorted(requested):
        raise ValueError("missing or ambiguous SPK segment at dates")
    with open(output,"w+b") as f:
        write_excerpt(src, f, start_jd, end_jd, summaries)
    print(hashlib.sha256(Path(output).read_bytes()).hexdigest(), Path(output).name)


def offsets(kernel, reception_tdb, planet, target=75052):
    """Apparent astrometric displacement relative to the Horizons planet center."""
    with SPK.open(str(kernel)) as spk:
        moon, body = spk[7,target], spk[7,799]
        delta = planet[:,4]*AU
        los = direction(planet[:,0],planet[:,1])
        emit = reception_tdb-delta/C/86400.
        diff = moon.compute(emit).T-body.compute(emit).T
        # Differential light time, relative to planet (typically <0.2 s).
        moon_emit = emit - (diff*los).sum(axis=1)/C/86400.
        return delta[:,None]*los + moon.compute(moon_emit).T - body.compute(emit).T


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("start_jd_tdb",type=float)
    p.add_argument("end_jd_tdb",type=float)
    p.add_argument("output")
    p.add_argument("cache")
    p.add_argument("--check-body", type=int)
    a = p.parse_args()
    excerpt(a.start_jd_tdb, a.end_jd_tdb, a.output, a.cache, a.check_body)
