"""Restore the large DE440s kernel and the two external comparison data products."""

import hashlib
import json
import urllib.request
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "work"
KERNEL_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp"
PRIOR_URL = "https://raw.githubusercontent.com/katvolk/TNO-classifications/4b1b081fa0b847034a94c1b154e971506629b36b/december-2023-allTNO-classifications.csv"


def save(url, name):
    dest = WORK / name
    if not dest.exists():
        with urllib.request.urlopen(url, timeout=120) as fp:
            dest.write_bytes(fp.read())
    return dest


if __name__ == "__main__":
    for url, name in [(KERNEL_URL, "de440s.bsp"), (PRIOR_URL, "volk.csv")]:
        f = save(url, name)
        expected = json.loads((ROOT / "checksums.json").read_text())[name]
        if hashlib.sha256(f.read_bytes()).hexdigest() != expected:
            raise SystemExit("Unexpected bytes for " + name)
