"""Read and hash public survey files for the extended KMT sightlines."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
KMT = HERE.parent
sys.path[:0] = [str(KMT), str(KMT.parent), str(KMT.parent / 'ogle')]
from freeze import get
from fit import engine
from source_size import cadence_colour

PORTAL = 'https://kmtnet.kasi.re.kr/search/?'


def read(cache, tag, url):
    file = cache / (tag + '.cache')
    meta = cache / (tag + '.receipt')
    if not (file.exists() and meta.exists()):
        buf, rec = get(url)
        file.write_bytes(buf)
        meta.write_text(json.dumps(rec))
    buf, rec = file.read_bytes(), json.loads(meta.read_text())
    if hashlib.sha256(buf).hexdigest() != rec['sha256'] or rec['url'] != url:
        raise ValueError('cached URL or digest mismatch: '+tag)
    return buf, rec


def median_abs_std(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


