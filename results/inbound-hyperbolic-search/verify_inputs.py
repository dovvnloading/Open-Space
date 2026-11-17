"""Verify deposited tracklet and decompressed comet-observation receipt bytes."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def verify():
    inputs = json.loads((ROOT / 'tracklets/run/input.json').read_text())
    slices = {'f51.mpc': 'ps1_sha256', 'f51-all.mpc': 'ps1_all_sha256', 'x05.mpc': 'selected_sha256'}
    checked = 0
    for name, key in slices.items():
        data = (ROOT / 'tracklets/run' / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != inputs[key]:
            raise ValueError(f'Tracklet receipt mismatch: {name}')
        checked += 1
    receipt = json.loads((ROOT / 'orbits/data/acquisition.json').read_text())
    for name, item in receipt['objects'].items():
        path = ROOT / 'orbits/data' / (name + '.obs.gz')
        data = gzip.decompress(path.read_bytes())
        if hashlib.sha256(data).hexdigest() != item['sha256']:
            raise ValueError(f'Observation receipt mismatch: {name}')
        checked += 1
    return checked


if __name__ == '__main__':
    print(f'{verify()} frozen observation inputs verified')
