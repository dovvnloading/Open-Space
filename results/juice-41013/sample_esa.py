"""Recreate the released half-hour JUICE state vectors from ESA SPICE kernels."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import spiceypy as spice


def write_states(directory, output):
    manifest = json.loads((Path(__file__).parent/'data/sources.json').read_text())
    writer = csv.writer(output)
    writer.writerow(['revision','jd_tdb','x_km','y_km','z_km','vx_km_s','vy_km_s','vz_km_s'])
    for filename, receipt in manifest.items():
        if not filename.endswith('.bsp'):
            continue
        path = directory / filename
        with path.open('rb') as f:
            digest = hashlib.sha256(f.read()).hexdigest()
        if digest != receipt['sha256']:
            raise ValueError(f'Unrecognized kernel bytes for {filename}')
        revision = 'crema5' if 'crema' in filename else filename.split('_')[2]
        spice.furnsh(str(path))
        try:
            for step in range(-4,38):
                jd = 2461735.0 + (step/2)/24
                vector, _ = spice.spkezr('JUICE',(jd-2451545)*86400,'J2000','NONE','SUN')
                writer.writerow([revision,repr(jd),*[format(v,'.16g') for v in vector]])
        finally:
            spice.unload(str(path))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kernels',type=Path)
    write_states(parser.parse_args().kernels,sys.stdout)
