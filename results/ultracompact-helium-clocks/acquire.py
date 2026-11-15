"""Download the exact source coordinates, observing sector, and passband metadata."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile
from astropy.io import fits
import numpy as np

def array_digest(path):
    with fits.open(path) as h:
        d=h[1].data
        b=sha256()
        for key,dtype in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4')]:
            b.update(np.ascontiguousarray(d[key],dtype=dtype).tobytes())
        return b.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--target-dir', default='inputs')
    args=p.parse_args()
    data=Path(args.target_dir);data.mkdir(exist_ok=True)
    manifest=json.loads(Path('sources.json').read_text())
    for row in manifest['cutouts']:
        out=data/row['name']
        if not out.exists():
            with urlopen(row['url'],timeout=180) as r, ZipFile(BytesIO(r.read())) as z:
                out.write_bytes(z.read(next(s for s in z.namelist() if s.endswith('.fits'))))
        if array_digest(out)!=row['time_quality_flux_sha256']:
            raise RuntimeError(f'MAST arrays changed in {out}. Compare reductions before continuing.')
        print(out)
if __name__=='__main__':main()
