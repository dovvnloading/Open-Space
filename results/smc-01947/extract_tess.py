"""Regenerate the quality-zero TESS time and net-flux series from the MAST cutouts."""
import argparse
import gzip
import io
import json
import zipfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
POSITION = (12.668879641, -73.215620831)


def extract(name, record, output, cache):
    sector = name.split('-')[1].split('.')[0]
    saved = cache / f's{sector}.fits' if cache else None
    if saved and saved.exists():
        blob = saved.read_bytes()
    else:
        with urlopen(record['url'], timeout=180) as response:
            with zipfile.ZipFile(io.BytesIO(response.read())) as archive:
                files = [p for p in archive.namelist() if p.endswith('.fits')]
                if len(files) != 1:
                    raise ValueError(f'Expected one FITS file in sector {sector}')
                blob = archive.read(files[0])
    with fits.open(io.BytesIO(blob)) as h:
        tab = h[1].data
        if h[0].header['SECTOR'] != int(sector):
            raise ValueError('Unexpected sector')
        xy = WCS(h[2].header, naxis=2).world_to_pixel_values(*POSITION)
        x, y = np.rint(xy).astype(int)
        t = tab['TIME'].astype(float)
        cube = tab['FLUX'].astype(float)
        valid = (tab['QUALITY'] == 0) & np.isfinite(t) & np.isfinite(cube).all(axis=(1,2))
        cube = cube[valid]
        flux = cube[:, y, x] - np.median(cube[:, [0,0,6,6], [0,6,0,6]],axis=1)
        out = io.StringIO()
        np.savetxt(out,np.column_stack((t[valid],flux)),fmt=['%.12f','%.8f'],
                   header='btjd,net_flux',comments='',delimiter=',')
    with gzip.open(output/name,'wb') as compressed:
        compressed.write(out.getvalue().encode())
    with gzip.open(ROOT/'data'/name,'rb') as original:
        if out.getvalue().encode() != original.read():
            raise ValueError(f'Numerical input changed in sector {sector}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination',type=Path)
    parser.add_argument('--fits-dir',type=Path)
    args = parser.parse_args()
    args.destination.mkdir(parents=True,exist_ok=True)
    if args.destination.resolve() == (ROOT/'data').resolve():
        parser.error('Select a separate directory for the comparison')
    for filename, meta in json.loads((ROOT/'sources.json').read_text()).items():
        if filename.startswith('tess-'):
            extract(filename, meta, args.destination, args.fits_dir)
