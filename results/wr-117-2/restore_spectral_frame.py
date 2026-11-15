"""Restore one pinned full frame and verify the pixels behind an offline measurement."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from replay import RA, DEC, ROOT


def crop(source, row):
    with fits.open(source, memmap=False) as hdus:
        hdr = hdus['IMAGE'].header
        if hdr['OBSID'] != row['obs_id'] or hdr['DETECTOR'] != int(row['detector']):
            raise ValueError('Archive identifier mismatch')
        x, y = WCS(hdr, hdus).all_world2pix([[RA, DEC]], 0)[0]
        ny, nx = hdus['IMAGE'].shape
        i, j = int(np.floor(x)), int(np.floor(y))
        x0, x1 = max(0, i-10), min(nx, i+12)
        y0, y1 = max(0, j-10), min(ny, j+12)
        result = fits.HDUList([fits.PrimaryHDU(header=hdus[0].header.copy())])
        for name in ('IMAGE', 'FLAGS', 'VARIANCE'):
            header = hdus[name].header.copy()
            for axis, shift in [(1, x0), (2, y0)]:
                for suffix in ('', 'A', 'W'):
                    key = f'CRPIX{axis}{suffix}'
                    if key in header:
                        header[key] -= shift
            result.append(fits.ImageHDU(hdus[name].section[y0:y1, x0:x1], header=header, name=name))
        result.append(hdus['WCS-WAVE'].copy())
        stream = io.BytesIO()
        result.writeto(stream, output_verify='silentfix')
        return stream.getvalue()


def restore(detector, obs, destination):
    with (ROOT/'spectral-inputs.csv').open() as handle:
        row = next((r for r in csv.DictReader(handle)
                    if r['obs_id'] == obs and int(r['detector']) == detector), None)
    if row is None:
        raise ValueError(f'No archived frame for {obs} D{detector}')
    with urlopen(row['aws_url'], timeout=120) as response:
        raw = response.read()
    if hashlib.sha256(raw).hexdigest() != row['full_sha256']:
        raise ValueError('Full-frame checksum mismatch')
    stamp = crop(io.BytesIO(raw), row)
    if hashlib.sha256(stamp).hexdigest() != row['cutout_sha256']:
        raise ValueError('Cropped FITS bytes differ; use Astropy 7.2.3 for byte-exact writing')
    destination.write_bytes(stamp)


if __name__ == '__main__':
    from astropy import log
    log.setLevel('ERROR')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('detector', type=int, choices=range(1, 6))
    parser.add_argument('obs_id')
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    restore(args.detector, args.obs_id, args.destination)
