"""Extract a lossless spatial subset of the public SPHEREx DR1 all-sky tile."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
TILE_URL = ('https://nasa-irsa-spherex.s3.amazonaws.com/dr1/hips5_cube/'
            'hips_20261026/Norder5/Dir10000/Npix10542_cube.fits')
TARGET = (242.61055777065044, -51.357039296633985)


def main(path):
    digest = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    if digest != 'cbd76bc310e4063b30f66205d097f570c832d2d4982bd97a8b4469eede1e3e89':
        raise ValueError('Unexpected DR1 cube SHA-256: ' + digest)
    with fits.open(path, memmap=False) as hdus:
        wcs = WCS(hdus['IMAGE'].header, hdus, naxis=[1, 2]); radius = 19
        x, y = np.rint(wcs.all_world2pix([TARGET], 0)[0]).astype(int)
        x0, y0 = x - radius, y - radius
        cut = np.s_[:, y0:y+radius+1, x0:x+radius+1]
        spatial = wcs.to_header()
        spatial['CRPIX1'] -= x0
        spatial['CRPIX2'] -= y0
        channels = hdus['SPECTRAL_CHANNELS'].data
        np.savez_compressed(ROOT/'inputs'/'spherex.npz', image=hdus['IMAGE'].data[cut],
                            hits=hdus['NHITS'].data[cut], flags=hdus['FLAGS'].data[cut],
                            wavelength=channels['WAVELENGTH'], header=spatial.tostring())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tile', type=Path, help='Original Npix10542_cube.fits')
    main(parser.parse_args().tile)
