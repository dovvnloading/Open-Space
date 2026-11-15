"""Compute and apply the published JPL EFCC18 position and proper-motion corrections."""
import hashlib
import io
import tarfile

import astropy.units as u
import numpy as np
from astropy.time import Time
from astropy_healpix import HEALPix

ARCHIVE_SHA = '288132cefea7f1bdf0989782db1173948dd7dfef7da31910be90b16b73cfea6c'
MAPPED = {'USNOA2': 'c', 'USNOB1': 'o', 'UCAC2': 'r', 'UCAC4': 'q'}
SKIP = {'Gaia1', 'Gaia2', 'Gaia3E', 'Gaia3', 'ATLAS2'}


def corrected(rows, archive):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA:
        raise ValueError('EFCC18 source digest differs')
    with tarfile.open(archive, 'r:gz') as bundle:
        content = bundle.extractfile('bias.dat').read()
    lines = content.decode().splitlines()
    codes = lines[4].removeprefix('!').split()
    dataset = np.loadtxt(io.BytesIO(content), comments='!').reshape(49152, len(codes), 4)
    hp = HEALPix(nside=64, order='ring')
    corrections = []
    result = []
    for original in rows:
        row = original.copy()
        result.append(row)
        if row['catalog'] in SKIP:
            continue
        cat = MAPPED[row['catalog']]
        ra, dec = float(row['ra_deg']), float(row['dec_deg'])
        pixel = hp.lonlat_to_healpix(ra*u.deg, dec*u.deg)
        years = Time(row['utc'].removesuffix('Z'), scale='utc').tt.jyear - 2000.
        dra, ddec, pmra, pmdec = dataset[pixel, codes.index(cat)]
        east, north = dra + years * pmra / 1000, ddec + years * pmdec / 1000
        row['ra_deg'] = str(ra - east / (3600*np.cos(np.radians(dec))))
        row['dec_deg'] = str(dec - north / 3600)
        corrections.append({'obsid': row['obsid'], 'catalog': row['catalog'], 'pixel': int(pixel),
                            'subtracted_east_arcsec': round(east, 7),
                            'subtracted_north_arcsec': round(north, 7)})
    return result, corrections
