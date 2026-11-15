"""Check the detector footprint of ULX-1 against the Chandra bad-pixel file."""
from hashlib import sha256
from io import BytesIO
from gzip import decompress
from urllib.request import urlopen

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys

PATH = ('https://heasarc.gsfc.nasa.gov/FTP/chandra/data/byobsid/9/23479/'
        'primary/acisf23479_000N002_bpix1.fits.gz')
HASH = 'af8abc6d8fcaca03aa8573329caa3c7ecc51b0c721f33947b8cae787337021db'
with fits.open(sys.argv[1]) as hdus:
    evt = hdus['EVENTS']
    sky = WCS(naxis=2)
    cols = [evt.columns.names.index(c)+1 for c in ('x','y')]
    for field, attr in [('TCRVL','crval'),('TCRPX','crpix'),('TCDLT','cdelt'),('TCTYP','ctype')]:
        setattr(sky.wcs, attr, [evt.header[f'{field}{i}'] for i in cols])
    x, y = sky.all_world2pix([[202.4164166667,47.2101666667]], 1)[0]
    scale = abs(sky.wcs.cdelt[0])*3600
    d = evt.data
    source = np.hypot(d['x']-x,d['y']-y)*scale <= 1.8
    assert np.all(d['ccd_id'][source]==7)
    print('CHIPX, CHIPY observed footprint',[(min(d[c][source]),max(d[c][source])) for c in ('chipx','chipy')])
with urlopen(PATH, timeout=60) as resp:
    raw = resp.read()
assert sha256(raw).hexdigest() == HASH
with fits.open(BytesIO(decompress(raw))) as hdus:
    regions = next(h.data for h in hdus if h.header.get('CCD_ID')==7)
    hits = []
    for row in regions:
        x, y = row['CHIPX'], row['CHIPY']
        if row['SHAPE']=='point':
            x, y = x[:1], y[:1]
        elif row['SHAPE']!='rectangle':
            raise ValueError(row['SHAPE'])
        if min(x)<=425 and max(x)>=370 and min(y)<=736 and max(y)>=680:
            hits.append((x.tolist(),y.tolist()))
    print('bad regions overlapping expanded detector box',hits)
