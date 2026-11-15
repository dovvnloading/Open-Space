"""Check an independent 2014 Chandra date for another X039 flare."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from astropy.io import fits

from analyze import RA, DEC, sky_wcs
from screen import retrieve

URL = 'https://cxc.cfa.harvard.edu/cdaftp/byobsid/4/16024/primary/acisf16024N004_evt2.fits.gz'
SHA256 = '355fde1e6377bfcedacd62e8a51bbde939d61c820ae8ebb1fca6450ecd73b18f'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('x039-cache'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    path = retrieve(args.cache/'16024.fits.gz', URL)
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA256:
        raise ValueError('Unexpected 2014 Chandra file digest')
    with fits.open(path) as hdus:
        events = hdus['EVENTS']
        header = events.header
        xy = sky_wcs(header).all_world2pix([[RA, DEC]], 1)[0]
        offaxis = np.linalg.norm(xy - [header['TCRPX11'],header['TCRPX12']]) * .492/60
        radius = max(1.75, 1.3 + .07*offaxis**2)
        table = events.data
        distance = np.hypot(table['x']-xy[0],table['y']-xy[1])*.492
        good = (table['ccd_id']==1) & (table['energy']>=350) & (table['energy']<=8000)
        source = good & (distance<radius)
        annulus = good & (distance>5) & (distance<=12)
        print(f"ObsID {header['OBS_ID']} on {header['DATE-OBS']}, {header['EXPOSURE']:.1f} s exposure")
        print(f'{source.sum()} source-aperture photons within {radius:.3f} arcsec, {annulus.sum()} in 5–12 arcsec annulus')
        print(f'Estimated background in source aperture: {annulus.sum()*radius**2/(12**2-5**2):.2f} photons')


if __name__ == '__main__':
    main()
