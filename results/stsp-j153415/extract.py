"""Extract each stellar spectrum and same-petal sky controls from immutable DESI files."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.utils import iers

HERE = Path(__file__).resolve().parent
TARGET = 39633011049761464
BASE = 'https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/exposures'
VISITS = [(20220321, 127032, 5), (20220426, 131968, 9)]
DIGESTS = {
    'b5': '46a0018ac953d58d6b604fe4adec81791faffca4b15ef0695243b97633519059',
    'r5': 'e9a12f0ff55f7ed33f95382de3a6821b835af10a66b4f1431367dc1caefe6c01',
    'z5': 'b3bed5a675735cf74526455da7395160f71ff4e24dc69dceeba091aa9fdac7b6',
    'b9': 'e8d73f1af0e8fc90e043e18932c87159e17406b695633a015493c610be0c4873',
    'r9': '56c234095fcab545a5133185a32310c2a0e95a498ee3899cca33afce085002c9',
    'z9': 'b81f4f9891c1d4bb30a44d349bf3996a350352daf6b920b6d3680b7e1e698ec6',
}


def extract(cache):
    cache.mkdir(parents=True, exist_ok=True)
    iers.conf.auto_download = False
    observatory = EarthLocation.from_geodetic(-111.6003, 31.9634, 2120)
    star = SkyCoord(233.564885, 36.604730, unit='deg')
    rows, sources = [], []
    arrays = {}
    for night, exposure, petal in VISITS:
        for band in 'brz':
            camera = f'{band}{petal}'
            name = f'cframe-{camera}-{exposure:08}.fits.gz'
            url = f'{BASE}/{night}/{exposure:08}/{name}'
            path = cache / name
            if not path.exists():
                temp = path.with_suffix('.partial')
                urlretrieve(url, temp)
                temp.replace(path)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != DIGESTS[camera]:
                raise ValueError(f'Digest mismatch for {path}')
            sources.append({'url': url, 'sha256': digest})
            with fits.open(path) as h:
                fibers = h['FIBERMAP'].data
                row = np.flatnonzero(fibers['TARGETID'] == TARGET).item()
                sky = (fibers['OBJTYPE'] == 'SKY') & (fibers['FIBERSTATUS'] == 0)
                key = f'{exposure}_{band}_'
                wave = h['WAVELENGTH'].data.copy()
                flux = h['FLUX'].data
                ivar = h['IVAR'].data
                mask = h['MASK'].data
                arrays[key+'wave'] = wave
                for tag, arr in [('flux', flux), ('ivar', ivar), ('mask', mask)]:
                    arrays[key+tag] = arr[row].copy()
                # Retain individual sky fibers in the line regions as empirical controls.
                line_window = (
                    (abs(wave - 6564.61) < 22) | (abs(wave - 4862.68) < 22)
                    | (abs(wave - 4341.68) < 22) | (abs(wave - 8191) < 33)
                )
                arrays[key+'sky_wave'] = wave[line_window]
                arrays[key+'sky_flux'] = flux[sky][:,line_window].copy()
                arrays[key+'sky_mask'] = mask[sky][:,line_window].copy()
                if band == 'r':
                    meta = h[0].header
                    mid = Time(meta['MJD-OBS'] + meta['EXPTIME']/172800.,
                               format='mjd', scale='utc', location=observatory)
                    hjd = mid.utc.jd + mid.light_travel_time(star, kind='heliocentric').to_value(u.day)
                    rows.append(dict(hjd_utc_mid=float(hjd), night=night, expid=exposure, targetid=str(TARGET),
                                     mjd_start=meta['MJD-OBS'], seconds=meta['EXPTIME'],
                                     seeing_arcsec=meta['SEEING'], heliocor=meta['HELIOCOR'],
                                     fiber=int(fibers[row]['FIBER']), status=int(fibers[row]['FIBERSTATUS']),
                                     ra=float(fibers[row]['FIBER_RA']), dec=float(fibers[row]['FIBER_DEC']),
                                     sky_fibers=int(sum(sky))))
    np.savez_compressed(HERE/'data'/'spectra.npz', **arrays)
    with (HERE/'data'/'visits.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    (HERE/'data'/'desi-products.json').write_text(json.dumps(sources, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    extract(parser.parse_args().cache)
