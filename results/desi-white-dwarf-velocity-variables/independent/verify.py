"""Run next to the bundled inputs. Requires numpy, pandas, scipy and astropy."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import LombScargle
from astropy.wcs import WCS

parser = argparse.ArgumentParser()
parser.add_argument('--scan', action='store_true', help='Repeat the full blind g/r frequency scan.')
args = parser.parse_args()
here = Path(__file__).resolve().parent
manifest = json.loads((here / 'manifest.json').read_text())
for filename, product in manifest.items():
    path = here / filename
    packed = path.with_suffix(path.suffix + '.gz')
    if path.exists():
        raw = path.read_bytes()
    elif packed.exists():
        raw = gzip.decompress(packed.read_bytes())
    else:
        raw = urlopen(product['url'], timeout=150).read()
        path.write_bytes(raw)
    if hashlib.sha256(raw).hexdigest() != product['sha256']:
        raise ValueError(f'Input changed: {filename}. Use the frozen DR24 input.')

for name, frequency in [('0421', 54.323502), ('0238', 47.353499)]:
    path = here / f'{name}_ztf.csv'
    table = pd.read_csv(path if path.exists() else path.with_suffix('.csv.gz'))
    clean = table.query('catflags == 0 and .003 < magerr < .15 and chi < 4 and abs(sharp) < .5 and airmass < 2.5')
    clean = clean.sort_values('magerr').drop_duplicates(['expid', 'filtercode']).copy()
    clean['residual'] = clean.mag - clean.groupby('oid').mag.transform('median')
    clean['error'] = np.hypot(clean.magerr, .01)
    gr = clean[clean.filtercode.isin(['zg', 'zr'])]
    if args.scan:
        search = LombScargle(gr.hjd.values - 2458000, gr.residual.values, gr.error.values)
        freqs, power = search.autopower(minimum_frequency=.2, maximum_frequency=144, samples_per_peak=7, method='fast')
        index = np.argmax(power)
        fap = search.false_alarm_probability(power[index], minimum_frequency=.2, maximum_frequency=144)
        print(name, 'frequency/day, power, Baluev FAP', freqs[index], power[index], fap)
    for key, group in [('g', gr[gr.filtercode == 'zg']), ('r', gr[gr.filtercode == 'zr'])]:
        phase = 2 * np.pi * frequency * (group.hjd.values - 2458000)
        design = np.stack([np.ones(len(group)), np.cos(phase), np.sin(phase)], axis=1)
        weighted = design / group.error.values[:, None]
        beta = np.linalg.lstsq(weighted, group.residual.values / group.error.values, rcond=None)[0]
        residual = (group.residual.values - design @ beta) / group.error.values
        covariance = np.linalg.inv(weighted.T @ weighted) * (residual @ residual) / (len(group) - 3)
        amplitude = np.linalg.norm(beta[1:])
        uncertainty = np.sqrt(beta[1:] @ covariance[1:, 1:] @ beta[1:]) / amplitude
        print(name, key, 'N, semi-amplitude and local error [mmag]', len(group), amplitude * 1000, uncertainty * 1000)

header = fits.getheader(here / 'lamost-257903093.fits.gz')
print('LAMOST', {k: header[k] for k in ['DATE-OBS', 'NEXP', 'EXPTIME', 'SNRG', 'Z', 'Z_ERR', 'SUBCLASS', 'HELIO', 'VACUUM']})
print('LAMOST cz, error [km/s]', 299792.458 * header['Z'], 299792.458 * header['Z_ERR'])

filename = here / 'ztf-139326751.fits'
data = fits.getdata(filename)
wcs = WCS(fits.getheader(filename), naxis=2)
center = np.array(wcs.world_to_pixel_values(39.539388, 27.557849))
y, x = np.indices(data.shape)
selected = (np.hypot(x-center[0], y-center[1]) > 6) & (data-np.nanmedian(data) > 100)
points = np.column_stack([x[selected], y[selected]])
_, basis = np.linalg.eigh(np.cov(points.T))
print('trail N, length, cross-track rms, WD offset [pixel]', len(points),
      np.ptp(points @ basis[:, 1]), np.std(points @ basis[:, 0]), abs((center-points.mean(axis=0)) @ basis[:, 0]))

# Use Gaia's epoch-2016.0 proper motion to isolate the target in each PS1 exposure.
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers
iers.conf.auto_download = False
panstarrs = EarthLocation.from_geodetic(-156.2569*u.deg, 20.7075*u.deg, 3055*u.m)
objects = [('0421',65.48040660076,-13.42532402194,34.897,11.456,54.323502),
           ('0238',39.53931386723,27.55791404235,-1.229,-81.414,47.353499)]
for name, ra, dec, pmra, pmdec, freq in objects:
    data = pd.read_csv(here / f'{name}_ps1_detection.csv')
    year = 2000 + (data.obsTime-51544.5)/365.25
    east = (data.ra-ra)*np.cos(np.deg2rad(dec))*3600 - (year-2016)*pmra/1000
    north = (data.dec-dec)*3600 - (year-2016)*pmdec/1000
    good = (np.hypot(east,north)<1) & (data.psfQfPerfect>.9) & (data.psfFlux>0) & (data.psfFluxErr/data.psfFlux<.1)
    data = data[good].copy()
    data['magnitude'] = -2.5*np.log10(data.psfFlux)
    residual = (data.magnitude-data.groupby('filterID').magnitude.transform('median')).to_numpy()
    error = np.hypot(.01,1.086*data.psfFluxErr/data.psfFlux).to_numpy()
    direction = SkyCoord(ra*u.deg, dec*u.deg)
    times = Time(data.obsTime.to_numpy(), format='mjd', scale='utc', location=panstarrs)
    hjd = times.utc.jd + times.light_travel_time(direction,kind='heliocentric').to_value('d')
    phase = 2*np.pi*freq*(hjd-2458000)
    design = np.column_stack([np.ones(len(phase)),np.cos(phase),np.sin(phase)])
    weighted = design/error[:,None]
    beta = np.linalg.lstsq(weighted,residual/error,rcond=None)[0]
    chi = (residual-design@beta)/error
    cov = np.linalg.inv(weighted.T@weighted)*(chi@chi)/(len(chi)-3)
    amp = np.linalg.norm(beta[1:])
    err = np.sqrt(beta[1:]@cov[1:,1:]@beta[1:])/amp
    print(name, 'PS1 N, amplitude, conditional error [mmag]',len(data),amp*1000,err*1000)
dr5 = fits.getheader(here / 'lamost-257903093-dr5.fits.gz')
print('Same observation in DR5',dr5['SUBCLASS'],299792.458*dr5['Z'],299792.458*dr5['Z_ERR'])
