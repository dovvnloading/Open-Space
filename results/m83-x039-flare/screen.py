"""Rebuild the 12-observation time-clustering discovery screen (~0.3 GB)."""
import argparse
import csv
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from scipy.spatial import cKDTree
from scipy.special import betainc

from analyze import sky_wcs

OBSIDS = (793, 2064, 12992, 12993, 12994, 12995, 12996, 13202, 13241, 13248, 14332, 14342)
RUN_SIZES = (6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256)
CATALOG = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?' + urlencode({
    '-source': 'J/ApJS/212/21/xray', '-out': 'ID,RAJ2000,DEJ2000', '-out.max': '500',
})


def retrieve(path, url):
    if not path.exists():
        with urlopen(url, timeout=180) as response:
            content = response.read()
        temporary = path.with_suffix('.part')
        temporary.write_bytes(content)
        temporary.replace(path)
    return path


def source_catalog(cache):
    lines = retrieve(cache / 'positions.tsv', CATALOG).read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith('ID\t'))
    rows = [line.split('\t') for line in lines[header+1:] if re.match(r'^X\d{3}\t', line)]
    if len(rows) != 458:
        raise ValueError(f'Expected 458 positions, received {len(rows)}')
    identifiers, ra, dec = zip(*rows)
    positions = SkyCoord(list(ra), list(dec), unit=(u.hourangle, u.deg))
    return identifiers, np.column_stack([positions.ra.deg, positions.dec.deg])


def event_file(cache, obsid):
    path = cache / f'{obsid}.fits.gz'
    if path.exists():
        return path
    url = f'https://cxc.cfa.harvard.edu/cdaftp/byobsid/{obsid % 10}/{obsid}/primary/'
    filenames = re.findall(r'href="([^"/]+_evt2.fits.gz)"', urlopen(url, timeout=45).read().decode())
    if len(filenames) != 1:
        raise ValueError(f'Ambiguous event file for {obsid}')
    return retrieve(path, url + filenames[0])


def scan(path, obsid, identifiers, sky_coordinates):
    with fits.open(path) as hdus:
        events = hdus['EVENTS']
        header = events.header
        if int(header['OBS_ID']) != obsid or 'M83' not in header['OBJECT'].replace(' ', '') and '5236' not in header['OBJECT']:
            raise ValueError(f'Input field or ObsID mismatch for {obsid}')
        positions = sky_wcs(header).all_world2pix(sky_coordinates, 1)
        off_axis = np.linalg.norm(positions - [header['TCRPX11'],header['TCRPX12']], axis=1) * .492 / 60
        radius = np.maximum(1.75, 1.3 + .07 * off_axis**2) / .492
        nearest = cKDTree(positions).query(positions, k=2)[0][:, 1]
        radius = np.minimum(radius, .45*nearest)
        table = events.data
        mask = (table['energy'] >= 350) & (table['energy'] <= 8000)
        coords = np.column_stack([table['x'][mask], table['y'][mask]])
        apertures = cKDTree(coords).query_ball_point(positions, radius)
        time = table['time'][mask]
        chips = table['ccd_id'][mask]
        gtis = {hdu.header['CCD_ID']: hdu.data for hdu in hdus if hdu.name == 'GTI'}
        rows = []
        for source_id, indices in zip(identifiers, apertures):
            n = len(indices)
            if n < min(RUN_SIZES):
                continue
            t = np.sort(time[indices])
            chip = np.bincount(chips[indices]).argmax()
            gti = gtis[chip]
            duration = gti['STOP'][-1] - gti['START'][0]
            best = (1., 0, 0.)
            for m in RUN_SIZES:
                if m > n:
                    continue
                dt = float(np.min(t[m-1:] - t[:n-m+1]))
                fraction = min(1., (dt + 3.241) / duration)
                bound = min(1., float(betainc(m-1, n-m+2, fraction)) * (n-m+1))
                best = min(best, (bound, m, dt))
            rows.append((source_id, obsid, n, *best))
        return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('screen-cache'))
    parser.add_argument('--output', type=Path, default=Path('screen.csv'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    identifiers, positions = source_catalog(args.cache)
    rows = []
    for obsid in OBSIDS:
        rows.extend(scan(event_file(args.cache, obsid), obsid, identifiers, positions))
    with args.output.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_id','obsid','photons','within_run_length_p_bound','run_length','span_s','family_p_bound'])
        for source, obsid, n, bound, m, dt in sorted(rows, key=lambda r: r[3]):
            writer.writerow((source, obsid, n, f'{bound:.9g}', m, f'{dt:.3f}', min(1, bound*len(rows)*len(RUN_SIZES))))
    print(f'{len(rows)} eligible source-observations')


if __name__ == '__main__':
    main()
