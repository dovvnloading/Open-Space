"""Reproduce the Swift J1745-26 measurements. The default run uses the included tables."""
import argparse
import csv
import hashlib
import math
import statistics
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TARGET = '13284310018508'
RADIO_RA = 266.295204166667
RADIO_DEC = -26.4035


def read_table(name):
    with (ROOT / name).open() as stream:
        return list(csv.DictReader(stream))


def weighted(rows, key, error):
    weights = [1 / float(row[error]) ** 2 for row in rows]
    return (sum(float(row[key]) * w for row, w in zip(rows, weights)) / sum(weights),
            1 / math.sqrt(sum(weights)))


def angular_offset(ra, dec):
    return 3600 * math.hypot((ra - RADIO_RA) * math.cos(math.radians(RADIO_DEC)),
                             dec - RADIO_DEC)


def catalog():
    photometry = read_table('photometry.csv')
    frame = {r['sourceid']: r for r in read_table('field.csv')}
    clean = [r for r in photometry if r['phot_flag'] == '0' and r['ambiguous_match'] == '0']
    target = [r for r in clean if r['sourceid'] == TARGET]
    grouped = defaultdict(list)
    for row in clean:
        grouped[row['sourceid'], int(float(row['mjdobs']))].append(row)
    early = grouped[TARGET, 56183]
    magnitude, error = weighted(early, 'mag', 'emag')
    print(f'Pretrigger Ks {magnitude:.3f} ± {error:.3f} mag (formal errors)')
    print('Centroid offsets (arcsec)', [round(angular_offset(float(r['ra']),float(r['de'])), 3) for r in early])
    mean_epoch = statistics.mean(float(r['mjdobs']) for r in early)
    print(f'Lead over BAT trigger {56186.38618 - mean_epoch:.4f} days')
    changes = {}
    for source, reference in frame.items():
        if not reference['phot_ks_mean_mag'] or not 14 <= float(reference['phot_ks_mean_mag']) <= 18:
            continue
        before, after = grouped[source, 56154], grouped[source, 56183]
        if before and after:
            changes[source] = statistics.median(float(r['mag']) for r in before) - statistics.median(float(r['mag']) for r in after)
    amplitude = changes.pop(TARGET)
    print(f'Initial brightening {amplitude:.3f} mag; {len(changes)} local controls, greatest absolute change {max(map(abs, changes.values())):.3f} mag')
    for day in (56497, 56511, 56570):
        value, sigma = weighted(grouped[TARGET, day], 'mag', 'emag')
        print(f'MJD {day}: Ks {value:.3f} ± {sigma:.3f}')
    bat = read_table('bat.csv')
    for start, stop in [(56183,56185), (56186,56186), (56497,56506), (56570,56583)]:
        used = [r for r in bat if start <= int(r['mjd']) <= stop and r['data_flag'] == '0']
        flux, err = weighted(used,'rate','error')
        print(f'BAT {start}–{stop}: {flux:.6f} ± {err:.6f} counts/cm²/s, {len(used)} days')
    assert len(target) == 44 and len(early) == 2 and len(changes) == 16
    assert round(magnitude,3) == 15.772 and round(amplitude,3) == 1.493


def pixel_replay(directory):
    import warnings
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS, FITSFixedWarning

    directory.mkdir(parents=True, exist_ok=True)
    stars = read_table('field.csv')
    coords = np.array([[float(r['ra']), float(r['de'])] for r in stars])
    target = next(i for i, r in enumerate(stars) if r['sourceid'] == TARGET)
    coords[target] = RADIO_RA, RADIO_DEC
    separation = 3600 * np.hypot((coords[:, None, 0] - coords[:, 0]) * np.cos(np.radians(RADIO_DEC)),
                                 coords[:, None, 1] - coords[:, 1])
    np.fill_diagonal(separation, np.inf)
    standards = [i for i, r in enumerate(stars) if r['ks_p2'] and r['ks_p98']
                 and 12.8 < float(r['phot_ks_mean_mag']) < 15.7
                 and float(r['ks_p98']) - float(r['ks_p2']) < 0.30 and separation[i].min() > 1.15]
    assert len(standards) == 10
    # Some archive headers use the historical spelling of WCS keywords.
    warnings.simplefilter('ignore', FITSFixedWarning)
    for row in read_table('images.csv'):
        path = directory / (row['catid'] + '.fits')
        if not path.exists():
            with urllib.request.urlopen(row['url'], timeout=120) as request:
                data = request.read()
            if not data.startswith(b'SIMPLE  ='):
                raise ValueError(f'ESO returned an invalid FITS response for {path.name}')
            path.write_bytes(data)
        original = hashlib.sha256(path.read_bytes()).hexdigest() == row['original_cutout_sha256']
        with fits.open(path) as hdus:
            image = hdus[1].data.astype(float)
            positions = WCS(hdus[1].header).all_world2pix(coords, 0)
        yy, xx = np.indices(image.shape)
        offsets = []
        for i in standards:
            x, y = positions[i]
            distance = np.hypot(xx-x, yy-y)
            sky = np.median(image[(distance >= 3) & (distance <= 5)])
            weight = np.maximum(image - sky, 0) * (distance < 2.5)
            offsets.append([(xx*weight).sum()/weight.sum()-x, (yy*weight).sum()/weight.sum()-y])
        positions += np.median(offsets, axis=0)
        flux = {}
        for i in [target, *standards]:
            x, y = positions[i]
            distance = np.hypot(xx-x, yy-y)
            sky = np.median(image[(distance > 3.5) & (distance < 5)])
            fraction = np.zeros(image.shape)
            for dx in np.linspace(-0.45, 0.45, 10):
                for dy in np.linspace(-0.45, 0.45, 10):
                    fraction += ((xx-x+dx)**2 + (yy-y+dy)**2 < 1.7**2) / 100
            flux[i] = np.sum((image-sky) * fraction)
        zp = np.median([float(stars[i]['phot_ks_mean_mag']) + 2.5*np.log10(flux[i]) for i in standards])
        mag = zp - 2.5*np.log10(flux[target])
        # SODA can regenerate non-science FITS header cards on later requests.
        print(f'{path.name}: fixed-aperture Ks = {mag:.3f}, original bytes = {original}')
        if abs(mag - float(row['aperture_ks'])) > 0.003:
            raise ValueError(f'Photometry differs from the frozen result for {path.name}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--images', type=Path, help='cache and independently measure the 12 ESO FITS cutouts')
    args = parser.parse_args()
    catalog()
    if args.images:
        pixel_replay(args.images)
