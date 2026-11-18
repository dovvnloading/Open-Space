"""Independent ASPIICS continuum counterfactual and primary-output comparison."""
import argparse
import datetime as dt
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from astropy.io import fits


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def wideband_parents(header):
    return {str(header[key]) for key in header if key.startswith('PARENT')
            and str(header[key]).startswith('aspiics_wb_l2_')}


def read_input(raw, manifest, item, orbit):
    path = raw / item['name']
    if sha256(path) != manifest[path.name]['sha256']:
        raise ValueError(f'Checksum mismatch: {path.name}')
    with fits.open(path) as hdul:
        header = hdul[0].header.copy()
        image = np.array(hdul[0].data, dtype=np.float64)
    expected = {'CYCLE_ID': item['cycle'], 'QUALITY': item['quality'],
                'ORBIT_ID': int(orbit), 'FILENAME': item['name']}
    for key, value in expected.items():
        if header[key] != value:
            raise ValueError(f'{path.name}: {key} mismatch')
    if not path.name.endswith('_v03.fits'):
        raise ValueError(f'Unexpected version: {path.name}')
    return header, image


def analyze(raw, manifest, orbit, sample):
    inputs = [read_input(raw, manifest, sample[key], orbit)
              for key in ('line', 'donor', 'next_brightness')]
    (header, line), (donor_header, donor), (next_header, following) = inputs
    if len(wideband_parents(header)) != 3 or (
            wideband_parents(header) != wideband_parents(donor_header)):
        raise ValueError('Wideband parent set mismatch')
    for key in ('CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2'):
        if any(abs(float(h[key]) - float(header[key])) > 1e-8
               for h in (donor_header, next_header)):
            raise ValueError(f'Image geometry mismatch: {key}')
    history = (str(entry) for entry in header['HISTORY'])
    conversion = next(float(re.search(r'MSB equals\s+([\deE+.-]+)', entry)[1])
                      for entry in history if 'MSB equals' in entry)
    times = [dt.datetime.fromisoformat(h['DATE-AVG'])
             for h in (header, donor_header, next_header)]
    weight = (times[0] - times[1]).total_seconds() / (
        times[2] - times[1]).total_seconds()
    if not 0 < weight < 1:
        raise ValueError('Line time is not bracketed')
    y, x = np.ogrid[:line.shape[0], :line.shape[1]]
    x = (x + 1 - header['CRPIX1']) * header['CDELT1'] + header['CRVAL1']
    y = (y + 1 - header['CRPIX2']) * header['CDELT2'] + header['CRVAL2']
    radius = np.sqrt(x*x + y*y) / header['RSUN_OBS']
    continuum = (np.isfinite(donor) & (donor > 0) & np.isfinite(following)
                 & (following > 0))
    positive = continuum & np.isfinite(line) & (line > 0)
    with np.errstate(invalid='ignore'):
        change = weight * (donor - following)
        flips = positive & (line + conversion * change <= 0)
        robust = flips.copy()
        for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            neighbor = np.roll(following, (dy, dx), axis=(0, 1))
            support = np.roll(np.isfinite(following) & (following > 0),
                              (dy, dx), axis=(0, 1))
            # The annuli are far from image edges, so roll wrap is never used.
            shifted = 0.5 * (following + neighbor)
            robust &= support & (line + conversion * weight *
                                 (donor - shifted) <= 0)
    annuli = []
    for low, high in ((1.2, 1.4), (1.4, 1.7)):
        ring = (radius >= low) & (radius < high)
        fractions = np.abs(change[continuum & ring] / donor[continuum & ring])
        annuli.append({
            'geometric_pixels': int(ring.sum()),
            'continuum_valid_pixels': int((continuum & ring).sum()),
            'line_positive_pixels': int((positive & ring).sum()),
            'flips': int((flips & ring).sum()),
            'robust_flips': int((robust & ring).sum()),
            'abs_fractional_K_median': float(np.median(fractions)),
            'abs_fractional_K_p95': float(np.percentile(fractions, 95)),
        })
    return {'orbit': int(orbit), 'line': sample['line']['name'], 'annuli': annuli}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection', type=Path, required=True)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compare', type=Path)
    args = parser.parse_args()
    selection = json.loads(args.selection.read_text(encoding='utf-8'))
    manifest = json.loads(args.manifest.read_text(encoding='utf-8'))
    result = [analyze(args.raw, manifest, orbit, sample)
              for orbit in ('795', '772') for sample in selection[orbit]['sample']]
    differences = []
    if args.compare:
        primary = {row['line']: row for row in
                   json.loads(args.compare.read_text(encoding='utf-8'))}
        for row in result:
            for index, (actual, expected) in enumerate(zip(
                    row['annuli'], primary[row['line']]['annuli'], strict=True)):
                for key, value in actual.items():
                    tolerance = 1e-15 if key.startswith('abs_fractional') else 0
                    if abs(value - expected[key]) > tolerance:
                        differences.append([row['line'], index, key, value,
                                            expected[key]])
        if len(result) != len(primary):
            differences.append(['sample count', len(result), len(primary)])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'results': result, 'differences': differences},
                                      indent=2) + '\n', encoding='utf-8', newline='\n')
    print(f'{len(result)} triples, {sum(len(x["annuli"]) for x in result)} annuli; '
          f'{len(differences)} comparison differences')
    if differences:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
