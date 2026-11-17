"""Compile CK04 + Polletta SWIRE spectra to relative Fnu for the blind search."""
import argparse
import concurrent.futures
import csv
import hashlib
import io
import json
import tarfile
import time
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
CK = 'https://ssb.stsci.edu/trds/grid/ck04models/ckp00/ckp00_{}.fits'
SWIRE = 'https://www.iasf-milano.inaf.it/~polletta/templates/data/swire_library.tar.gz'
TEMPERATURES = (3500, 4000, 4500, 5000, 5750, 6500, 7500, 9000, 12000)
GRID = np.geomspace(0.20, 6.0, 1200)


def compact(wave, flam):
    """Both libraries use Angstroms and F_lambda; normalization cancels c."""
    usable = np.isfinite(wave) & np.isfinite(flam) & (flam > 0)
    order = np.argsort(wave[usable])
    wave = wave[usable][order] * 1e-4
    fnu = flam[usable][order] * wave**2
    if wave[0] > GRID[0] or wave[-1] < GRID[-1]:
        raise ValueError('Template does not span the stored wavelength grid')
    values = np.interp(GRID, wave, fnu)
    return values / np.interp(2.2, GRID, values)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, default=HERE / 'localcache/template-inputs')
    ap.add_argument('--out', type=Path, default=HERE / 'templates')
    args = ap.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    args.out.mkdir(exist_ok=True, parents=True)
    receipt_path = args.out / 'sources.json'
    expected = {row['url']: row['sha256'] for row in json.loads(receipt_path.read_text())['inputs']} if receipt_path.exists() else {}

    def get(url):
        destination = args.cache / url.rsplit('/', 1)[1]
        if not destination.exists():
            for attempt in range(3):
                try:
                    with urlopen(url, timeout=60) as response:
                        raw = response.read()
                    destination.write_bytes(raw)
                    break
                except OSError:
                    if attempt == 2:
                        raise
                    time.sleep(2**attempt)
        raw = destination.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if url in expected and digest != expected[url]:
            raise ValueError(f'Input changed: {url}')
        return raw, {'url': url, 'bytes': len(raw), 'sha256': digest}

    urls = [CK.format(t) for t in TEMPERATURES] + [SWIRE]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        inputs = list(pool.map(get, urls))
    columns = {}
    for temp, (raw, _) in zip(TEMPERATURES, inputs):
        with fits.open(io.BytesIO(raw), memmap=False) as hdus:
            tab = hdus[1].data
            for gravity in ('g25', 'g45'):
                if np.any(tab[gravity] > 0):
                    columns[f'CK04_{temp}_{gravity}'] = compact(tab['WAVELENGTH'], tab[gravity])
    with tarfile.open(fileobj=io.BytesIO(inputs[-1][0]), mode='r:gz') as archive:
        for item in sorted(archive.getmembers(), key=lambda entry: entry.name):
            if item.name.endswith('_template_norm.sed'):
                wave, flam = np.loadtxt(archive.extractfile(item), unpack=True)
                columns['SWIRE_' + item.name.split('_')[0]] = compact(wave, flam)
    out = args.out / 'relative_fnu.csv'
    with out.open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['wave_um', *columns])
        for row in zip(GRID, *columns.values()):
            writer.writerow([f'{value:.7g}' for value in row])
    receipt = {'inputs': [r for _, r in inputs], 'output': out.name,
               'output_sha256': hashlib.sha256(out.read_bytes()).hexdigest(),
               'references': ['https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas',
                              'https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html',
                              'https://doi.org/10.1086/518113']}
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
