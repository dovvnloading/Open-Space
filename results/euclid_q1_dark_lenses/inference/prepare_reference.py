"""Build the solar AB calibration and transcribe the spectroscopic control tables.

Usage: python prepare_reference.py INPUT_CACHE OUTPUT_DIR
Requires NumPy, lxml and Python 3. Downloads pinned public sources if absent.
"""
import csv
import hashlib
import sys
import urllib.request
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from lxml import html

BANDS = {'vis': 'VIS.vis', 'y': 'NISP.Y', 'j': 'NISP.J', 'h': 'NISP.H'}
SOLAR_URL = 'https://ssb.stsci.edu/trds/calspec/sun_reference_stis_002.fits'
F_URL = 'https://zenodo.org/api/records/19260200/files/lens_table.csv/content'
ROJAS_URL = 'https://arxiv.org/html/2503.15325v1'
FILES = {
    'solar.fits': SOLAR_URL, 'rojas.html': ROJAS_URL, 'lens_table.csv': F_URL,
    **{b + '.xml': 'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=Euclid/' + f for b, f in BANDS.items()}
}
HASHES = {
    'solar.fits': 'a50b70a5c6515322cee14a3fe8060d9e4d7260912ae8d3e85b7fa88302130ab3',
    'rojas.html': '010e5d49f51c418bfdb17d79f3cf5dab89fc0c9423a0db1007d0cab6c3c946c8',
    'lens_table.csv': '108016603f04cc00409277c71292d3b152cdceafa81705ef7551a0ff5c5fa619',
    'vis.xml': 'c4f9aedd15455adfea3c4598831f6c717f03382fab9c5dba82b5d23b90220468',
    'y.xml': '3fe7828235552cf4d7c1013988881ca8b702fe482f4671b4a48cea8645dac88b',
    'j.xml': 'e7420cd1e5529973a91b87172acd8fac7dc3f4be557d39dcc8f28df2405c3b30',
    'h.xml': 'cff7be57e9021eabddc5a68e50e11788cc3fc738cf6b8d190e4700a39441a623'
}


def prepare(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    for filename, url in FILES.items():
        f = cache / filename
        if not f.exists():
            f.write_bytes(urllib.request.urlopen(url, timeout=90).read())
        if hashlib.sha256(f.read_bytes()).hexdigest() != HASHES[filename]:
            raise ValueError('Unrecognised input bytes: ' + filename)

    # Pinned FITS binary table: 1467 rows of 5 big-endian float32 after 2 header blocks.
    data = np.frombuffer((cache / 'solar.fits').read_bytes(), dtype='>f4',
                         offset=5760, count=1467 * 5).reshape(-1, 5).astype(float)
    wave, flux = data[:, :2].T   # Angstrom, erg s^-1 cm^-2 Angstrom^-1 at 1 au
    z_grid = np.arange(0, 401) * 0.01
    solar = {}
    for band in BANDS:
        tree = ET.parse(cache / (band + '.xml'))
        entries = [list(map(float, [c.text for c in row])) for row in tree.findall('.//{*}TR')]
        w_obs, response = np.array(entries).T
        mags = []
        for z in z_grid:
            w = w_obs / (1 + z)
            if w.min() < wave.min() or w.max() > wave.max():
                mags.append(float('nan'))
                continue
            # Integrate at the union of solar and filter sampling wavelengths.
            x = np.unique(np.concatenate((w, wave[(wave > w.min()) & (wave < w.max())])))
            throughput = np.interp(x, w, response)
            flam = np.interp(x, wave, flux) * (149597870700 / 3.085677581491367e17) ** 2
            fnu = np.trapezoid(flam * x * throughput, x) / (2.99792458e18 * np.trapezoid(throughput / x, x))
            mags.append(-2.5 * np.log10(fnu / 3.631e-20))
        solar[band] = mags
    with (out / 'sun_shifted_band_ab.csv').open('w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['z', *BANDS])
        for idx, z in enumerate(z_grid):
            w.writerow([f'{z:.2f}', *[f'{solar[b][idx]:.5f}' for b in BANDS]])

    h = html.fromstring((cache / 'rojas.html').read_bytes())
    tables = {t.getparent().get('id'): t for t in h.xpath('//table[contains(@class,"ltx_tabular")]')}
    redshifts = {}
    preceding = None
    for tr in tables['S4.T3'].xpath('.//tr')[1:]:
        cells = [' '.join(c.text_content().split()) for c in tr.xpath('./td')]
        name = cells[0] or preceding
        preceding = name
        redshifts.setdefault(name, {})
        if cells[3]:
            redshifts[name]['zl'] = cells[3]
        if cells[4]:
            redshifts[name].update(zs=cells[4], quality='Palomar_' + cells[5].split('/')[-1])
    redshifts['EUCL J174613.92+662840.2'] = {'zs': '1.303', 'quality': 'DESI_likely_single_OII'}
    with (out / 'rojas_spectroscopy.csv').open('w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['name', 'ra_deg', 'dec_deg', 'zl_spec', 'zs_spec', 'zs_quality'])
        for table_id in ('S3.T1', 'S3.T2'):
            for tr in tables[table_id].xpath('.//tr')[1:]:
                cells = [' '.join(c.text_content().split()) for c in tr.xpath('./td')]
                name = cells[0]
                if 'J17' not in name and 'J18' not in name:
                    continue   # exclude COSMOS pre-Q1 fields
                spec = redshifts.get(name, {})
                w.writerow([name, cells[1], cells[2], spec.get('zl', cells[3]),
                            spec.get('zs', ''), spec.get('quality', '')])

    with (out / 'ecker_bright_q1.csv').open('w') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['name','ra_deg','dec_deg','zl_phot','zl_spec','VI','mag_IE','sersic_n','log_Mstar','theta_E'])
        for row in (cache / 'lens_table.csv').read_text().splitlines()[1:]:
            # Some published rows quote the entire CSV line; only its first ten simple fields are used.
            cells = row.lstrip('"').split(',', 10)[:10]
            w.writerow(['' if c.strip() in ('...', '–', '') else c.strip() for c in cells])


if __name__ == '__main__':
    prepare(Path(sys.argv[1]), Path(sys.argv[2]))
