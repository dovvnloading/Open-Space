"""Rebuild the compact measurement tables from the public archives."""
import argparse
import csv
import hashlib
import gzip
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parent
SOURCE = '363580870076290816'
CONTROLS = ['684205069', '684203229', '684203029']
TARGET = '684208204'
LOW_RES = ['392902126', '491706223', '1189315023']


def write_table(path, columns, rows):
    with path.open('w', newline='') as out:
        writer = csv.writer(out, lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(rows)


def build(cache, out):
    cache.mkdir(exist_ok=True, parents=True)
    out.mkdir(exist_ok=True, parents=True)
    receipts = []

    def read(name, url):
        path = cache / name
        if not path.exists():
            with urlopen(url, timeout=180) as response:
                path.write_bytes(response.read())
        data = path.read_bytes()
        receipts.append({'name': name, 'url': url, 'sha256': hashlib.sha256(data).hexdigest()})
        return data

    def lam_sql(sql):
        return 'https://www.lamost.org/dr12/v2.0/sql/q?' + urlencode({'sql': sql, 'output.fmt': 'csv'})

    candidates = {}
    cat_cols = ['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'rv_amplitude_robust',
                'fit_period', 'fit_period_errlow', 'fit_period_errup', 'fit_companion_mass',
                'fit_companion_mass_errlow', 'fit_companion_mass_errup', 'mass_significance', 'ruwe']
    for kind in ('MS', 'RGB'):
        url = 'https://zenodo.org/api/records/19181131/files/' + kind + '+BH_candidates_from_GaiaDR3_summary_diagnostics.csv/content'
        for row in csv.DictReader(io.StringIO(read(kind + '.csv', url).decode())):
            if row['flag_quality'] == 'True' and float(row['fit_companion_mass']) > 3.0:
                candidates[row['source_id']] = (kind, row)
    matches = {}
    for name, table, limit in [('kovalev', 'J/MNRAS/527/521/tableb1', 20000), ('guo', 'J/ApJS/278/46/table1', 1000)]:
        url = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?' + urlencode({'-source': table, '-out.all': '', '-out.max': limit})
        lines = [line for line in read(name + '.tsv', url).decode().splitlines() if line and not line.startswith('#')]
        columns = lines[0].split('\t')
        ids = {r.split('\t')[columns.index('GaiaDR3')].strip() for r in lines[3:]}
        matches[name] = ids
    write_table(out / 'screen.csv', ['sample', *cat_cols, 'kovalev_sb2', 'guo_orbit'],
                [[kind, *[row[k] for k in cat_cols], sid in matches['kovalev'], sid in matches['guo']]
                 for sid, (kind, row) in sorted(candidates.items())])

    # Same observing plate, rank bright comparison spectra with template-like stellar parameters by S/N.
    sql = ("select * from med_stellar where planid='NT010031N374729M01' and coadd=1 and band='B' "
           "and teff_lasp between 5000 and 6200 and logg_lasp>3.7 and snr>30 order by snr desc limit 10")
    controls = list(csv.DictReader(io.StringIO(read('bluecontrols.csv', lam_sql(sql)).decode('utf-8-sig'))))[:3]
    if [r['obsid'] for r in controls] != CONTROLS:
        raise ValueError('Control selection changed in archive')
    fields = ['obsid', 'gaia_source_id', 'snr', 'rv_b1', 'teff_lasp', 'logg_lasp']
    write_table(out / 'controls.csv', fields, [[r[k] for k in fields] for r in controls])
    spectra = {}
    for obs in [TARGET, *CONTROLS, *LOW_RES]:
        service = 'spectrum' if obs in LOW_RES else 'medspectrum'
        raw = read(obs + '.fits.gz', 'https://www.lamost.org/dr12/v2.0/' + service + '/fits/' + obs)
        with fits.open(io.BytesIO(gzip.decompress(raw))) as hdus:
            extensions = ['COADD'] if obs in LOW_RES else ['COADD_B']
            if obs == TARGET:
                extensions += ['B-84126121', 'B-84126145']
            for name in extensions:
                row = hdus[name].data[0]
                w = row['WAVELENGTH'].astype(float)
                flux = row['FLUX'].astype(float)
                mask = 'ANDMASK' if 'ANDMASK' in hdus[name].columns.names else 'PIXMASK'
                good = (row['IVAR'] > 0) & (row[mask] == 0) & np.isfinite(w + flux)
                continuum = median_filter(np.interp(w, w[good], flux[good]), size=99)
                absorption = np.where(good, np.clip(1 - flux / continuum, -.6, .9), 0)
                keep = (w > 3850) & (w < 5850) if obs in LOW_RES else (w > 4870) & (w < 5340)
                spectra[obs + '_' + name] = np.stack((w[keep], absorption[keep])).astype('<f8')
    np.savez_compressed(out / 'spectra.npz', **spectra)

    lights = []
    for sector in (17, 57, 84):
        filename = f'hlsp_qlp_tess_ffi_s{sector:04}-0000000259010246_tess_v01_llc.fits'
        uri = f'mast:HLSP/qlp/s{sector:04}/0000/0002/5901/0246/{filename}'
        raw = read(filename, 'https://mast.stsci.edu/api/v0.1/Download/file?uri=' + uri)
        with fits.open(io.BytesIO(raw)) as hdus:
            data = hdus[1].data
            col = 'KSPSAP_FLUX' if sector == 17 else 'DET_FLUX'
            t, v = data['TIME'], data[col]
            mask = (data['QUALITY'] == 0) & np.isfinite(t + v)
            t, v = t[mask], v[mask]
            mask = abs(v - np.median(v)) < .035
            t, v = t[mask], v[mask] / np.median(v[mask])
            bins = np.round(24 * t).astype(int)
            for hour in np.unique(bins):
                use = bins == hour
                lights.append((sector, f'{np.mean(t[use]):.9f}', f'{np.mean(v[use]):.9f}', int(sum(use))))
    write_table(out / 'tess.csv', ['sector', 'btjd', 'flux', 'n'], lights)
    rows = []
    for obj, filename, x, y in [('target', 'msztf.csv', '15.773', '36.4274'),
                                ('neighbor', 'nearztf.csv', '15.7779733874', '36.4253118749')]:
        url = ('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20'
               + x + '%20' + y + '%200.0007&FORMAT=csv&BAD_CATFLAGS_MASK=32768')
        for row in csv.DictReader(io.StringIO(read(filename, url).decode())):
            rows.append([obj, *[row[k] for k in ('oid', 'expid', 'hjd', 'mag', 'magerr', 'catflags', 'filtercode', 'ra', 'dec')]])
    write_table(out / 'ztf.csv', ['position', 'oid', 'expid', 'hjd', 'mag', 'magerr', 'catflags', 'filtercode', 'ra', 'dec'], rows)
    (out / 'sources.json').write_text(json.dumps(receipts, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True, help='Folder for original downloads')
    parser.add_argument('--out', type=Path, required=True, help='Folder for regenerated evidence')
    args = parser.parse_args()
    build(args.cache, args.out)
