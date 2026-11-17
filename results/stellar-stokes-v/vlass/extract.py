"""Check all V footprints at DATE-OBS, then archive corresponding V and I cutouts."""
import argparse
import csv
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
DOWNLOAD = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VLASS/'
SAMPLES = ('hosts', 'gaia_nearby', 'cns5', 'dwarfs')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_bytes(name, lo, hi):
    url = DOWNLOAD + quote(name, safe='')
    headers = {'Range': f'bytes={lo}-{hi}'}
    for attempt in range(4):
        try:
            with urlopen(Request(url, headers=headers), timeout=100) as res:
                if res.status != 206 or not res.headers['Content-Range'].startswith(f'bytes {lo}-{hi}/'):
                    raise ValueError('Response range mismatch')
                data = res.read()
                if len(data) != hi - lo + 1:
                    raise ValueError('Truncated payload')
                return data
        except (OSError, ValueError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def load_header(name):
    data = read_bytes(name, 0, 2880 * 128 - 1)
    end = next(i + 80 for i in range(0, len(data), 80) if data[i:i + 8] == b'END     ')
    size = ((end + 2879) // 2880) * 2880
    hdr = fits.Header.fromstring(data[:end].decode('ascii'))
    if hdr['BUNIT'].strip().lower() != 'jy/beam' or hdr['BITPIX'] not in (-32,-64):
        raise ValueError(f'Unsupported FITS units/encoding in {name}')
    if any(hdr[f'NAXIS{i}'] != 1 for i in range(3, hdr['NAXIS'] + 1)):
        raise ValueError('Ambiguous nondegenerate higher axis')
    return hdr, size, sha(data[:size])


def footprint(hdr, sample):
    wcs = WCS(hdr).celestial
    epoch = Time(hdr['DATE-OBS']).jyear
    # Exact great-circle tangent proper motion (radial velocity and parallax omitted).
    cat = SkyCoord(sample.ra_deg.to_numpy(), sample.dec_deg.to_numpy(), unit='deg')
    pmra, pmdec = sample.pmra_masyr, sample.pmdec_masyr
    movable = np.isfinite(sample[['ref_epoch','pmra_masyr','pmdec_masyr']]).all(axis=1)
    pa = np.arctan2(pmra.fillna(0), pmdec.fillna(0)).to_numpy()
    dist = np.where(movable, (epoch - sample.ref_epoch) * np.hypot(pmra,pmdec), 0)
    at = cat.directional_offset_by(pa*u.rad, dist*u.mas)
    near = at.separation(SkyCoord(*wcs.wcs.crval, unit='deg')).deg < 1.1
    rows = []
    for idx in np.flatnonzero(near):
        x, y = wcs.world_to_pixel(at[idx]); nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
        if -0.5 <= x < nx-0.5 and -0.5 <= y < ny-0.5:
            s = sample.iloc[idx]
            rows.append(dict(sample=s['sample'], target=s.hostname, gaia_dr3_id=s.gaia_dr3_id,
                             x=float(x), y=float(y), ra=at[idx].ra.deg, dec=at[idx].dec.deg,
                             astrometry=s.astrometry, dist_pc=s.dist_pc, pm_applied=bool(movable.iloc[idx])))
    return rows


def cutout(name, hinfo, sky, outfile):
    hdr, header_n, header_sha = hinfo
    wcs = WCS(hdr).celestial
    x, y = np.rint(wcs.world_to_pixel_values(*sky)).astype(int)
    step = np.abs(np.linalg.det(wcs.pixel_scale_matrix))**0.5 * 3600
    side = int(np.ceil(48 / step))
    nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
    x0, y0, x1, y1 = max(0,x-side), max(0,y-side), min(nx,x+side+1), min(ny,y+side+1)
    if x0 >= x1 or y0 >= y1:
        raise ValueError('No overlap with corresponding image')
    bpp = abs(hdr['BITPIX']) // 8
    lo, hi = header_n + y0*nx*bpp, header_n + y1*nx*bpp - 1
    data = read_bytes(name, lo, hi)
    arr = np.frombuffer(data, dtype=f'>f{bpp}').reshape(y1-y0,nx)[:,x0:x1].copy()
    h = wcs.to_header()
    h['CRPIX1'] -= x0; h['CRPIX2'] -= y0
    for key in ('BUNIT','BMAJ','BMIN','BPA','DATE-OBS','RESTFRQ','OBJECT'):
        if key in hdr:
            h[key] = hdr[key]
    fits.PrimaryHDU(arr, header=h).writeto(ROOT/outfile, overwrite=True)
    return dict(path=outfile, product=name, header_bytes=header_n, header_sha256=header_sha,
                row_byte_range=[int(lo),int(hi)], rows_sha256=sha(data), data_sha256=sha(arr.tobytes()),
                original_xy0=[int(x0),int(y0)], shape=list(arr.shape))


def extract_field(info, samples):
    name = info['uri'].split('/')[-1]
    vheader = load_header(name)
    hdr, header_n, header_sha = vheader
    if not name.endswith('V_mean_minmaxrej.subim.fits'):
        raise ValueError(f'Unrecognised V product {name}')
    field = name.split('.')[4]
    if hdr['OBJECT'].strip() not in ('',field):
        raise ValueError('FITS field and filename disagree')
    grouped = {}
    for row in footprint(hdr, samples):
        # Keep exact aliases from separate catalogues; adjacent physical systems can share a beam.
        group = (round(row['x']), round(row['y']))
        grouped.setdefault(group, []).append(row)
    field_receipt = dict(v_product=name,field=field, obs_date=hdr['DATE-OBS'],
                         ref_hz=hdr['RESTFRQ'],shape=[hdr['NAXIS2'],hdr['NAXIS1']],
                         beam_maj_arcsec=3600*hdr['BMAJ'],beam_min_arcsec=3600*hdr['BMIN'],
                         header_bytes=header_n, header_sha256=header_sha, groups=[])
    if grouped:
        i_name = name.replace('.cc.','.se.').replace('.V_mean_minmaxrej.subim.fits','.v1.I.iter3.image.pbcor.tt0.subim.fits')
        iheader = load_header(i_name)
        if abs((Time(hdr['DATE-OBS'])-Time(iheader[0]['DATE-OBS'])).to_value('sec')) > 1:
            raise ValueError('Different I and V observation dates')
        for (x,y), rows in sorted(grouped.items()):
            tag = f'{field}_{x}_{y}'
            sky = rows[0]['ra'], rows[0]['dec']
            field_receipt['groups'].append(dict(rows=rows,
                V=cutout(name, vheader, sky, f'cutouts/{tag}_V.fits.gz'),
                I=cutout(i_name, iheader, sky, f'cutouts/{tag}_I.fits.gz')))
    return field_receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samples', type=Path, required=True)
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    (ROOT/'cutouts').mkdir(exist_ok=True)
    parts, frozen = [], {}
    for label in SAMPLES:
        path = args.samples / (label+'.csv')
        if not path.exists(): path = path.with_suffix('.csv.gz')
        t = pd.read_csv(path, dtype={'gaia_dr3_id':str}).fillna({'gaia_dr3_id':'','astrometry':''})
        numeric = ['ra_deg','dec_deg','ref_epoch','pmra_masyr','pmdec_masyr','dist_pc']
        t[numeric] = t[numeric].apply(pd.to_numeric, errors='coerce')
        parts.append(t[np.isfinite(t[['ra_deg','dec_deg']]).all(axis=1)].assign(sample=label))
        frozen[label] = dict(n=len(t), n_sky=len(parts[-1]), n_moving=t[numeric[:-1]].notna().all(axis=1).sum().item())
    samples = pd.concat(parts, ignore_index=True)
    fields = list(csv.DictReader((ROOT/'products.csv').open()))
    with ThreadPoolExecutor(args.workers) as ex:
        footprints = list(ex.map(lambda info: extract_field(info,samples),fields))
    (ROOT/'extraction.json').write_text(json.dumps(dict(samples=frozen, fields=footprints), indent=2)+'\n')


if __name__ == '__main__':
    main()
