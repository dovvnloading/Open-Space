"""Check catalog triggers on the individual exposures contributing to each sky position."""
from scipy import ndimage as nd, optimize as opt
import pandas as pd
import argparse
from astropy.wcs.utils import proj_plane_pixel_scales
from pathlib import Path
from urllib.request import urlopen, Request
from io import BytesIO
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import json
import sys
import time
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pixel_check import byte_range, stamp, XX, YY, STEP, aperture
from fetch_catalogs import product_url


def get_prod(prod, RA, DEC):
    url = product_url(prod, 'i2d.fits')

    def read(start, size):
        return byte_range(url, start, size)[0]
    offset = 0
    hinfo = {}
    q = {}
    while True:
        hdrbuf = read(offset, 57600)
        while not any(hdrbuf[i:i + 8] ==
                      b'END     ' for i in range(0, len(hdrbuf), 80)):
            hdrbuf += read(offset + len(hdrbuf), 57600)
        ff = BytesIO(hdrbuf)
        h = fits.Header.fromfile(ff)
        p = offset + ff.tell()
        name = h.get('EXTNAME', 'PRIMARY')
        dims = tuple(h.get(f'NAXIS{i + 1}') for i in range(h['NAXIS']))
        nbytes = abs(h['BITPIX']) // 8 * (int(np.prod(dims))
                                          if dims else 0) + h.get('PCOUNT', 0)
        offset = p + 2880 * ((nbytes + 2879) // 2880)
        if name == 'SCI':
            w = WCS(h, naxis=2)
            xx, yy = [int(np.rint(float(v)))
                      for v in w.world_to_pixel_values(RA, DEC)]
            hinfo = {'x': xx, 'y': yy, 'mjd': h['MJD-BEG']}
        if name in ('CON', 'WHT', 'ERR'):
            dtype = {32: '>i4', -32: '>f4'}[h['BITPIX']]
            for plane in range(dims[2] if len(dims) > 2 else 1):
                raw = read(p + 4 * ((plane * dims[1] + yy) * dims[0] + xx), 4)
                q[name + str(plane)] = np.frombuffer(raw,
                                                     dtype=dtype)[0].item()
        if name == 'HDRTAB':
            hb = BytesIO(hdrbuf[:ff.tell()] +
                         read(p, nbytes) + bytes((-nbytes) % 2880))
            tab = fits.BinTableHDU.readfrom(hb)
            keys = [
                n for n in tab.columns.names if 'NAME' in n or n in (
                    'DETECTOR',
                    'EXPOSURE',
                    'ACT_ID',
                    'PATT_NUM',
                    'VISIT',
                    'FILTER',
                    'OBSERVTN',
                    'VISITGRP',
                    'SEQ_ID',
                    'PROGRAM',
                    'MJD-BEG')]
            rows = []
            for i, row in enumerate(tab.data):
                if q['CON' + str(i // 32)] & (1 << (i % 32)):
                    rows.append({'row': i, **{k: str(row[k]) for k in keys}})
            hinfo['exposures'] = rows
            break
        if name == 'ASDF':
            break
    return dict(product=prod, **hinfo, **q)


def calc(arg):
    p, f, RA, DEC, stamps, case, gain = arg
    fn = f['FILENAME']
    stem = fn.split('_o')[0]
    visit = stem.split('_')[0]
    url = f"https://stpubdata.s3.amazonaws.com/jwst/public/{visit[:7]}/{visit}/{stem}_i2d.fits"
    hraw, rc = byte_range(url, 0, 57600)
    buf = BytesIO(hraw)
    h0 = fits.Header.fromfile(buf)
    h = fits.Header.fromfile(buf)
    im, receipt = stamp((url, h, buf.tell(), rc), RA, DEC)
    im = im.astype('float64')
    ref = np.load(stamps / f'{case}_{p}.npz')['pixels'].astype('float64')
    # Fit source-scene translations independently for each dither relative to
    # own epoch.
    g = np.isfinite(ref) & np.isfinite(im) & (np.hypot(XX, YY) < .55)
    av = np.nanmedian(im)
    work = np.where(np.isfinite(im), im, av)

    def shift(v): return nd.shift(
        work, (v[1], v[0]), order=1, prefilter=False, mode='nearest')

    def loss(v):
        a = (shift(v) - ref)[g]
        return np.mean(np.minimum((a - np.median(a))**2, .15))
    sol = opt.minimize(
        loss, [
            0., 0.], method='Powell', bounds=(
            (-8, 8), (-8, 8))).x
    registered = shift(sol)
    band = p.split('-')[-1].lower()
    row = dict(epoch=p, product=stem, filter=band, dither=f['PATT_NUM'], dx_arcsec=sol[0] * STEP, dy_arcsec=sol[1] * STEP,
               covered016=float(np.mean(np.isfinite(im)[np.hypot(XX, YY) < .16])), median_bg=av, center_valid=np.isfinite(im[51, 51]), mjd=f['MJD-BEG'], flux016_njy=aperture(registered, 0, 0, .16) * gain, flux028_njy=aperture(registered, 0, 0, .28) * gain,
               ra=RA, dec=DEC)
    return row, dict(url=url, receipts=receipt)


def inspect(arg):
    p, f, RA, DEC = arg
    fn = f['FILENAME']
    v = fn.split('_')[0]
    url = f'https://stpubdata.s3.amazonaws.com/jwst/public/{v[:7]}/{v}/{fn}'
    row = dict(product=fn, filter=f['FILTER'], dither=f['PATT_NUM'])
    offset = 0
    while True:
        b, rc = byte_range(url, offset, 57600)
        buf = BytesIO(b)
        h = fits.Header.fromfile(buf)
        start = offset + buf.tell()
        dims = [h[f'NAXIS{k + 1}'] for k in range(h['NAXIS'])]
        nbytes = (int(np.prod(dims)) if dims else 0) * \
            (abs(h['BITPIX']) // 8) + h.get('PCOUNT', 0)
        offset = start + 2880 * ((nbytes + 2879) // 2880)
        which = h.get('EXTNAME', '')
        width = dims[0] if dims else None
        if which == 'SCI':
            w = WCS(h, naxis=2)
            xf, yf = map(float, w.world_to_pixel_values(RA, DEC))
            scale = np.mean(proj_plane_pixel_scales(w)) * 3600
            row.update(x_native=xf,
                       y_native=yf,
                       scale=scale,
                       nearest_detector_edge_arcsec=scale * min(xf,
                                                                yf,
                                                                width - 1 - xf,
                                                                dims[1] - 1 - yf))
        if which in ('SCI', 'DQ'):
            y0 = max(0, int(yf) - 20)
            y1 = min(dims[1], int(yf) + 21)
            r, rec = byte_range(
                url, start + y0 * width * 4, (y1 - y0) * width * 4)
            arr = np.frombuffer(r, dtype='>f4' if which ==
                                'SCI' else '>i4').reshape(-1, width)
            yy, xx = np.indices(arr.shape)
            dist = scale * np.hypot(xx - xf, yy + y0 - yf)
            mask = dist <= .16
            if which == 'SCI':
                row['finite_fraction'] = float(np.isfinite(arr[mask]).mean())
                row['sci_median'] = float(np.nanmedian(arr[mask]))
                row['sci_max'] = float(np.nanmax(arr[mask]))
            else:
                vals = (
                    arr.astype('i8') *
                    h.get(
                        'BSCALE',
                        1) +
                    h.get(
                        'BZERO',
                        0)).astype('u4')[mask]
                row.update(
                    dq_do_not_use=float(
                        np.mean(
                            (vals & 1) != 0)),
                    dq_jump_det=float(
                        np.mean(
                            (vals & 4) != 0)),
                    dq_outlier=float(
                        np.mean(
                            (vals & 16) != 0)),
                    dq_values='|'.join(
                        f'{u}:{n}' for u,
                        n in zip(
                            *
                            np.unique(
                                vals,
                                return_counts=True))))
                return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "candidate",
        choices=[
            "ceers_june-002",
            "cosmos_medium-009"])
    parser.add_argument(
        "stamps",
        type=Path,
        help="Directory populated by pixel_check.py")
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--dq",
        action="store_true",
        help="Also download detector-level CRF science/DQ patches")
    a = parser.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    base = HERE.parent
    hits = pd.read_csv(base / "coincidences.csv").set_index("candidate")
    hit = hits.loc[a.candidate]
    bands = ["f140m",
             "f210m"] if a.candidate == "ceers_june-002" else ["f162m",
                                                               "f182m"]
    pair = pd.read_csv(base / "pairs.csv")
    gains = pd.read_csv(base / "screen_counts.csv")
    config = pair.merge(gains, on=["pair", "band"])
    selected = config[(config.pair == hit["pair"]) & config.band.isin(bands)]
    products = {p: g for row in selected.itertuples()
                for p, g in [(row.reference, row.flux_scale), (row.search, 1.)]}
    if len(products) != 4:
        raise ValueError("Expected two matched-band pairs")
    with ThreadPoolExecutor(4) as pool:
        cov = list(pool.map(lambda p: get_prod(p, hit.ra, hit.dec), products))
    (a.output / "coverage.json").write_text(json.dumps(cov, indent=2) + "\n")
    exposures = [(p["product"], f) for p in cov for f in p["exposures"]]
    tasks = [(p, f, hit.ra, hit.dec, a.stamps, a.candidate, products[p])
             for p, f in exposures]
    with ThreadPoolExecutor(4) as pool:
        res = list(pool.map(calc, tasks))
    pd.DataFrame([r for r, _ in res]).to_csv(
        a.output / "dithers.csv", index=False)
    (a.output /
     "receipts.json").write_text(json.dumps([r for _, r in res], indent=2) +
                                 "\n")
    if a.dq:
        with ThreadPoolExecutor(4) as pool:
            n = list(
                pool.map(
                    inspect, [
                        (p, f, hit.ra, hit.dec) for p, f in exposures]))
        pd.DataFrame(n).to_csv(a.output / "native.csv", index=False)


if __name__ == "__main__":
    main()
