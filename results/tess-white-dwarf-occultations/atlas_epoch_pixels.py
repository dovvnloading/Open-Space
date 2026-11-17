"""Extract independent ATLAS Sutherland image photometry at two WD positions."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from astropy.wcs import WCS, FITSFixedWarning

from ffi_acquire import propagated

HERE = Path(__file__).resolve().parent
TICS = ('220578510', '91529455')
START, STOP = '2022-01-01', '2025-03-01'
BASE = 'https://catch-api.astro.umd.edu/fixed?'
SITE = EarthLocation.from_geodetic(20.81, -32.38, 1790)
iers.conf.auto_download = False


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(path, url):
    if not path.exists():
        with urlopen(url,timeout=120) as response:
            path.write_bytes(response.read())


def measure(obs, star, cache):
    name = obs['product_id'].split(':')[-1]
    u = obs['cutout_url'].replace('0.16deg','240arcsec')
    if 'size=240arcsec' not in u:
        raise ValueError('Unexpected CATCH cutout size')
    path = cache/name
    result = dict(product_id=obs['product_id'], cutout_url=u, band=obs['filter'],
                  mjd_utc=.5*(obs['mjd_start']+obs['mjd_stop']))
    fetch(path,u)
    result['sha256'] = sha(path)
    with fits.open(path, memmap=False) as ff:
        h = ff[0].header.copy()
        d = np.asarray(ff[0].data,dtype=float)
    assert h['FILTER'].strip() == obs['filter'] and abs(result['mjd_utc']-(h['MJD-OBS']+h['EXPTIME']/172800.)) < .0001
    # The images supply both TPV and alternate SIP distortion coefficients.
    # Apply the declared TPV WCS once. Combining both double-counts distortion.
    if h['CTYPE1'] != 'RA---TPV' or h['CTYPE2'] != 'DEC--TPV':
        raise ValueError('Unexpected primary WCS')
    for key in list(h):
        if re.fullmatch(r'(A|B|AP|BP)_(ORDER|\d+_\d+)',key):
            del h[key]
    epoch = 2000.+(result['mjd_utc']-51544.5)/365.25
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        x, y = WCS(h,naxis=2).all_world2pix(*propagated(star,epoch),0)
    fwhm, zp = float(h['FWHM']), float(h['MAGZPT'])
    result.update(x=float(x), y=float(y), fwhm=fwhm, zp=zp, zp_rms=h['ZPRMS'])
    if not (1. <= fwhm <= 6. and 15. <= zp <= 28. and h['ZPNPT'] >= 50 and h['ZPRMS'] < .2):
        return dict(result,status='image calibration cut')
    yy, xx = np.indices(d.shape)
    radius = np.hypot(xx-x,yy-y)
    ap = radius <= max(3.,1.55*fwhm)
    ann = (radius >= max(9.,2.8*fwhm)) & (radius <= max(14.,4.*fwhm))
    if min(x,y,d.shape[1]-1-x,d.shape[0]-1-y) < max(14.,4.*fwhm):
        return dict(result,status='edge')
    if np.any(~np.isfinite(d[ap])) or np.any((d[ap] >= h['SATURATE']) | (d[ap] == 0)):
        return dict(result,status='saturated or invalid aperture pixel')
    sky = np.nanmedian(d[ann])
    signal = np.sum(d[ap]-sky)
    w = np.maximum(0.,d[ap]-sky)
    cent = float(np.hypot(np.sum(w*(xx[ap]-x)),np.sum(w*(yy[ap]-y)))/w.sum()) if w.sum()>0 else None
    # Preserve signed flux (including potential dropouts) with a common AB zero point.
    f25 = signal/h['EXPTIME']*10**(.4*(25-zp))
    return dict(result,status='measured', f25=float(f25),sky=float(sky),centroid_offset_pix=cent)


def main(cache,out,workers):
    cache.mkdir(exist_ok=True, parents=True)
    out.mkdir(exist_ok=True, parents=True)
    parent = list(csv.DictReader((HERE/'extension/extension_parent.csv').open()))
    checks = json.loads(gzip.decompress((HERE/'extension/extension_image_checks.json.gz').read_bytes()))
    summary, epochs = [], []
    rng = np.random.default_rng(149110)
    for tic in TICS:
        star = next(r for r in parent if r['tic_requested'] == tic)
        ra, dec = float(star['RA_ICRS']),float(star['DE_ICRS'])
        url = BASE+urlencode(dict(ra=round(ra,8),dec=round(dec,8),sources='atlas_sutherland',
                                   start_date=START,stop_date=STOP))
        path = cache/(tic+'.json')
        fetch(path,url)
        catch = json.loads(path.read_text())
        assert catch['count'] == len(catch['data']) and catch['query']['sources']==['atlas_sutherland']
        rows = sorted(catch['data'], key=lambda r:r['mjd_start'])[-200:]
        assert len({r['product_id'] for r in rows}) == len(rows)
        with ThreadPoolExecutor(workers) as pool:
            recs = list(pool.map(lambda r:measure(r,star,cache),rows))
        obs = [r for r in recs if r['status'] == 'measured']
        time = Time([r['mjd_utc'] for r in obs],format='mjd',scale='utc',location=SITE)
        btjd = (time.tdb+time.light_travel_time(SkyCoord(ra,dec,unit='deg'))).jd-2457000.
        for r,b in zip(obs,btjd):
            r['btjd_tdb'] = float(b)
        for band in sorted({r['band'] for r in obs}):
            group = [r for r in obs if r['band']==band]
            base = np.median([r['f25'] for r in group])
            for r in group:
                r['relative_flux'] = r['f25']/base
        y = np.array([r['relative_flux'] for r in obs])
        night = np.floor([r['mjd_utc']+.5 for r in obs]).astype(int)
        candidate = next(r for r in checks if r['tic']==tic)
        folds = []
        for trig in candidate['image_checks']:
            if 'period_days' not in trig:
                continue
            p, e, w = trig['period_days'], trig['epoch_btjd'], trig['duration_min']
            inside = abs((btjd-e+.5*p)%p-.5*p) < .5*w/1440.
            precision = np.array([r['zp_rms'] < .05 and r['zp'] > 20.5 for r in obs])
            valid_in, valid_out = inside & precision, ~inside & precision
            visits = np.unique(night)
            shifts = []
            if len(set(night[inside])) >= 2:
                for _ in range(1000):
                    drawn = rng.choice(visits, size=len(visits), replace=True)
                    idx = np.concatenate([np.flatnonzero(night == n) for n in drawn])
                    groups = inside[idx]
                    if groups.any() and (~groups).any():
                        shifts.append(np.median(y[idx][~groups])-np.median(y[idx][groups]))
            folds.append(dict(period_days=p,epoch_btjd=e,duration_min=w, n_in=int(inside.sum()),
                              n_out=int((~inside).sum()), nights_in=len(set(night[inside])),
                              stricter_calibration_in=int(valid_in.sum()),
                              stricter_calibration_out=int(valid_out.sum()),
                              stricter_depth_median=float(np.median(y[valid_out])-np.median(y[valid_in])) if valid_in.any() else None,
                              night_bootstrap_95=list(np.percentile(shifts,[2.5,97.5])) if shifts else None,
                              depth_median=float(np.median(y[~inside])-np.median(y[inside])) if inside.any() else None,
                              out_of_window_mad=float(1.4826*np.median(abs(y[~inside]-np.median(y[~inside]))))))
        summary.append(dict(tic=tic, gaia=star['GaiaEDR3'],catch_url=url, catch_sha256=sha(path),
                            catch_count=catch['count'], attempted=len(recs),measured=len(obs), first_mjd_utc=time.mjd.min(),
                            last_mjd_utc=time.mjd.max(), median_centroid_offset_pix=np.median([r['centroid_offset_pix'] for r in obs]),
                            folds=folds))
        epochs.extend([dict(tic=tic,**r) for r in recs])
        print(tic,len(obs),flush=True)
    data=''.join(json.dumps(r,separators=(',',':'),allow_nan=False)+'\n' for r in epochs)
    (out/'epochs.jsonl.gz').write_bytes(gzip.compress(data.encode(),mtime=0))
    (out/'folds.json').write_text(json.dumps(dict(epochs_sha256=sha(out/'epochs.jsonl.gz'),stars=summary), indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('--workers',type=int,default=6)
    a=p.parse_args()
    main(a.cache,a.out,a.workers)
