"""Extract SPHEREx fluxes and compare disjoint wavelength-matched triplets."""
import csv
from functools import lru_cache
import hashlib
import io
from itertools import combinations
from pathlib import Path
import zipfile

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from replay import DEC, RA, ROOT, tablewrite

PIXEL_ARCSEC = 6.15
JY_PER_PIXEL = 1e6 * (PIXEL_ARCSEC / 206264.806) ** 2
OFFSET = 5.5 * PIXEL_ARCSEC / 3600
COORDS = [[RA, DEC], [RA, DEC + OFFSET], [RA, DEC - OFFSET],
          [RA + OFFSET / np.cos(np.deg2rad(DEC)), DEC],
          [RA - OFFSET / np.cos(np.deg2rad(DEC)), DEC]]
CONTROLS = ('north', 'south', 'east', 'west')


def extract(content, obs_id, detector):
    with fits.open(io.BytesIO(content), memmap=False) as hdus:
        image = hdus['IMAGE'].data.astype(float)
        var = hdus['VARIANCE'].data
        flags = hdus['FLAGS'].data
        h = hdus['IMAGE'].header
        if h['OBSID'] != obs_id or h['DETECTOR'] != detector:
            raise ValueError(f'Wrong spectrum: {obs_id} D{detector}')
        xy = WCS(h, hdus).all_world2pix(COORDS, 0)
        x, y = xy[0]
        wave, width = WCS(h, hdus, key='W').wcs_pix2world([[x, y]], 0)[0]
        yy, xx = np.indices(image.shape)
        radius = np.hypot(xx - x, yy - y)
        good = ((flags & ~(1 << 21)) == 0) & np.isfinite(image) & np.isfinite(var) & (var > 0)
        sky = image[(radius > 4) & (radius < 7) & good]
        bg = np.median(sky) if len(sky) else np.nan
        scatter = 1.4826 * np.median(np.abs(sky - bg)) if len(sky) else np.nan
        # Reject an incomplete aperture instead of rescaling missing pixels.
        def aperture(cx, cy, r):
            mask = np.hypot(xx-cx, yy-cy) < r
            complete = (r-.5 <= cx <= image.shape[1]-.5-r and r-.5 <= cy <= image.shape[0]-.5-r)
            if not complete or not np.all(good[mask]) or len(sky) < 40:
                return np.nan, np.nan
            n = mask.sum()
            flux = JY_PER_PIXEL * (image[mask].sum() - n * bg)
            error = JY_PER_PIXEL * np.sqrt(var[mask].sum(dtype=float) + (n*1.253314*scatter)**2 / len(sky))
            return flux, error

        mjd = h['MJD-AVG']
        visit = 1 if mjd < 61000 else 2 if mjd < 61200 else 3
        ap = radius < 2.5
        result = dict(detector=detector, obs_id=obs_id, visit=visit, mjd=mjd,
                      wavelength_um=wave, bandwidth_um=width, sky_pixels=len(sky),
                      aperture_pixels=ap.sum(), bad_aperture_pixels=(ap & ~good).sum())
        result['flux_jy'], result['formal_error_jy'] = aperture(x,y,2.5)
        for field, r in [('core_jy',1.5),('small_jy',2.),('large_jy',3.)]:
            result[field] = aperture(x,y,r)[0]
        for name, (cx,cy) in zip(CONTROLS,xy[1:]):
            result[name+'_jy'] = aperture(cx,cy,2.5)[0]
        return result


def match(frames, detector, tolerance, field='flux_jy'):
    visits = [sorted((r for r in frames if r['detector'] == detector and r['visit'] == visit
                      and np.isfinite(r[field]) and r['bandwidth_um'] > 0),
                     key=lambda r: (r['wavelength_um'], r['obs_id'])) for visit in (1,2,3)]

    @lru_cache(None)
    def align(i,j,k):
        if any(n == len(v) for n,v in zip((i,j,k),visits)):
            return 0, 0., ()
        # Order-preserving alignment. Maximize disjoint triplets, breaking ties
        # by the sum of their largest pair separations in units of bandwidth.
        choices = [align(i+1,j,k), align(i,j+1,k), align(i,j,k+1)]
        trio = (visits[0][i],visits[1][j],visits[2][k])
        cost = max(abs(a['wavelength_um']-b['wavelength_um']) /
                   min(a['bandwidth_um'],b['bandwidth_um']) for a,b in combinations(trio,2))
        if cost <= tolerance:
            n, penalty, tail = align(i+1,j+1,k+1)
            choices.append((n+1,penalty-cost,((i,j,k),)+tail))
        return max(choices,key=lambda c: (c[0],c[1]))

    return [[v[n] for v,n in zip(visits,indices)] for indices in align(0,0,0)[2]]


def write_csv(path, rows):
    def cell(x):
        if isinstance(x,(float,np.floating)):
            return f'{x:.10g}' if np.isfinite(x) else ''
        return x
    keys = list(rows[0])
    tablewrite(path,keys,([cell(r[k]) for k in keys] for r in rows))


def summarize(frames, write=True):
    matches, summary = [], []
    for field, tolerance in [('flux_jy',0.5),('flux_jy',0.25),('core_jy',0.5)]:
        for d in range(1,7):
            trips = match(frames,d,tolerance,field)
            row = dict(detector=d, aperture_flux_column=field, tolerance_bandwidth=tolerance, triplets=len(trips))
            for v in (1,2,3):
                subset = [f for f in frames if f['detector']==d and f['visit']==v]
                row[f'all_v{v}'] = len(subset)
                row[f'clean_v{v}'] = sum(np.isfinite(r[field]) for r in subset)
            flux = np.array([[r[field] for r in trio] for trio in trips]).reshape(-1,3)
            row['rise_count'] = sum(flux[:,1]>flux[:,0])
            row['fade_count'] = sum(flux[:,1]>flux[:,2])
            wave = [r['wavelength_um'] for trio in trips for r in trio]
            row['min_um'] = min(wave,default=np.nan)
            row['max_um'] = max(wave,default=np.nan)
            for v in [1,2]:
                delta = flux[:,v]-flux[:,0]
                row[f'median_v{v+1}_minus_first_jy'] = np.median(delta) if len(delta) else np.nan
                # Ratios are unstable across zero. Publish signed differences there.
                row[f'median_v{v+1}_over_first'] = np.median(flux[:,v]/flux[:,0]) if len(flux) and np.all(flux[:,0]>0) else np.nan
            for c in CONTROLS:
                key = c+'_jy'
                delta = [t[1][key]-t[0][key] for t in trips if np.isfinite(t[1][key]) and np.isfinite(t[0][key])]
                row[f'{c}_pairs'] = len(delta)
                row[f'{c}_middle_minus_first_jy'] = np.median(delta) if delta else np.nan
            shape_checked = [t for t in trips if all(r['flux_jy'] > 0 and
                             0.75 <= r['small_jy'] / r['flux_jy'] <= 1.25 for r in t)]
            row['stable_aperture_triplets'] = len(shape_checked)
            for i in [1, 2]:
                ratios = [t[i]['flux_jy']/t[0]['flux_jy'] for t in shape_checked]
                row[f'stable_aperture_v{i+1}_over_first'] = np.median(ratios) if ratios else np.nan
            summary.append(row)
            if field=='flux_jy' and tolerance==0.5:
                for t in trips:
                    record = dict(detector=d)
                    for label, item in zip(('first','middle','late'),t):
                        for suffix, key in [('id','obs_id'),('um','wavelength_um'),('jy','flux_jy'),('formal_error_jy','formal_error_jy')]:
                            record[label+'_'+suffix] = item[key]
                    matches.append(record)
    if write:
        write_csv(ROOT/'six-detector-flux.csv',frames)
        write_csv(ROOT/'six-detector-matches.csv',matches)
        write_csv(ROOT/'six-detector-summary.csv',summary)
    return summary


def load_frames():
    with (ROOT/'spectral-inputs.csv').open() as handle:
        manifest = list(csv.DictReader(handle))
    from contextlib import ExitStack
    with ExitStack() as stack:
        archives = {d: stack.enter_context(zipfile.ZipFile(ROOT/(f'spectral-cutouts-d{d}.zip' if d<6 else 'spherex-cutouts.zip'))) for d in range(1,7)}
        for r in manifest:
            d,obs = int(r['detector']),r['obs_id']
            filename = f'{obs}D{d}.fits' if d<6 else obs+'.fits'
            data = archives[d].read(filename)
            if hashlib.sha256(data).hexdigest()!=r['cutout_sha256']:
                raise ValueError(f'Checksum mismatch for {filename}')
            yield extract(data,obs,d)


if __name__ == '__main__':
    from astropy import log
    log.setLevel('ERROR')
    summarize(list(load_frames()))
