"""Recover source and Gaussian parent fields, validating against image WCS and astrometry.

Usage: python map_fields.py DATA_DIRECTORY
Install NumPy, pandas, SciPy, and Astropy. Row numbers are zero-based FITS table row indices.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u


def blocks(ids):
    bounds = np.r_[0, np.where(np.diff(ids) < 0)[0] + 1, len(ids)]
    return list(zip(bounds[:-1], bounds[1:]))


def coordinates(rows):
    return SkyCoord(rows['RA'] * u.deg, rows['DEC'] * u.deg)


def associate(root):
    qa = pd.read_csv(root/'fields.csv').set_index('field')
    sources = fits.getdata(root/'sources.fits')
    gauss = fits.getdata(root/'gaussians.fits')
    wcs = {f: WCS(fits.Header.fromfile(root/'headers'/(f+'.hdr'))).celestial for f in qa.index}
    gblocks = {}
    gfield = np.full(len(gauss), '', dtype='U10')
    checks = []
    for first, last in blocks(gauss['Gaus_id']):
        block = gauss[first:last]
        if not np.array_equal(block['Gaus_id'], np.arange(len(block))):
            raise ValueError(f'Noncontiguous Gaussian IDs at row {first}')
        sample = block[::max(1, len(block)//100)]
        residual = {}
        for field, w in wcs.items():
            x, y = w.world_to_pixel_values(sample['RA'], sample['DEC'])
            d = np.hypot(x-sample['Xposn'], y-sample['Yposn'])
            residual[field] = np.nanmedian(d) if np.isfinite(d).any() else np.inf
        best, second = sorted(residual, key=residual.get)[:2]
        if residual[best] > 3 or residual[second] < 100 or best in gblocks:
            raise ValueError(f'Ambiguous field for Gaussian row {first}: {best} {residual[best]}, {second} {residual[second]}')
        gfield[first:last] = best
        gblocks[best] = block
        checks.append(dict(field=best, gauss_start=first, gauss_end_exclusive=last,
                           wcs_median_residual_pix=residual[best], next_wcs_median_pix=residual[second]))
    if set(gblocks) != set(qa.index):
        raise ValueError('Gaussian catalogue and release pointings differ')
    output = pd.DataFrame({'row': np.arange(len(sources)), 'field': '', 'separation_deg': np.nan})
    for first, last in blocks(sources['Source_id']):
        block = sources[first:last]
        single = block[block['S_Code'] == 'S']
        if not len(single):
            raise ValueError(f'No simple sources in block {first}')
        sky = coordinates(single)
        matches = {}
        for field, gs in gblocks.items():
            # A simple source should have exactly one Gaussian component with equal ID and position.
            idx = np.argsort(gs['Source_id'])
            ids = gs['Source_id'][idx]
            pos = np.minimum(np.searchsorted(ids, single['Source_id']), len(ids)-1)
            gg = gs[idx[pos]]
            valid = (gg['Source_id'] == single['Source_id']) & (gg['Isl_id'] == single['Isl_id'])
            matches[field] = int(np.sum(valid & (sky.separation(coordinates(gg)).arcsec < 1)))
        best, second = sorted(matches, key=matches.get, reverse=True)[:2]
        if matches[best]/len(single) < .95 or matches[second] > .05*len(single):
            raise ValueError(f'Unvalidated source block {first}: {best},{matches[best]}/{len(single)}; {second}, {matches[second]}')
        q = qa.loc[best]
        output.loc[first:last-1, 'field'] = best
        output.loc[first:last-1, 'separation_deg'] = coordinates(block).separation(SkyCoord(q.ra_deg*u.deg, q.dec_deg*u.deg)).deg
        check = next(c for c in checks if c['field'] == best)
        if 'source_start' in check:
            raise ValueError(f'Duplicate source block for {best}')
        check.update(source_start=first, source_end_exclusive=last, singles=len(single),
                     matched_singles=matches[best], other_max_matches=matches[second])
    output['quality'] = output['field'].map(qa['quality'])
    output.to_csv(root/'source_fields.csv', index=False)
    checks = pd.DataFrame(checks).join(qa['quality'], on='field')
    checks.to_csv(root/'field_membership.csv', index=False, float_format='%.6f')
    pd.DataFrame({'row': np.arange(len(gauss)), 'field': gfield}).to_csv(root/'gaussian_fields.csv', index=False)
    print(output['quality'].value_counts(), checks[['wcs_median_residual_pix','next_wcs_median_pix']].agg(['min','max']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    associate(parser.parse_args().root)
