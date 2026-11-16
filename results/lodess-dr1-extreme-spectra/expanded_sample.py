"""Relax the parent LoLSS selection and verify the 24 largest corrected in-band slopes.

Reuses the primary LoLSS and comparison scripts. Distances are in arcseconds except
the good-field radius (degrees). LoLSS input fluxes are mJy, comparison output Jy.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

from comparison.fetch_lotss import fetch
from comparison.match_catalogues import lotss_match, survey_match
from lolss_sample import FREQUENCIES, fit_slope, nearest, read_table
from release.photometry import RadioImage


def run(release, lolss_cache, comparison_cache, out):
    out.mkdir(parents=True, exist_ok=True)
    src = read_table(lolss_cache/'lolss.fits')
    fields = pd.read_csv(release/'fields.csv').query("quality=='good'").rename(
        columns={'ra_deg':'RA','dec_deg':'DEC'})
    d, ix = nearest(src,fields)
    src['field'] = fields.field.to_numpy()[ix[:,0]]
    src['radius'] = d[:,0]/3600
    d, _ = nearest(src,src)
    src['nn'] = d[:,1]
    ratio = src.Total_flux/src.Peak_flux
    src['new'] = (src.Total_flux < 100) | (ratio >= 1.35)
    src = src[(src.radius<=3.5) & (src.Total_flux>=60) & (src.nn>=60) &
              (ratio<1.65) & (src.S_Code=='S')].copy()
    n_initial = len(src)
    for i, f in enumerate(FREQUENCIES):
        ch = read_table(lolss_cache/f'chan{i}.fits')
        d, ix = nearest(src,ch)
        for col, key in [('Total_flux','s'),('E_Total_flux','e'),('Peak_flux','pk'),('Isl_rms','rms')]:
            src[f'{key}{f}'] = np.where(d[:,0]<5,ch[col].to_numpy()[ix[:,0]],np.nan)
    src['all_snr'] = np.min([src[f'pk{f}']/src[f'rms{f}'] for f in FREQUENCIES],axis=0)
    src = src[src.all_snr>=8].copy()
    src[['slope','eslope','chisq']] = np.array([fit_slope(r) for _,r in src.iterrows()])
    src['stat'] = (src.slope-.27)/src.eslope
    src.sort_values('stat',ascending=False,inplace=True)
    src.index.name = 'lolrow'
    src.to_csv(out/'ranked.csv.gz')
    print(f'{n_initial} parent sources, {len(src)} with usable six-band fits')
    sample = src.head(24).reset_index()
    coords = sample[['RA','DEC']].rename(columns={'RA':'ra','DEC':'dec'})
    fetch(Table.from_pandas(sample[['RA','DEC']]),0,out)
    lotss_match(out,coords,out)
    for survey in ['gmrtas150m','vlssr','wenss','nvss']:
        survey_match(survey,comparison_cache,coords,out)
    for f in sorted(out.glob('match_*.csv.gz')):
        sample = sample.join(pd.read_csv(f,index_col='row'))
    phot = {}
    for field, group in sample.groupby('field'):
        with RadioImage(release/'images'/(field+'-I.fits')) as image:
            for i,r in group.iterrows():
                phot[i] = image.measure(r.RA,r.DEC) or {}
    sample.join(pd.DataFrame.from_dict(phot,orient='index')).to_csv(out/'top24.csv',index=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['release','lolss_cache','comparison_cache','output']:
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    run(args.release,args.lolss_cache,args.comparison_cache,args.output)
