"""Re-fit the archived small cutouts, or retrieve the wide fields for the sky controls."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from screen import ROOT, download
from validate import analyze, fit_field

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fields",type=Path,help="Fetch the original wide cutouts to this cache and repeat the sky comparisons")
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    fits_table=pd.read_csv(ROOT/'beam-fits.csv')
    visits=pd.read_csv(ROOT/'followup-visits.csv')
    published=fits_table.merge(visits,on=['name','epoch'])
    max_error=0.
    for row in published.itertuples(index=False):
        result,_=fit_field(ROOT/'inputs'/f'{row.name}-{row.epoch}.fits',row.ra,row.dec)
        max_error=max(max_error,abs(result['forced']-row.forced))
        if not np.isclose(result['forced'],row.forced,atol=2e-5):
            raise ValueError(f'{row.name} epoch {row.epoch} disagrees with the wide-field fit')
        if args.fields:
            import hashlib
            data=download(row,.065,args.fields/f'{row.name}-{row.epoch}.fits')
            if hashlib.sha256(data).hexdigest()!=row.field_sha256:
                raise ValueError(f'Wide cutout checksum differs for {row.name} epoch {row.epoch}')
    print(f'{len(published)} archived forced fluxes verified. Maximum difference = {max_error:.3g} mJy/beam.')
    if args.fields:
        dest=args.output or args.fields
        dest.mkdir(parents=True,exist_ok=True)
        analyze(visits,args.fields,dest)

if __name__=='__main__':
    main()
