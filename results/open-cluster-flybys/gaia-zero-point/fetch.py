"""Refresh Gaia fields for the frozen UCC memberships and recompute zero points."""
import argparse
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import pandas as pd
from zero_point import zpt

ROOT = Path(__file__).resolve().parent
TAP = 'https://gea.esac.esa.int/tap-server/tap/sync'
FIELDS = ('source_id,ra,dec,parallax,parallax_error,pmra,pmdec,radial_velocity,'
          'radial_velocity_error,ruwe,phot_g_mean_mag,nu_eff_used_in_astrometry,'
          'pseudocolour,ecl_lat,astrometric_params_solved')
ZPT_FIELDS = ['phot_g_mean_mag','nu_eff_used_in_astrometry','pseudocolour',
              'ecl_lat','astrometric_params_solved']
CLUSTERS = ['ascc114','ubc168','upk90','upk93','ascc99','upk20','ngc6633','upk41']

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'refreshed.csv.gz')
    args=parser.parse_args()
    m=pd.read_csv(ROOT.parent/'inputs'/'ucc_members.csv',dtype={'Source':str})
    m=m.loc[m.name.isin(CLUSTERS),['name','Source','probs']]
    ids=m.Source.unique()
    blocks=[]
    for start in range(0,len(ids),80):
        subset=','.join(ids[start:start+80])
        query=f'SELECT {FIELDS} FROM gaiadr3.gaia_source WHERE source_id IN ({subset})'
        params=urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query))
        with urlopen(TAP+'?'+params,timeout=90) as response:
            blocks.append(pd.read_csv(BytesIO(response.read()),dtype={'source_id':str}))
    dr3=pd.concat(blocks,ignore_index=True).rename(columns={'source_id':'Source'})
    table=m.merge(dr3,on='Source',validate='many_to_one')
    if len(table)!=len(m):
        raise ValueError('Incomplete Gaia response')
    zpt.load_tables()
    table['zpt_mas']=zpt.get_zpt(*(table[c].to_numpy() for c in ZPT_FIELDS),_warnings=False)
    table.sort_values(['name','Source']).to_csv(args.output,index=False,float_format='%.11g')

if __name__=='__main__':
    main()
