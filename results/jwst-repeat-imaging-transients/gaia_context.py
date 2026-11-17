"""Evaluate proximity to Gaia DR3 point sources at the epoch of imaging."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
import pandas as pd

QUERY="""SELECT TOP 10000 source_id,ra,dec,pmra,pmdec,parallax,parallax_error,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',268.46,65.2,0.18))"""

def run(pos,out,cache):
    url='https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=QUERY))
    if not cache.exists(): cache.write_bytes(urlopen(url,timeout=100).read())
    g=pd.read_csv(cache,dtype={'source_id':str})
    # Linear propagation 2016.0 to 2026.6. Missing motions leave the Gaia
    # catalog position at its reference epoch and are marked in the output.
    dec=g.dec+g.pmdec.fillna(0)*10.6/3.6e6
    ra=g.ra+g.pmra.fillna(0)*10.6/3.6e6/np.cos(np.deg2rad(g.dec))
    rows=[]
    for r in pd.read_csv(pos).itertuples():
        sep=np.hypot((r.ra_deg-ra)*np.cos(np.deg2rad(r.dec_deg)),r.dec_deg-dec)*3600
        j=sep.idxmin()
        rows.append(dict(id=r.id,nearest_dr3=g.source_id[j],distance_arcsec=sep[j],pmra_masyr=g.pmra[j],
                         pmdec_masyr=g.pmdec[j],phot_g_mean_mag=g.phot_g_mean_mag[j]))
    pd.DataFrame(rows).to_csv(out,index=False)
    out.with_suffix('.json').write_text(json.dumps(dict(query=QUERY,url=url,rows=len(g),sha256=hashlib.sha256(cache.read_bytes()).hexdigest()),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('positions',type=Path);p.add_argument('out',type=Path);p.add_argument('cache',type=Path)
    a=p.parse_args();run(a.positions,a.out,a.cache)
