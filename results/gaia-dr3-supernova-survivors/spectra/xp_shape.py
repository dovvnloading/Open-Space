"""Compare observed XP continuum colours to a reddened 20,000 K blackbody.

Downloads the Gaia DR3 GAVO 400–800 nm sampled spectra, Gaia Apsis solutions,
and IRSA all-sightline reddening for every XP-detected Gaia-screen hit.
"""
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).parent
ESA='https://gea.esac.esa.int/tap-server/tap/sync'
GAVO='https://dc.g-vo.org/tap/sync'
IRSA='https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust'
GRID=np.arange(400,801,10)
RV=3.1


def get(url,outfile):
    with urlopen(url,timeout=90) as response:
        data=response.read()
    if outfile:
        (ROOT/outfile).write_bytes(data)
    return data,dict(url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def k_lambda(w):
    # Cardelli, Clayton & Mathis (1989), optical a(x), b(x).
    y=1000/w - 1.82
    a=np.polynomial.polynomial.polyval(y,[1,0.17699,-0.50447,-0.02427,0.72085,0.01979,-0.77530,0.32999])
    b=np.polynomial.polynomial.polyval(y,[0,1.41338,2.28305,1.07233,-5.38434,-0.62251,5.30260,-2.09002])
    return RV*a+b


def colour(f):
    return np.mean(f[(GRID>=430)&(GRID<=480)])/np.mean(f[(GRID>=680)&(GRID<=730)])


def dust(row):
    url=IRSA+'?'+urlencode(dict(locstr=f"{float(row['ra']):.6f} {float(row['dec']):+.6f} equ j2000"))
    raw,receipt=get(url,None)
    doc=ET.fromstring(raw)
    if doc.get('status')!='ok':
        raise ValueError(doc.findtext('message'))
    ebv=float(doc.findtext('./result/statistics/refPixelValueSandF').split()[0])
    return row['source_id'],ebv,receipt


def main():
    rows=[r for r in csv.DictReader((ROOT.parent/'search/gaia.csv').open())
          if r['has_xp_continuous']=='true']
    ids=','.join(sorted(r['source_id'] for r in rows))
    def query(server,table,fields,filename):
        adql=f'SELECT {fields} FROM {table} WHERE source_id IN ({ids}) ORDER BY source_id'
        url=server+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=adql))
        return get(url,filename)
    raw,a=query(GAVO,'gdr3spec.spectra','source_id,flux,flux_error','xp_sampled.csv')
    text=raw.decode()
    fluxes=list(csv.DictReader(io.StringIO(text[text.index('source_id,'): ])))
    raw,b=query(ESA,'gaiadr3.astrophysical_parameters',
                'source_id,teff_gspphot,logg_gspphot,mh_gspphot,azero_gspphot,ag_gspphot,'
                'classprob_dsc_combmod_whitedwarf,classprob_dsc_combmod_star','xp_apsis.csv')
    with ThreadPoolExecutor(3) as pool:
        extinction=list(pool.map(dust,rows))
    ebvs={i:e for i,e,_ in extinction}
    k=k_lambda(GRID)
    # Shape of Planck radiance per unit wavelength. Absolute normalisation cancels.
    hot=GRID.astype(float)**-5/np.expm1(1.438776877e7/(GRID*20000.))
    with (ROOT/'xp_continuum.csv').open('w') as fp:
        out=csv.writer(fp)
        out.writerow(['source_id','blue_red_ratio','ratio_20k_at_map_EBV','EBV_required_20k','EBV_map_to_infinity'])
        for row in fluxes:
            flux=np.asarray(json.loads(row['flux']))
            if len(flux)!=len(GRID):
                raise ValueError('Unexpected XP sampled grid')
            ratio=colour(flux)
            reddened=lambda ebv:colour(hot*10**(-0.4*ebv*k))
            required=brentq(lambda e:reddened(e)-ratio,-1,6)
            out.writerow([row['source_id'],ratio,reddened(ebvs[row['source_id']]),required,ebvs[row['source_id']]])
    if len(fluxes)!=len(rows):
        raise ValueError('Missing GAVO XP product')
    receipts={'sampled':a,'apsis':b,'dust':{i:{**r,'EBV':e} for i,e,r in extinction}}
    (ROOT/'xp_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__':
    main()
