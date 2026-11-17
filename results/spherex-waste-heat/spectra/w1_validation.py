"""Approximate W1 synthetic flux check for A–J, without applying source-specific colour corrections."""
import csv
import gzip
import hashlib
import io
import json
import sys
from pathlib import Path
from urllib.request import urlopen
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'wise'))
from fetch import position,separation

def main():
    url='https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/figures/RSR-W1.txt'
    with urlopen(url,timeout=45) as response: raw=response.read()
    wavelength,resp,_=np.loadtxt(io.BytesIO(raw)).T
    weight=resp/wavelength
    (ROOT/'w1_response.json').write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(raw).hexdigest()),indent=2)+'\n')
    spec=list(csv.DictReader((ROOT/'channel_spectra.csv').open()))
    result=[]
    for target in csv.DictReader((ROOT/'positions.csv').open()):
        label=target['label']
        if '_cmp' in label: continue
        with gzip.open(ROOT.parent/'wise/archive'/f'{label}_allwise.csv.gz','rt') as src:
            ranked=sorted((separation(float(r['ra']),float(r['dec']),*position(target,2010.5)),r) for r in csv.DictReader(src))
        cat_mjy=309540*10**(-.4*float(ranked[0][1]['w1mpro'])) if ranked and ranked[0][0]<=2 else None
        # A diagnostic band average for approximately flat Fnu; no W1 response extrapolation outside DR1.
        source=[r for r in spec if r['label']==label and r['outlier_included_mJy']]
        wl=np.array([float(r['wavelength_um']) for r in source])
        flux=np.array([float(r['outlier_included_mJy']) for r in source])
        grid=(wavelength>=min(wl))&(wavelength<=max(wl))
        num=np.trapezoid(np.interp(wavelength[grid],wl,flux)*weight[grid],wavelength[grid])
        den=np.trapezoid(weight[grid],wavelength[grid])
        spx=num/den
        result.append(dict(label=label,spx_W1avg_mJy=round(spx,5),allwise_W1_mJy=round(cat_mjy,5) if cat_mjy else '',
                           spx_over_allwise=round(spx/cat_mjy,4) if cat_mjy else '',
                           response_covered=round(den/np.trapezoid(weight,wavelength),5)))
    with (ROOT/'w1_validation.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=result[0]); w.writeheader(); w.writerows(result)

if __name__=='__main__': main()
