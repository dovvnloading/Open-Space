"""Scale the same PHOENIX grid to W3 and independently to Ks for the MIRI host."""
from pathlib import Path
import sys,csv,json,hashlib
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from mcc124_sed import get_band,band_integral,get_inputs

def main(cache):
    get_inputs(cache)
    for ref in json.loads((HERE/'bandpasses.json').read_text()):
        path=cache/ref['filename']
        if not path.exists():path.write_bytes(urlopen(ref['url'],timeout=70).read())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=ref['sha256']:raise ValueError(path)
    names=('ks','w3','f1065c','f1550c')
    bands=[get_band(cache/(n+'.xml')) for n in names]
    def photo(name,identifier,field):
        rows=list(csv.DictReader((HERE/(name+'.csv')).open()))
        return float(next(s[field] for s in rows if s['designation']==identifier))
    km=photo('twomass','05312734-0340356','k_m')
    wm=photo('allwise','J053127.91-034059.8','w3mpro')
    with fits.open(cache/'vega.fits') as hd:
        vega=np.array([band_integral(hd[1].data['WAVELENGTH'],hd[1].data['FLUX'],b) for b in bands])
    rows=[]
    for t in (3500,3700,3800,3900):
        with fits.open(cache/f'ph{t}.fits') as hd:
            tab=hd[1].data
            for g in (45,50):
                fl=np.array([band_integral(tab['WAVELENGTH'],tab[f'g{g}'],b) for b in bands])
                mags=-2.5*np.log10(fl/vega)
                r=dict(teff_K=t,logg=g/10)
                for j in (2,3):
                    for k,anchor in enumerate((km,wm)):
                        r[names[j]+'_Jy_from_'+names[k]]=bands[j][3]*10**(-.4*(anchor+mags[j]-mags[k]))
                rows.append(r)
    with (HERE/'photometry.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
    for key in rows[0]:
        if '_Jy_' in key:
            vals=[r[key] for r in rows];print(key,min(vals),np.median(vals),max(vals))
if __name__=='__main__':main(Path(sys.argv[1]))
