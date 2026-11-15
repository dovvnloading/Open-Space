"""Repeat the archived TESS-FFI target list and blind period search.

Usage: python repeat_screen.py --cache /path/for/downloads > replay.csv
"""
import argparse
import csv
from hashlib import sha256
from io import BytesIO
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile
import numpy as np
from scipy.ndimage import median_filter
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.timeseries import LombScargle

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cache',type=Path,required=True)
    cache=parser.parse_args().cache;cache.mkdir(parents=True,exist_ok=True)
    cat=cache/'catalog.fits'
    if not cat.exists():
        cat.write_bytes(urlopen('https://zenodo.org/records/18842123/files/amcvn_catalogue.fits?download=1',timeout=90).read())
    if sha256(cat.read_bytes()).hexdigest()!='2b11c042119bf0ed8a0d257ee07fe35488ff5ebe96462ee4797aae30929b96e3':
        raise ValueError('Catalog differs from the screening release')
    catalog=Table.read(cat)
    writer=csv.writer(sys.stdout);writer.writerow(['catalog_row','sector','n','minutes','peak_to_median'])
    for row in csv.DictReader(open('field-screen.csv')):
        i,sector=int(row['catalog_row']),int(row['sector'])
        s=catalog[i];ra,dec=float(s['RA']),float(s['Dec'])
        path=cache/f'{i}_{sector}.fits'
        if not path.exists():
            # The exact original query for two reference cutouts used four decimal places.
            digits=4 if (i,sector) in [(101,62),(101,89)] else 7
            q=urlencode(dict(ra=round(ra,digits),dec=round(dec,digits),x=7,y=7,units='px',sector=sector))
            raw=urlopen('https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+q,timeout=150).read()
            with ZipFile(BytesIO(raw)) as z:path.write_bytes(z.read(next(v for v in z.namelist() if v.endswith('.fits'))))
        with fits.open(path) as h:
            source=WCS(h[2].header).world_to_pixel_values(ra,dec)
            d=h[1].data;t=d['TIME'].astype(float);cube=d['FLUX'].astype(float)
            ok=(d['QUALITY']==0)&np.isfinite(t)&np.isfinite(cube).all(axis=(1,2))
            t,cube=t[ok],cube[ok]
        yy,xx=np.indices(cube.shape[1:]);rr=(xx-source[0])**2+(yy-source[1])**2
        bg=np.median(cube[:,rr>=9],axis=1)
        flux=(cube[:,rr<4]-bg[:,None])@np.exp(-rr[rr<4]/(2*.8**2))
        y=flux-median_filter(flux,size=73,mode='nearest')
        ok=abs(y-np.median(y))<5*1.4826*np.median(abs(y-np.median(y)))
        t,y=t[ok],y[ok]
        hz=min(144,.495/np.median(np.diff(t)))
        f,p=LombScargle(t,y).autopower(minimum_frequency=18,maximum_frequency=hz,samples_per_peak=7,method='fast')
        # The inventory excludes boundary bins from its peak selection.
        from scipy.signal import find_peaks
        peaks,_=find_peaks(p,distance=7);k=peaks[np.argmax(p[peaks])]
        writer.writerow([i,sector,len(t),1440/f[k],p[k]/np.median(p)])
if __name__=='__main__':main()
