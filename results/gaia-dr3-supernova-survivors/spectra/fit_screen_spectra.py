"""Independently check Balmer absorption in survey spectra of two Gaia-screen targets.

Needs numpy, scipy and astropy. Wavelengths must be vacuum, heliocentric or
barycentric survey extractions to compare their apparent line shifts.
"""
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

from astropy.io import fits
import numpy as np
from fit_control import velocity_fit, vacuum, factor

ROOT=Path(__file__).parent
CACHE=ROOT/'.cache'
CACHE.mkdir(exist_ok=True)
FILES=[
    ('2588254701112130048','sdss_2588','https://dr20.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/spectra/daily/lite/102XXX/102260/60645/spec-102260-60645-63050395067745527.fits'),
    ('2588254701112130048','sdssold_2588','https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/1898/spec-1898-53260-0343.fits'),
    ('2588254701112130048','lamost_2588','https://www.lamost.org/dr12/v2.0/spectrum/fits/1185006031'),
    ('297028668197924352','sdss_297','https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/7683/spec-7683-58133-0772.fits')]


def fetch(item):
    _,key,url=item
    p=CACHE/key
    if not p.exists():
        with urlopen(url,timeout=90) as response:
            p.write_bytes(response.read())
    raw=p.read_bytes()
    return p,dict(url=url,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())


def read_spectrum(p,lamost):
    with fits.open(p) as hdus:
        data=hdus[1].data
        if lamost:
            w=data['WAVELENGTH'][0].astype(float)
            f=data['FLUX'][0].astype(float)
            iv=data['IVAR'][0].astype(float)
            mask=data['ORMASK'][0]
            cls=hdus[0].header.get('SUBCLASS')
        else:
            w=10.**data['loglam'].astype(float)
            f=data['flux'].astype(float)
            iv=data['ivar'].astype(float)
            mask=data['or_mask']
            cls=hdus[2].data[0]['SUBCLASS']
    return np.column_stack([w,f,np.zeros_like(w),iv,mask]),str(cls)


def pseudo_ew(data,rest,rv):
    # Local continuum either side of Halpha. The observed spectrum is
    # logarithmically sampled; integrate each pixel's measured width.
    center=vacuum(rest)*factor(rv)
    w,f,_,iv,mask=data.T
    good=(mask==0)&np.isfinite(f)&(iv>0)
    x=w-center
    side=good&(np.abs(x)>=18)&(np.abs(x)<=30)
    line=good&(np.abs(x)<=10)
    dx=np.gradient(w)[line]
    if dx.sum()<18 or side.sum()<10:
        return None
    des=np.column_stack((np.ones_like(w),x))
    s=des[side]*np.sqrt(iv[side,None])
    pars=np.linalg.lstsq(s,f[side]*np.sqrt(iv[side]),rcond=None)[0]
    c=des[line]@pars
    grad=(dx*f[line]/c**2)@des[line]
    residual=(des[side]@pars-f[side])*np.sqrt(iv[side])
    var=(dx**2/(iv[line]*c**2)).sum()+grad@np.linalg.inv(s.T@s)@grad
    sigma=np.sqrt(var*max(1.,np.sum(residual**2)/(len(residual)-2)))
    return np.sum(dx*(1-f[line]/c)),sigma


def main():
    receipts={}
    fits_out=[]
    for item in FILES:
        gaia,key,url=item
        p,receipt=fetch(item)
        data,classification=read_spectrum(p,key.startswith('lamost'))
        # Widths of 0.45–5 A fit the Balmer cores only. Wings differ by gravity.
        lines=[4340.472,4861.333,6562.800]
        good=(data[:,4]==0)&(data[:,3]>0)&np.isfinite(data).all(axis=1)
        lines=[line for line in lines if np.sum(good&(abs(data[:,0]-vacuum(line))<13))>=12]
        fit=velocity_fit([data],lines)
        report_line=lines[-1]
        ew=pseudo_ew(data,report_line,fit['rv'])
        fits_out.append({'source_id':gaia,'spectrum':key,'pipeline_subclass':classification,
                         'fitted_air_A':'/'.join(map(str,lines)),**fit,'EW_line_air_A':report_line,
                         'pEW_A':ew[0] if ew else None, 'pEW_sigma_A':ew[1] if ew else None})
        receipts[key]=receipt
    with (ROOT/'screen_balmer.csv').open('w') as fp:
        writer=csv.DictWriter(fp,fieldnames=fits_out[0].keys())
        writer.writeheader()
        writer.writerows(fits_out)
    (ROOT/'screen_spectra_receipt.json').write_text(json.dumps(receipts,indent=2)+'\n')


if __name__=='__main__':
    main()
