"""Measure native-aperture chromatic consistency on the frozen sparse tracks."""
import argparse
import csv
from pathlib import Path

import numpy as np
from astropy.io import fits

from screen import ROOT,reduce


def aperture(image,dq,xy):
    # Fractional 2.5-pixel aperture, native electrons; use the median in a 5-to-8-pixel annulus.
    ix,iy=np.rint(xy).astype(int)
    patch=image[iy-8:iy+9,ix-8:ix+9]
    if patch.shape!=(17,17):return np.nan
    flags=dq[iy-8:iy+9,ix-8:ix+9];yy,xx=np.mgrid[-8:9,-8:9]
    xx=xx+ix-xy[0];yy=yy+iy-xy[1]
    inside=np.zeros(patch.shape)
    for u in [-.375,-.125,.125,.375]:
        for v in [-.375,-.125,.125,.375]:inside+=(xx+u)**2+(yy+v)**2<=2.5**2
    inside/=16.
    valid=(flags==0)&np.isfinite(patch)
    if (inside*valid).sum()<.95*inside.sum():return np.nan
    outer=(xx**2+yy**2>25)&(xx**2+yy**2<=64)&valid
    if outer.sum()<50:return np.nan
    return ((np.where(valid,patch,0)-np.median(patch[outer]))*inside*valid).sum()*inside.sum()/(inside*valid).sum()


def run(tag,cache):
    path=ROOT/'epochs'/tag
    obs=list(csv.DictReader((path/'sparse_vet_exposures.csv').open()))
    result=[]
    for name in sorted({r['rootname'] for r in obs}):
        sel=[r for r in obs if r['rootname']==name]
        with fits.open(cache/(name+'_flt.fits')) as hdus:
            hdu=next(s for s in hdus if s.name=='SCI' and s.header['CCDCHIP']==int(sel[0]['ccdchip']))
            image=hdu.data;dq=hdus['DQ',hdu.header['EXTVER']].data;hdr=hdu.header
            exposure=hdus[0].header['EXPTIME']
            zp=-2.5*np.log10(hdr['PHOTFLAM']*hdr['PHOTPLAM']**2/2.99792458e18)-48.60
            off=np.concatenate([np.column_stack([r*np.cos(t),r*np.sin(t)])
                                for r,n in [(13,12),(22,20),(31,28)]
                                for t in [np.arange(n)*2*np.pi/n]])
            for obsrow in sel:
                xy=np.array([float(obsrow['physical_x'])+hdr.get('LTV1',0),float(obsrow['physical_y'])+hdr.get('LTV2',0)])
                electrons=aperture(image,dq,xy)
                null=np.array([aperture(image,dq,xy+o) for o in off]);null=null[np.isfinite(null)]
                if len(null)>=15:
                    med=np.median(null);sigma=1.4826*np.median(abs(null-med));net=electrons-med
                else:sigma=net=np.nan
                fluxerr=np.sqrt(sigma**2+max(0,net)) if np.isfinite(net) else np.nan
                mag=zp-2.5*np.log10(net/exposure) if net>0 else np.nan
                lim=zp-2.5*np.log10(5*sigma/exposure) if sigma>0 else np.nan
                result.append(dict(template=obsrow['template'],filter=obsrow['filter'],rootname=name,visit=obsrow['visit'],
                        orbit_support=obsrow['supports'],exptime_s=exposure,AB_zp=zp,comparison_apertures=len(null),
                        net_electrons=net,background_MAD_electrons=sigma,poisson_plus_background_error=fluxerr,
                        aperture_AB=mag,background_5MAD_aperture_AB=lim))
    reduce.write_csv(path/'sparse_vet_apertures.csv',result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('epoch',nargs='+')
    a=p.parse_args()
    for t in a.epoch:run(t,a.cache)
