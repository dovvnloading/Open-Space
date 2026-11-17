"""PSF-matched forced surface photometry of coherent positive Euclid ridges.

These measurements integrate the same central 70% of each VIS segment in all
bands. They include unresolved objects, unrelated field light and optical
artefacts when present. Errors use translated apertures in each image.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi

from acquire import write_csv

FWHM={'VIS':.17,'Y':.30,'J':.32,'H':.36}
TARGET=.45


def observe(a,mask,p0,p1,scale,zeropoint,band):
    # Apply Gaussian PSF matching and remove aperture samples which touch holes.
    kernel=np.sqrt(TARGET**2-FWHM[band]**2)/(2.355*scale)
    sm=ndi.gaussian_filter(np.where(mask,a,0.),kernel,truncate=3.)
    sm[~ndi.minimum_filter(mask.astype('u1'),size=2*int(np.ceil(3*kernel))+1).astype(bool)]=np.nan
    dp=p1-p0;v=dp/np.linalg.norm(dp);perp=np.array([-v[1],v[0]])
    t=np.linspace(.15,.85,51)
    transverse=np.array([-.20,-.10,0,.10,.20,-.80,-.70,-.60,-.50,-.40,.40,.50,.60,.70,.80])/scale
    on=p0+t[:,None]*dp
    traces=on[:,None,:]+transverse[None,:,None]*perp
    offsets=np.array([(x,y) for x in (-5,-3,-1,1,3,5) for y in (-5,-3,-1,1,3,5)
                      if abs(np.dot([x,y],perp))>1.3])/scale
    grid=traces[None,:,:,:]+np.vstack([[0,0],offsets])[:,None,None,:]
    z=ndi.map_coordinates(sm,[grid[...,1].ravel(),grid[...,0].ravel()],order=1,mode='constant',cval=np.nan).reshape(grid.shape[:-1])
    contrast=z[:,:,:5].mean(axis=(1,2))-z[:,:,5:].mean(axis=(1,2))
    blank=contrast[1:][np.isfinite(contrast[1:])]
    out=dict(usable=0,blank_apertures=len(blank),flux_njy_arcsec2='',noise_njy_arcsec2='',sn='',mu_ab_arcsec2='')
    if len(blank)<12 or not np.isfinite(contrast[0]): return out
    factor=3631e9*10**(-.4*zeropoint)/(scale**2)
    bg=np.median(blank);fl=(contrast[0]-bg)*factor;sig=np.median(np.abs(blank-bg))*1.4826*factor
    if not sig>0:return out
    out.update(usable=1,flux_njy_arcsec2=float(fl),noise_njy_arcsec2=float(sig),sn=float(fl/sig),
               mu_ab_arcsec2=float(-2.5*np.log10(fl/(3631e9))) if fl>0 else '')
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['candidates','metrics','images','cache','out']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    byid={r['id']:r for r in csv.DictReader(a.candidates.open())}
    ids=[m['id'] for m in csv.DictReader(a.metrics.open()) if m['sign']=='positive' and m['band']=='VIS' and m['usable']=='1'
         and int(m['positive_bins'])>=4 and float(m['central_ridge_contrast_sigma'])>=3]
    rows=[];colors=[]
    for ident in ids:
        r=byid[ident];selected={}
        for band in FWHM:
            fp=a.images/(r['tile']+'_'+band+'.fits')
            if not fp.exists():fp=a.cache/(ident+'_'+band+'.fits')
            if not fp.exists():continue
            with fits.open(fp) as hdus:
                h=hdus[0].header;arr=hdus[0].data;w=WCS(h);mask=np.isfinite(arr)&(arr!=0)
                if band=='VIS':
                    with fits.open(a.images/(r['tile']+'_rms.fits')) as mf:
                        rms=mf[0].data;mask &= np.isfinite(rms)&(rms>0)&(rms<1e8)
                p0,p1=w.all_world2pix([[float(r[f'ra{i}_deg']),float(r[f'dec{i}_deg'])] for i in (0,1)],0)
                scale=3600*np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))
                m=observe(arr,mask,p0,p1,scale,float(h['MAGZERO']),band)
            selected[band]=m;rows.append(dict(id=ident,band=band,**m))
        col=dict(id=ident)
        for b0,b1 in [('VIS','Y'),('Y','J'),('J','H'),('VIS','H')]:
            key=b0+'_'+b1+'_AB';col[key]='';col[key+'_error']=''
            v0,v1=selected.get(b0,{}),selected.get(b1,{})
            if v0.get('usable') and v1.get('usable') and min(v0['sn'],v1['sn'])>=3:
                col[key]=v0['mu_ab_arcsec2']-v1['mu_ab_arcsec2']
                col[key+'_error']=1.085736*np.hypot(1/v0['sn'],1/v1['sn'])
        colors.append(col)
    if rows:write_csv(a.out/'matched_surface_photometry.csv',rows)
    if colors:write_csv(a.out/'matched_colours.csv',colors)


if __name__=='__main__':main()
