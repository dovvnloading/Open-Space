"""Localize a ranked contrast in the source's coadded image and shifted apertures."""
import argparse
import csv
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from acquire import ROOT,REPO,TILE
sys.path.insert(0,str(REPO/'results/spherex-waste-heat/spectra'))
from read_cubes import unpack
from photometry import dist_arcsec,MJYSR_TO_MJY
from extract import write

def contrast_image(st,ch,width,allowed):
    c=ch-1; r=4 if width==1 else 6;gap=r-2
    left=np.arange(c-r,c-gap+1); right=np.arange(c+gap,c+r+1)
    outer=np.r_[left,right];inner=np.arange(c-width//2,c+width//2+1)
    if any(st['DETECTOR'][outer]!=st['DETECTOR'][c]): raise ValueError('window outside detector')
    img=st['IMAGE'].astype(float)
    use=(st['NHITS']>0)&np.isfinite(img)&((st['FLAGS'].astype('u4')&np.uint32(0xffffffff^allowed))==0)
    y=np.where(use[outer],img[outer],0)
    mask=use[outer];valid=(mask[:3].sum(0)>=2)&(mask[3:].sum(0)>=2)&use[inner].all(0)
    # Ordinary per-pixel least squares evaluated at the mean feature wavelength.
    delta=(st['WAVELENGTH'][outer]-np.mean(st['WAVELENGTH'][inner]))[:,None,None]
    n=mask.sum(0); sx=(mask*delta).sum(0); sxx=(mask*delta**2).sum(0)
    sy=y.sum(0);sxy=(y*delta).sum(0)
    denom=n*sxx-sx**2
    with np.errstate(divide='ignore',invalid='ignore'):
        cont=np.where(valid,(sxx*sy-sx*sxy)/denom,np.nan)
    residual=np.mean(img[inner],axis=0)-cont
    return cont,residual

def run(rankfile):
    with rankfile.open() as f:ranks={t['source_id']:t for t in csv.DictReader(f)}
    with (ROOT/'sources.csv').open() as f:sources={t['source_id']:t for t in csv.DictReader(f)}
    selected={k:t for k,t in ranks.items() if int(t['template_rank'])<=5 or t['peer_rank'] in ('1','1.0')}
    rows=[]
    def record(row,st):
        key=row['label'];t=selected[key];ch=int(t['peak_channel']); width=int(t['peak_fwhm_ch'])
        x,y=st['xy']; xx0,yy0=st['origin'];ny,nx=st['IMAGE'].shape[1:]
        pixx,pixy=np.meshgrid(np.arange(nx)+xx0,np.arange(ny)+yy0)
        w=WCS(fits.Header.fromstring(str(st['wcs']),sep=''));w.wcs.bounds_check(False,False)
        ra,dec=w.pixel_to_world_values(x,y)
        lon,lat=w.pixel_to_world_values(pixx,pixy)
        bins=6
        xx,yy=np.meshgrid((np.arange(nx*bins)+.5)/bins+xx0-.5,(np.arange(ny*bins)+.5)/bins+yy0-.5)
        l,b=w.pixel_to_world_values(xx,yy)
        def weights(cx,cy):
            cr,cd=w.pixel_to_world_values(cx,cy); rad=dist_arcsec(l,b,cr,cd)
            def aperture(lo,hi):return ((rad>=lo)&(rad<hi)).reshape(ny,bins,nx,bins).mean((1,3))
            return aperture(0,19),aperture(0,32),aperture(36,60)>.5
        offsets=[(0,0)]+[(dx,dy) for dx in (-10,-5,0,5,10) for dy in (-10,-5,0,5,10) if dx or dy]
        areas=[weights(x+dx,y+dy) for dx,dy in offsets]
        for policy,allowed in [('strict',1<<21),('include_OUTLIER',(1<<19)|(1<<21))]:
            cont,res=contrast_image(st,ch,width,allowed)
            def net(im,ap,ring):
                valid=np.isfinite(im)
                if (ring&valid).sum()<12 or np.sum(ap*valid)<.97*np.sum(ap): return np.nan
                bg=np.nanmedian(im[ring])
                return np.nansum((im-bg)*ap)*MJYSR_TO_MJY
            flux=[net(res,ap,ring) for ap,_,ring in areas]
            sky=np.array(flux[1:]);sky=sky[np.isfinite(sky)]
            ap,core,ring=areas[0];oc=np.nanmedian(cont[ring]);orr=np.nanmedian(res[ring])
            # Centroid uses a 32-arcsec core. Clip opposite-sign light before
            # calculating its first moment, preserving sign for scalar contrasts.
            sign=np.sign(flux[0]);rcore=np.maximum(sign*np.nan_to_num(res-orr),0)*core
            ccore=np.maximum(np.nan_to_num(cont-oc),0)*core
            def centroid(a):
                if a.sum()<=0:return np.nan,np.nan
                return w.pixel_to_world_values(np.sum(a*pixx)/a.sum(),np.sum(a*pixy)/a.sum())
            rr,dd=centroid(rcore);rq,dq=centroid(ccore)
            finite_core=np.sum(core*np.isfinite(res))/np.sum(core)
            if rcore.sum()>0 and np.isfinite(rcore.sum()):
                iy,ix=np.unravel_index(np.argmax(rcore),rcore.shape)
                frac=rcore[iy,ix]/rcore.sum()
                shift=dist_arcsec(lon[iy,ix],lat[iy,ix],ra,dec)
                brightest_flag=int(st['FLAGS'][ch-1,iy,ix]);brightest_nhits=int(st['NHITS'][ch-1,iy,ix])
            else:
                frac=shift=np.nan;brightest_flag=brightest_nhits=''
            rows.append(dict(source_id=key,template_rank=t['template_rank'],peer_rank=t['peer_rank'],channel=ch,wave_um=t['peak_wave_um'],
                mask=policy,aperture19_signed_excess_mJy=flux[0],n_shifted=len(sky),
                shifted_median_mJy=np.median(sky) if len(sky) else np.nan,shifted_MAD_mJy=np.median(abs(sky-np.median(sky))) if len(sky) else np.nan,
                shifted_maxabs_mJy=np.max(abs(sky)) if len(sky) else np.nan,
                n_shifted_as_extreme=np.sum(abs(sky)>=abs(flux[0])) if np.isfinite(flux[0]) else '',
                usable_core32_fraction=finite_core,signed_excess_offset_arcsec=dist_arcsec(rr,dd,ra,dec),
                excess_to_cont_centroid_arcsec=dist_arcsec(rr,dd,rq,dq),continuum_offset_arcsec=dist_arcsec(rq,dq,ra,dec),
                single_pixel_fraction=frac,single_pixel_offset_arcsec=shift,single_pixel_flags=brightest_flag,single_pixel_nhits=brightest_nhits))
    jobs=[dict(label=k,ra_2025_9=sources[k]['extraction_ra'],dec_2025_9=sources[k]['extraction_dec']) for k in selected]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        unpack((ROOT/'localcache'/f'Npix{TILE}_cube.fits').read_bytes(),jobs,ROOT/'localcache',TILE,sink=record)
    write('spatial_localization.csv',rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('ranks',type=Path);a=p.parse_args();run(a.ranks)
