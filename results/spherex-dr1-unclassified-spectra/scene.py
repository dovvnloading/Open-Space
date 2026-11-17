"""Coadd spatial checks at each ranked peak, alternative masks/apertures and local continuum.

These checks reuse the same coadd. Local straight-line continua do not include spectral templates.
"""
import argparse
import collections
import csv
import math
import sys
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from acquire import ROOT, REPO, TILE
sys.path.insert(0,str(REPO/'results/spherex-waste-heat/spectra'))
from photometry import dist_arcsec, MJYSR_TO_MJY
from read_cubes import unpack
from extract import number, read, write

def run(rankpath):
    with rankpath.open() as f: ranks={t['source_id']:t for t in csv.DictReader(f)}
    with (ROOT/'sources.csv').open() as f:
        s={t['source_id']:t for t in csv.DictReader(f) if t['source_id'] in ranks}
    spec=collections.defaultdict(dict)
    with (ROOT/'spectra.csv').open() as f:
        for t in csv.DictReader(f):
            spec[t['source_id']][int(t['channel'])]=number(t['flux'])
    bright=read('bright'); br,bd=(np.array([float(t[k]) for t in bright]) for k in ('ra','dec'))
    catalog={t['designation']:t for t in read('phot')}
    out=[]
    def record(row,st):
        key=row['label']; src=s[key]; peak=ranks[key]
        ch=int(peak['peak_channel']); center=ch-1; width=int(peak['peak_fwhm_ch'])
        # Width-one: central channel vs offsets 2..4. Width-three: central
        # three channels vs offsets 4..6. All windows are inside one detector.
        span=4 if width==1 else 6; gap=span-2
        oi=np.r_[np.arange(center-span,center-gap+1),np.arange(center+gap,center+span+1)]
        ci=np.arange(center-width//2,center+width//2+1)
        if min(oi)<0 or max(oi)>=102 or any(st['DETECTOR'][oi]!=st['DETECTOR'][center]): raise ValueError('detector seam')
        wl=st['WAVELENGTH'];img=st['IMAGE'].astype(float); nh=st['NHITS'];fl=st['FLAGS']
        ny,nx=img.shape[1:]
        x,y=st['xy']
        w=WCS(fits.Header.fromstring(str(st['wcs']),sep=''));w.wcs.bounds_check(False,False)
        ra,dec=w.pixel_to_world_values(x,y)
        ox,oy=st['origin']
        bins=6
        xx,yy=np.meshgrid((np.arange(nx*bins)+.5)/bins+ox-.5,(np.arange(ny*bins)+.5)/bins+oy-.5)
        lon,lat=w.pixel_to_world_values(xx,yy)
        rad=dist_arcsec(lon,lat,ra,dec)
        xp,yp=np.meshgrid(np.arange(nx)+ox,np.arange(ny)+oy)
        def weight(lo,hi):
            return ((rad>=lo)&(rad<hi)).reshape(ny,bins,nx,bins).mean((1,3))
        dbright=dist_arcsec(br,bd,ra,dec)
        dbright=np.where(np.array([b['designation']!=key for b in bright]),dbright,np.inf)
        bi=int(np.argmin(dbright))
        band=[v for i,v in spec[key].items() if 1.2<=wl[i-1]<=1.3 and np.isfinite(v)]
        j=number(catalog[key]['j_m_2mass']);jflux=1594e3*10**(-.4*j)
        xymin=np.rint(st['xy']-st['origin']).astype(int)
        q0=nh[center,xymin[1],xymin[0]]
        ring=weight(36,60)>.5; h0=np.median(nh[center][ring])
        sizes=[(13,36,60),(19,36,60),(25,36,60),(19,60,90)]
        for maskname,allowed in [('strict',1<<21),('include_OUTLIER',(1<<21)|(1<<19))]:
            valid=np.isfinite(img)&(nh>0)&((fl.astype('u4')&np.uint32(0xffffffff^allowed))==0)
            for radius,rin,rout in sizes:
                ap=weight(0,radius);ring=weight(rin,rout)>.5
                flux=[];cloud=[];bg=[]
                for ii in [*oi,*ci]:
                    sky=np.nanmedian(np.where(ring&valid[ii],img[ii],np.nan))
                    support=np.sum(ap*valid[ii])/np.sum(ap)
                    ff=np.nansum(np.where(valid[ii],img[ii]-sky,0)*ap)*MJYSR_TO_MJY
                    if (ring&valid[ii]).sum()<12 or support<.97: ff=np.nan
                    flux.append(ff);bg.append(sky)
                    cloud.append(np.where(valid[ii],img[ii]-sky,np.nan))
                oflux,cflux=flux[:len(oi)],flux[len(oi):]
                valid_o=np.isfinite(oflux);valid_c=np.isfinite(cflux)
                contrast,continuum=np.nan,np.nan
                if sum(valid_o[:3])>=2 and sum(valid_o[3:])>=2 and np.all(valid_c):
                    coeff=np.polyfit(wl[oi][valid_o],np.asarray(oflux)[valid_o],1)
                    continuum=np.mean(np.polyval(coeff,wl[ci]))
                    contrast=np.mean(cflux)-continuum
                # Intensity centroid in the 25 arcsec core in the peak vs side
                # planes, only for positive net source signals. Not the
                # centroid of an intrinsically signed, possibly dipolar residual.
                def centroid(im):
                    z=np.maximum(np.where(np.isfinite(im),im,0),0)*weight(0,25)
                    if z.sum()<=0: return np.nan,np.nan
                    return w.pixel_to_world_values(np.sum(z*xp)/z.sum(),np.sum(z*yp)/z.sum())
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    contim=np.nanmedian(np.stack(cloud[:len(oi)]),axis=0)
                    peakim=np.nanmean(np.stack(cloud[len(oi):]),axis=0)
                    rc,dc=centroid(contim);rp,dp=centroid(peakim)
                out.append(dict(source_id=key,template_rank=int(peak['template_rank']),peer_rank=number(peak['peer_rank']),
                    peak_channel=ch,peak_wave_um=wl[center],width=width,score=number(peak['score']),mask=maskname,
                    radius_arcsec=radius,sky_inner=rin,sky_outer=rout,local_continuum_mJy=continuum,
                    peak_minus_continuum_mJy=contrast,fraction=contrast/continuum if continuum>0 else np.nan,
                    peak_aperture_mJy=np.mean(cflux),peak_center_nhits=q0,peak_ring_median_nhits=h0,
                    nhits_ratio=q0/h0 if h0>0 else np.nan,peak_center_flags=int(fl[center,xymin[1],xymin[0]]),
                    peak_centroid_offset_arcsec=dist_arcsec(rp,dp,ra,dec),
                    side_centroid_offset_arcsec=dist_arcsec(rc,dc,ra,dec),peak_side_shift_arcsec=dist_arcsec(rc,dc,rp,dp),
                    nearest_wise_arcsec=number(src['wise_sep_arcsec']),nearest_wise_W1=number(src['wise_neighbor_w1']),
                    n_other_wise19=src['wise_n19'],n_gaia19=src['gaia_n19'],closest_W1lt9_arcsec=dbright[bi],
                    closest_W1lt9_mag=number(bright[bi]['w1mpro']),tmass_jmag=j,j_proxy_ratio=np.median(band)/jflux if band and jflux>0 else np.nan))
    jobs=[dict(label=k,ra_2025_9=t['extraction_ra'],dec_2025_9=t['extraction_dec']) for k,t in s.items()]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        unpack((ROOT/'localcache'/f'Npix{TILE}_cube.fits').read_bytes(),jobs,ROOT/'localcache',TILE,sink=record)
    write('coadd_checks.csv',out)
if __name__=='__main__':
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('ranks',type=Path);args=p.parse_args();run(args.ranks)
