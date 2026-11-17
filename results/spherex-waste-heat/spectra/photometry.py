"""Aperture spectrophotometry and observer-frame feature contrasts on DR1 tiles."""
import argparse
import csv
from pathlib import Path
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT=Path(__file__).resolve().parent
# 512^2 equal-area pixels per order-5 tile, ICRS order 14 overall
MJYSR_TO_MJY = 1e9 * (4*np.pi / (12*16384**2))
FEATURES = {
    'PAH3.3': ((3.07,3.20),(3.23,3.36),(3.40,3.54)),
    'H2O3.0': ((2.60,2.78),(2.92,3.12),(3.42,3.57)),
    'CO2_4.27': ((4.13,4.21),(4.23,4.31),(4.34,4.43)),
    'CO4.67': ((4.40,4.54),(4.59,4.75),(4.82,4.98)),
    'PAH_E_z0.410': ((4.49,4.58),(4.61,4.66),(4.70,4.79)),
}

def dist_arcsec(lon,lat,ra,dec):
    x=np.deg2rad((lon-ra+180)%360-180)
    y=np.deg2rad(lat-dec)
    h=np.sin(y/2)**2 + np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(dec))*np.sin(x/2)**2
    return np.rad2deg(2*np.arcsin(np.sqrt(np.minimum(h,1))))*3600

def extract(stamp):
    img=stamp['IMAGE'].astype(float)
    nh,fl=stamp['NHITS'],stamp['FLAGS']
    w=WCS(fits.Header.fromstring(str(stamp['wcs']),sep=''))
    w.wcs.bounds_check(False,False)
    xy=stamp['xy']
    ox,oy=stamp['origin']
    ny,nx=img.shape[1:]
    # Subsample the celestial WCS for fractional aperture weights, including distorted HPX pixels.
    bins=6
    grid=(np.arange(bins)+.5)/bins-.5
    xx,yy=np.meshgrid(np.arange(nx*bins)/bins+grid[0]+ox,np.arange(ny*bins)/bins+grid[0]+oy)
    lon,lat=w.pixel_to_world_values(xx,yy)
    valid=(nh>0)&np.isfinite(img)&((fl==0)|(fl==(1<<21)))
    permissive=(nh>0)&np.isfinite(img)&((fl & np.uint32(~((1<<21)|(1<<19)) & 0xffffffff))==0)
    def aperture(x,y,mask):
        ra,dec=w.pixel_to_world_values(x,y)
        radii=dist_arcsec(lon,lat,ra,dec)
        def weights(inner,outer):
            return ((radii>=inner)&(radii<outer)).reshape(ny,bins,nx,bins).mean((1,3))
        ap=weights(0,19)
        ring=weights(36,60)
        # Reject cutoffs in the source aperture. Require enough independently located sky pixels.
        expected_area=np.pi*19**2/(MJYSR_TO_MJY/1e9*(180/np.pi*3600)**2)
        if ap.sum()<.95*expected_area or (ring>.5).sum()<15:
            return np.full(102,np.nan),np.zeros(102)
        bgmask=(ring>.5)&mask
        sky=np.nanmedian(np.where(bgmask,img,np.nan),axis=(1,2))
        area=(mask*ap).sum((1,2))
        frac=area/ap.sum()
        vals=np.nansum(np.where(mask,img-sky[:,None,None],0)*ap,axis=(1,2))*MJYSR_TO_MJY
        vals[(frac<.97)|(bgmask.sum((1,2))<12)]=np.nan
        return vals,frac
    f,cov=aperture(*xy,valid)
    loose,_=aperture(*xy,permissive)
    # Shift by identical integer pixel offsets to retain the target's subpixel phase.
    blanks=np.stack([aperture(xy[0]+dx,xy[1]+dy,permissive)[0]
                     for dx in [-12,-6,0,6,12] for dy in [-12,-6,0,6,12] if dx or dy])
    count=np.isfinite(blanks).sum(axis=0)
    center=np.nanmedian(blanks,axis=0)
    sigma=1.4826*np.nanmedian(np.abs(blanks-center),axis=0)
    sigma[count<8]=np.nan
    return f,loose,cov,blanks,sigma,count

def feature(y,wl,intervals):
    """Linear sideband interpolation in Fnu. Use the same center definition on sky apertures."""
    bands=[]
    for low,high in intervals:
        band=(wl>=low)&(wl<high)
        use=band & np.isfinite(y)
        if not use.any() or use.sum()<np.ceil(band.sum()/2): return np.nan,np.nan
        bands.append((np.mean(wl[use]),np.mean(y[use])))
    (l,fl),(m,fm),(r,fr)=bands
    base=fl+(fr-fl)*(m-l)/(r-l)
    return fm-base,base

def measure(files):
    spectral,feature_rows,summaries,indices=[],[],[],[]
    for path in sorted(files):
        label=path.stem
        with np.load(path) as stamp:
            wl=stamp['WAVELENGTH'].astype(float)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                f,loose,cov,blanks,s,n=extract(stamp)
            summaries.append(dict(label=label,available_ch=int(np.isfinite(loose).sum()),
                                  strict_ch=int(np.isfinite(f).sum()),usable_ch=int((np.isfinite(f)&np.isfinite(s)).sum())))
            for idx,l in enumerate(wl):
                spectral.append(dict(label=label,channel=idx+1,wavelength_um=l,
                                     aperture_mJy=f[idx],sky_sigma_mJy=s[idx],n_sky=n[idx],
                                     aperture_coverage=cov[idx],outlier_included_mJy=loose[idx]))
            for name,intervals in FEATURES.items():
                if name=='PAH_E_z0.410' and label.split('_')[0]!='E': continue
                for extraction,flux in [('strict',f),('include_OUTLIER',loose)]:
                    use=np.isfinite(flux)
                    val,base=feature(flux,wl,intervals)
                    sky_val=np.array([feature(np.where(use,b,np.nan),wl,intervals)[0] for b in blanks])
                    v=sky_val[np.isfinite(sky_val)]
                    zero,err=np.nan,np.nan
                    if len(v)>=8:
                        zero=np.median(v)
                        err=1.4826*np.median(abs(v-zero))
                    z=(val-zero)/err if err>0 else np.nan
                    feature_rows.append(dict(label=label,feature=name,extraction=extraction,contrast_mJy=val,
                                             continuum_mJy=base,blank_median_mJy=zero,sky_sigma_mJy=err,
                                             sky_z=z,n_blank=len(v)))
            for extraction,flux in [('strict',f),('include_OUTLIER',loose)]:
                band={}
                perturb={}
                for key,(lo,hi) in {'H':(1.51,1.72),'W1core':(3.16,3.62),'W2core':(4.42,4.90)}.items():
                    ind=(wl>=lo)&(wl<=hi)
                    valid=ind&np.isfinite(flux)
                    if valid.sum()<2 or valid.sum()<np.ceil(ind.sum()/2):
                        band[key]=np.nan
                        perturb[key]=np.full(blanks.shape[0],np.nan)
                    else:
                        band[key]=np.median(flux[valid])
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore',RuntimeWarning)
                            perturb[key]=np.nanmedian(flux[valid][None,:]+blanks[:,valid],axis=1)
                        perturb[key][np.isfinite(blanks[:,valid]).sum(axis=1)<np.ceil(valid.sum()/2)]=np.nan
                ratio=band['W2core']/band['H'] if band['H']>0 else np.nan
                with np.errstate(divide='ignore',invalid='ignore'):
                    ratios=perturb['W2core']/perturb['H']
                rgood=ratios[np.isfinite(ratios)&(perturb['H']>0)]
                sig=1.4826*np.median(np.abs(rgood-np.median(rgood))) if len(rgood)>=8 else np.nan
                indices.append(dict(label=label,extraction=extraction,H_mJy=band['H'],W1core_mJy=band['W1core'],
                                    W2core_mJy=band['W2core'],W2core_over_H=ratio,ratio_sky_sigma=sig,n_blank=len(rgood)))
    return spectral,feature_rows,summaries,indices

def save(name,rows):
    with (ROOT/name).open('w') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]); w.writeheader()
        for row in rows:
            w.writerow({k:('' if not np.isfinite(v) else f'{v:.7g}') if isinstance(v,(float,np.floating)) else v for k,v in row.items()})

def main():
    p=argparse.ArgumentParser()
    p.add_argument('cache',type=Path)
    a=p.parse_args()
    spectra,features,summaries,indices=measure(a.cache.glob('*.npz'))
    save('channel_spectra.csv',spectra); save('feature_tests.csv',features); save('coverage.csv',summaries); save('continuum.csv',indices)
if __name__=='__main__': main()
