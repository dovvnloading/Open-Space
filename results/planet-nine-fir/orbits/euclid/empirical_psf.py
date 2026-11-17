"""Estimate the registered VIS PSF from Gaia stars in the measurements themselves."""
import csv
import io
import json
from warnings import catch_warnings,simplefilter
import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates
from dense import BASE,HERE,CACHE,design
from mosaic import xy

R=11

def main():
    stars=list(csv.DictReader((HERE/'psf_gaia.response').open()))
    geo=json.loads((BASE/'geometry.json').read_text())
    names,ts,_,_=design()
    dx,dy=np.meshgrid(np.arange(-R,R+1),np.arange(-R,R+1))
    dist=np.hypot(dx,dy)
    norms=[]
    parts=[]
    details={}
    for name,t in zip(names,ts):
        d=np.load(BASE/(name+'_data.npy'),mmap_mode='r')
        w=np.load(BASE/(name+'_weight.npy'),mmap_mode='r')
        psfs=[]
        for s in stars:
            # Gaia DR3 positions refer to 2016.0. Over nine years high PM can matter.
            ra,de=float(s['ra']),float(s['dec'])
            years=(t-57388.5)/365.25
            ra+=float(s['pmra'] or 0)*years/3.6e6/np.cos(np.radians(de))
            de+=float(s['pmdec'] or 0)*years/3.6e6
            x,y=xy(ra,de);x-=geo['x0'];y-=geo['y0']
            if not(R+5<x<d.shape[1]-R-6 and R+5<y<d.shape[0]-R-6): continue
            ix,iy=int(round(x)),int(round(y))
            local=d[iy-3:iy+4,ix-3:ix+4]
            good=w[iy-3:iy+4,ix-3:ix+4]>0
            if not good[2:5,2:5].all() or good.mean()<.85:continue
            jy,jx=np.unravel_index(local.argmax(),local.shape)
            if np.hypot(jx-3,jy-3)>2:continue
            for _ in range(2):
                z=map_coordinates(d,[y+dy,x+dx],order=1,prefilter=False)
                good=map_coordinates((w[iy-R-3:iy+R+4,ix-R-3:ix+R+4]>0).astype(float),
                    [y-iy+R+3+dy,x-ix+R+3+dx],order=1,prefilter=False)>.999
                bg=np.median(z[(dist>8)&good]) if np.any((dist>8)&good) else 0
                zz=(z-bg)*good*(dist<=3)
                if zz.sum()<=0:break
                x+=np.sum(dx*zz)/zz.sum(); y+=np.sum(dy*zz)/zz.sum()
            if (z[dist<1].max()-bg)<.3 or not good[dist<1.5].all():continue
            if np.sum((z-bg)[(dist<3)&good])<=0:continue
            z=np.where(good,z-bg,np.nan)
            flux=np.nansum(z[dist<=6.5])
            if flux<.5:continue
            psfs.append(z/flux)
        parts.append(psfs)
        for fp in sorted(CACHE.glob(name+'_*SCI_hdr.bin')):
            h=fits.Header.fromfile(io.BytesIO(fp.read_bytes()))
            if h.get('APERCORR') and np.isfinite(float(h['APERCORR'])): norms.append(float(h['APERCORR']))
    pool=[s for plist in parts for s in plist]
    if len(pool)<3: raise ValueError('Insufficient unflagged calibrating star measurements')
    apcorr=float(np.median(norms))
    f13=10**(.4*apcorr)
    profiles=[]
    for name,ps in zip(names,parts):
        samples=ps if len(ps)>=3 else pool
        with catch_warnings():
            simplefilter('ignore',RuntimeWarning)
            profile=np.nanmedian(samples,axis=0)
        profile=np.maximum(0,np.nan_to_num(profile))
        profile-=np.median(profile[dist>9])
        profile=np.maximum(0,profile)
        profile*=f13/profile[dist<=6.5].sum()
        profiles.append(profile)
        details[name]=dict(local_stars=len(ps),used=len(samples),sum_stamp=float(profile.sum()),
                           core_fwhm_gaussian_equivalent_as=float(.1*2.35482*np.sqrt(
                               np.sum(profile[dist<=2.5]*dist[dist<=2.5]**2)/(2*np.sum(profile[dist<=2.5])))))
    np.save(CACHE/'measured_psf.npy',np.array(profiles).astype('f4'))
    (HERE/'measured_psf.json').write_text(json.dumps(dict(median_APERCORR=apcorr,
         enclosed_flux_13pix=f13,exposures=details),indent=2)+'\n')
    print(details)


if __name__=='__main__': main()
