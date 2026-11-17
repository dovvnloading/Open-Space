"""Reference-star subtraction, blind PSF search and source injection on G29-38 MIRI images."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage
from scipy.optimize import least_squares
from scipy.signal import correlate
from scipy.stats import t

from fetch import main as fetch_inputs

HERE = Path(__file__).resolve().parent
HALF = 98
YY, XX = np.mgrid[-HALF:HALF+1, -HALF:HALF+1]
R = np.hypot(XX, YY)
CUT = 12
PATCH = np.s_[HALF-CUT:HALF+CUT+1, HALF-CUT:HALF+CUT+1]
EDGES = (3, 4.5, 6, 7.5, 9)
FACTORS = (0.75, 1, 1.25, 1.5, 2, 3)
PEAK_KEYS = ('filter','east_arcsec','north_arcsec','ra','dec','snr','contrast')
def save_table(path, rows, keys=None):
    with path.open('w') as f:
        w=csv.DictWriter(f, fieldnames=keys or rows[0].keys())
        w.writeheader(); w.writerows(rows)

def prepare(path):
    with fits.open(path) as f:
        a=f['SCI'].data.astype(float); w=f['WHT'].data.astype(float)
        p=f[0].header; h=f['SCI'].header; z=WCS(h,fix=False)
        px,py=z.world_to_pixel_values(p['TARG_RA'],p['TARG_DEC'])
        # Restrict the host localization to the predicted archive position.
        x,y=int(round(float(px))),int(round(float(py)))
        gy,gx=np.mgrid[-8:9,-8:9]
        data=a[y-8:y+9,x-8:x+9]
        def cost(v):
            bg,d,xc,yc,s=v
            return (bg + d*np.exp(-((gx-xc)**2+(gy-yc)**2)/(2*s**2)) - data).ravel()
        bg=np.median(data); peak=data.max()-bg
        fit=least_squares(cost,[bg,peak,0,0,3],bounds=([-np.inf,0,-5,-5,1],[np.inf,np.inf,5,5,6]),
                          loss='soft_l1', f_scale=max(1,peak*.02)).x
        x+=fit[2];y+=fit[3]
        coords=(YY+y,XX+x)
        sky=np.nanmedian(a[w>.3*np.nanmax(w)])
        arr=ndimage.map_coordinates(np.nan_to_num(a-sky),coords,order=1)
        good=ndimage.map_coordinates((np.isfinite(a)&(w>.30*np.nanmax(w))).astype(float),coords,order=1)>.99
        mat=z.pixel_scale_matrix*3600
        return dict(image=arr,valid=good,x=x,y=y,wcs=z,mat=mat,scale=abs(np.linalg.det(mat))**.5,
                    solid=h['PIXAR_SR'],mjd=h['MJD-AVG'],coreweight=w[int(round(y)),int(round(x))]/np.nanmax(w),
                    gaussian_sigma=fit[-1],time=p['EFFEXPTM'])

def core(im):
    # Finite core/aperture convention held fixed for science and calibrators.
    sky=np.median(im[(R>=17)&(R<22)])
    data=(im[PATCH]-sky).copy()
    data[R[PATCH]>CUT]=0
    return data

def radial_reference(std):
    # Suppress off-axis stars and thermal stripes in the reference sight line.
    # Each calibration epoch contributes one independent, unit-core radial profile.
    ann=np.rint(R).astype(int)
    def profile(s):
        a=s['image']; norm=core(a).sum()
        radial=np.array([np.median(a[(ann==i)&s['valid']]) if np.any((ann==i)&s['valid']) else np.nan
                         for i in range(int(np.max(ann))+1)])
        valid=np.isfinite(radial)
        base=np.interp(np.arange(len(radial)),np.flatnonzero(valid),radial[valid])
        return ndimage.gaussian_filter1d(base,1)/norm
    smooth=np.mean([profile(s) for s in std],axis=0)
    return np.interp(R,np.arange(len(smooth)),smooth)

def matched(sci, ref, filt):
    fit_area=sci['valid']&(R<=14)
    image=sci['image']
    def cost(v):
        return image[fit_area] - v[0]*ref[fit_area] - v[1]
    v=least_squares(cost,[core(image).sum(),np.median(image[fit_area])],loss='soft_l1',f_scale=2).x
    residual=image-v[0]*ref-v[1]
    residual[~sci['valid']]=0
    return correlate(residual,filt,mode='same',method='fft')/np.sum(filt*filt)

def peak(data,x,y):
    def d(a,b,c):
        denom=a-2*b+c
        return float(np.clip((a-c)/(2*denom),-.5,.5)) if denom<0 else 0.
    return x+d(*data[y,x-1:x+2]), y+d(*data[y-1:y+2,x])

def main(cache,out):
    out.mkdir(exist_ok=True,parents=True)
    fetch_inputs(cache,False)
    rows=json.loads((HERE/'inputs.json').read_text())
    metadata=[]; curves=[]; trials=[]; detections=[]; checks=[]
    for band, fwhm_as in [('F2100W',.685),('F2550W',.83)]:
        rr=[r for r in rows if r['filters']==band]
        s=prepare(cache/rr[0]['filename'])
        std=[prepare(cache/r['filename']) for r in rr[1:]]
        # Matched kernel from the earlier calibration visit and injection PSF from later visit.
        # Calculate a unit-integral core with the same 12-pixel aperture photometry as G29-38.
        flux=core(s['image']).sum()
        kernel=core(std[0]['image']); kernel *= flux/kernel.sum()
        stamp=core(std[1]['image']); stamp *= flux/stamp.sum()
        inside=R[PATCH]<=CUT
        kernel[inside]-=np.mean(kernel[inside])
        ref=radial_reference(std)
        data=matched(s,ref,kernel)
        support=ndimage.minimum_filter(s['valid'].astype(int),2*CUT+5,mode='constant').astype(bool)
        local=data==ndimage.maximum_filter(data,size=7)
        scale=s['scale']; fwhm=fwhm_as/scale
        metadata.append(dict(filter=band,mjd=s['mjd'],x=s['x'],y=s['y'],scale_arcsec=scale,core_weight_fraction=s['coreweight'],
                             host_core_flux_mJy=flux*s['solid']*1e9,host_core_radius_arcsec=CUT*scale,
                             gaussian_sigma_pix=s['gaussian_sigma'],effexptm_s=s['time'],fit_ref_target='BD+60-1753'))
        rawpeaks=[]
        for lo,hi in zip(EDGES,EDGES[1:]):
            reg=(R*scale>=lo)&(R*scale<hi)
            valid=reg&support
            vals=data[valid]
            med=np.median(vals); sigma=1.4826*np.median(abs(vals-med))
            cov=valid.sum()/reg.sum()
            n=max(3,int(2*np.pi*lo*cov/fwhm_as)); fac=t.ppf(1-2.86651572e-7,n-1)*np.sqrt(1+1/n)
            boundary=med+fac*sigma
            for y,x in np.argwhere(valid&local&(data>=boundary)):
                rawpeaks.append(((data[y,x]-med)/sigma,x,y,data[y,x]-med))
            rad=(lo+hi)/2
            sites=[]
            for pa in np.arange(48)*360/48+2.4:
                ang=np.deg2rad(pa)
                delta=np.linalg.solve(s['mat'],rad*np.array([np.sin(ang),np.cos(ang)]))
                fx,fy=delta+HALF; x,y=np.rint([fx,fy]).astype(int)
                if not valid[y,x]:continue
                fake=ndimage.shift(np.pad(stamp,2),(fy-y,fx-x),order=3,mode='constant')
                response=correlate(fake,kernel,mode='full')/np.sum(kernel**2)
                k=len(response)//2
                sites.append((pa,x,y,fake,response[k-7:k+8,k-7:k+8]))
            if sites:
                pa,x,y,fake,patch=sites[0]
                changed=s['image'].copy()
                changed[y-CUT-2:y+CUT+3,x-CUT-2:x+CUT+3]+=fac*sigma*fake
                new=matched(dict(s,image=changed),ref,kernel)
                err=np.max(np.abs((new-data)[y-7:y+8,x-7:x+8]-fac*sigma*patch))
                if err>1e-8:raise ValueError((band,rad,err))
                checks.append(dict(filter=band,radius_arcsec=rad,max_full_vs_linear_error=err))
            yy,xx=np.mgrid[-7:8,-7:8]; near=np.hypot(xx,yy)<=fwhm/2
            def pass_at(a):
                return bool(np.any(near&(a==ndimage.maximum_filter(a,7))&(a>=boundary)))
            for factor in FACTORS:
                contrast=fac*sigma*factor
                rec=0
                for pa,x,y,_,patch in sites:
                    original=data[y-7:y+8,x-7:x+8]
                    ok=pass_at(original+contrast*patch) and not pass_at(original)
                    rec+=ok
                    trials.append(dict(filter=band,radius_arcsec=rad,pa_deg=pa,contrast=contrast,recovered=int(ok)))
                curves.append(dict(filter=band,inner_arcsec=lo,outer_arcsec=hi,radius_arcsec=rad,coverage_fraction=cov,
                                   independent_element_count=n,threshold_sigma=fac,threshold_contrast=fac*sigma,
                                   contrast=contrast,core_flux_uJy=contrast*flux*s['solid']*1e12,recovered=rec,tried=len(sites)))
        kept=[]
        for snr,x,y,c in sorted(rawpeaks,reverse=True):
            if any(np.hypot(x-xo,y-yo)<=fwhm for xo,yo in kept):continue
            kept.append((x,y)); xx,yy=peak(data,x,y)
            e,n=s['mat']@np.array([xx-HALF,yy-HALF])
            ra,dec=s['wcs'].pixel_to_world_values(s['x']+xx-HALF,s['y']+yy-HALF)
            detections.append(dict(filter=band,east_arcsec=e,north_arcsec=n,ra=float(ra),dec=float(dec),snr=snr,contrast=c))
    save_table(out/'visits.csv',metadata);save_table(out/'contrast_recovery.csv',curves)
    save_table(out/'injection_trials.csv',trials);save_table(out/'peaks.csv',detections,PEAK_KEYS)
    save_table(out/'linearity_check.csv',checks)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();main(a.cache,a.out)
