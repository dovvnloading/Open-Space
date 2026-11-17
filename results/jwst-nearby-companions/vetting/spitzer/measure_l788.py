"""Forced empirical-PSF fits at fixed-sky and co-moving positions in the two IRAC2 mosaics."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates, maximum_filter
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
# Host sky position reported by MAST for the JWST observation; east is RA*cos(dec).
JW_RA,JW_DEC=335.781477,-17.612468
PM_E,PM_N=.308712,-.718387
EAST,NORTH=3.18,.22

def offset(ra,dec,e,n):
    return ra+e/(3600*np.cos(np.deg2rad(dec))),dec+n/3600

def sample(image,x,y):
    return map_coordinates(image,[y,x],order=1,mode='constant',cval=np.nan,prefilter=False)

def centroid(image,pred):
    # Match the saturated star to its 180-degree reflection using unsaturated wings.
    x0,y0=pred; y,x=np.mgrid[-18:19,-18:19]; x=x.flatten()+round(x0);y=y.flatten()+round(y0)
    a=image[y,x]; rad=np.hypot(x-x0,y-y0); choose=np.isfinite(a)&(rad>=5)&(rad<=14)&(a>0)
    a=a[choose];x=x[choose];y=y[choose]
    def loss(c):
        b=sample(image,2*c[0]-x,2*c[1]-y)
        return np.nanmedian(np.abs((a-b)/(np.maximum(a,0.2)+np.maximum(b,0.2))))
    grid=[(x0+dx,y0+dy) for dx in np.arange(-1,1.1,.25) for dy in np.arange(-1,1.1,.25)]
    best=min(grid,key=loss)
    return minimize(loss,best,method='Nelder-Mead',options={'initial_simplex':[best,(best[0]+.08,best[1]),(best[0],best[1]+.08)],'xatol':.005}).x

def psf(image,coverage,origin):
    # Use the brightest isolated unsaturated sources near the primary on the same mosaic.
    cx,cy=np.rint(origin).astype(int); w=235
    crop=image[cy-w:cy+w,cx-w:cx+w]; cov=coverage[cy-w:cy+w,cx-w:cx+w]
    filt=np.nan_to_num(crop,nan=-1e9)
    peaks=np.argwhere((filt==maximum_filter(filt,25))&(filt>.7)&(cov>1.5))
    dy,dx=np.mgrid[-18:19,-18:19]; rad=np.hypot(dx,dy)
    stamps=[];coords=[]
    for y,x in sorted(peaks,key=lambda p:filt[tuple(p)],reverse=True):
        x+=cx-w;y+=cy-w
        if np.hypot(x-cx,y-cy)<40: continue
        small=image[y-4:y+5,x-4:x+5];iy,ix=np.mgrid[-4:5,-4:5];wt=np.maximum(small-np.nanmedian(small[np.hypot(ix,iy)>3]),0)
        if not np.all(np.isfinite(wt)): continue
        sx=x+(ix*wt).sum()/wt.sum();sy=y+(iy*wt).sum()/wt.sum()
        stamp=sample(image,sx+dx,sy+dy); bg=np.nanmedian(stamp[(rad>13)&(rad<18)]);stamp-=bg
        flux=np.sum(stamp[rad<=10]);peak=stamp[18,18]
        if not np.isfinite(flux) or flux<=0 or peak/flux<.025 or peak/flux>.2: continue
        stamps.append(stamp/flux);coords.append([float(sx),float(sy)])
        if len(stamps)==8: break
    if len(stamps)<3: raise ValueError('Too few field PSF references')
    profile=np.nanmedian(stamps[1:],axis=0);profile/=np.sum(profile[rad<=10]);return profile,coords,stamps[0]

def flux(image,center,loc,template,scale,method,inject=None):
    # Fit only within 3 pixels of the requested locus. The antisymmetric residual
    # cancels a centrosymmetric halo while the constant + gradient accounts for local sky.
    px,py=loc; y,x=np.mgrid[round(py)-6:round(py)+7,round(px)-6:round(px)+7];dx=x-px;dy=y-py
    def at_pixels(stamp):
        p=sample(stamp,18+dx,18+dy)
        if method=='reflected_linear':
            p-=np.nan_to_num(sample(stamp,18+2*center[0]-x-px,18+2*center[1]-y-py))
        return p
    val=image[y,x].astype(float)
    if method=='reflected_linear':
        val-=sample(image,2*center[0]-x,2*center[1]-y)
    pred=at_pixels(template)
    if inject is not None:
        # Add a distinct held-out field-star PSF with 25 microJy inside the same 6-arcsec aperture.
        val+=25*at_pixels(inject)/scale
    fit_radius=3.7 if method=='reflected_linear' else 5.0
    choose=(np.hypot(dx,dy)<=fit_radius)&np.isfinite(val)&np.isfinite(pred);n=np.count_nonzero(choose)
    if n<19:return np.nan,n
    vals=val[choose];basis=np.column_stack([pred[choose],np.ones(n),dx[choose],dy[choose]])
    if method=='native_quadratic':
        basis=np.column_stack([basis,dx[choose]**2,dy[choose]**2,(dx*dy)[choose]])
    return np.linalg.lstsq(basis,vals,rcond=None)[0][0]*scale,n

def visit(path):
    image,h=fits.getdata(path,header=True);w=WCS(h).celestial
    aid=path.name.split('_')[2];cov=fits.getdata(next(path.parent.glob('*'+aid+'*mcov.fits')))
    scale=abs(np.linalg.det(w.pixel_scale_matrix))*(np.pi/180)**2*1e12   # microJy / pixel per MJy/sr
    dt=(Time(h['MJD_OBS'],format='mjd')-Time('2025-10-29T09:55:00')).to_value('yr')
    nominal=np.array(w.world_to_pixel_values(*offset(JW_RA,JW_DEC,PM_E*dt,PM_N*dt)),dtype=float)
    center=centroid(image,nominal);template,coords,injection=psf(image,cov,center)
    bg=np.array(w.world_to_pixel_values(*offset(JW_RA,JW_DEC,EAST,NORTH)))
    co=np.array(w.world_to_pixel_values(*offset(JW_RA,JW_DEC,EAST+PM_E*dt,NORTH+PM_N*dt)))+center-nominal
    output=[]
    for label,loc,method in [(k,l,m) for k,l in [('fixed_sky',bg),('co_moving_fixed_offset',co)] for m in ('reflected_linear','native_quadratic')]:
        f,n=flux(image,center,loc,template,scale,method);r=np.linalg.norm(loc-center)
        injected,_=flux(image,center,loc,template,scale,method,inject=injection)
        # Exclude the test position and its opposite; paired opposite controls are not independent.
        pa=np.arctan2(*(loc-center)[::-1]);num=max(6,int(np.floor(2*np.pi*r/(1.72/.6))))
        angles=pa+2*np.pi*np.arange(num)/num
        use=(abs(np.sin(angles-pa))*r>1.72/.6)
        controls=np.array([flux(image,center,center+r*np.array([np.cos(a),np.sin(a)]),template,scale,method)[0] for a in angles[use]])
        controls=controls[np.isfinite(controls)];sigma=1.4826*np.median(np.abs(controls-np.median(controls))) if len(controls)>3 else np.nan
        output.append({'aor':aid,'date':h['DATE_OBS'],'hypothesis':label,'method':method,'recovered_increment_for_25_uJy':round(injected-f,2),'x':round(loc[0],3),'y':round(loc[1],3),
            'host_x':round(center[0],3),'host_y':round(center[1],3),'radius_arcsec':round(r*.6,3),
            'nominal_to_fit_pixels':round(np.linalg.norm(nominal-center),3),'finite_pixels':n,'n_controls':len(controls),
            'forced_uJy_6as':round(f,2),'azimuth_mad_uJy':round(sigma,2),'control_min_uJy':round(min(controls,default=np.nan),2),
            'control_max_uJy':round(max(controls,default=np.nan),2),'control_fluxes_uJy':[round(t,1) for t in controls],
            'psf_ref_pixels':coords,'mjy_sr_to_ujy_pixel':round(scale,5)})
    return output

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);args=ap.parse_args()
    out=[]
    for file in sorted(args.cache.glob('*maic.fits')):out.extend(visit(file))
    if len(out)!=8:raise ValueError('expected four locus measurements for each of two fits')
    (HERE/'forced-photometry.json').write_text(json.dumps(out,indent=2,default=lambda v:v.item(),allow_nan=False)+'\n')
    for row in out:print({k:v for k,v in row.items() if k not in ['control_fluxes_uJy','psf_ref_pixels']})
