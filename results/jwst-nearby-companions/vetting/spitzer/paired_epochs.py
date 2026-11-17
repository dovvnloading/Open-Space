"""Subtract the other IRAC epoch as a reference PSF and test both predicted sky tracks."""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.optimize import least_squares

from measure_l788 import HERE, EAST, NORTH, JW_RA, JW_DEC, offset, psf, sample
from tie_wcs import ties


def mosaic(cache, aid, rows):
    image, header = fits.getdata(next(cache.glob('*'+aid+'*maic.fits')), header=True)
    cov = fits.getdata(next(cache.glob('*'+aid+'*mcov.fits')))
    row = next(r for r in rows if r['aor']==aid and r['hypothesis']=='fixed_sky')
    c = np.array([row['host_x'],row['host_y']]); w = WCS(header).celestial
    template, refs, held = psf(image,cov,c)
    delta, scatter, tie_rows = ties(image,w,Time(header['MJD_OBS'],format='mjd').jyear,template,c)
    co = next(r for r in rows if r['aor']==aid and r['hypothesis']=='co_moving_fixed_offset')
    return {'image':image,'center':c,'psf':template,'held':held,'wcs':w,'tie_rows':tie_rows,
            'tie_offset_pix':delta,'tie_scatter_pix':scatter,'fixed_sky':np.array([row['x'],row['y']])+delta,
            'co_moving_fixed_offset':np.array([co['x'],co['y']])}


def rotate(x,y,a):
    return x*np.cos(a)-y*np.sin(a),x*np.sin(a)+y*np.cos(a)


def subtracted(a,b):
    y,x=np.mgrid[-35:36,-35:36];xx,yy=x+np.rint(a['center'][0]).astype(int),y+np.rint(a['center'][1]).astype(int)
    dx=xx-a['center'][0];dy=yy-a['center'][1];v=a['image'][yy,xx];r=np.hypot(dx,dy)
    eligible=(r>=4)&(r<=28)&np.isfinite(v)&(v<45)
    # Exclude the two trial locations of a real source in both epochs while fitting the halo.
    for locus in ['fixed_sky','co_moving_fixed_offset']:
        X,Y=a[locus];eligible &= np.hypot(xx-X,yy-Y)>4.0
        X,Y=b[locus]-b['center'];eligible &= np.hypot(dx-X,dy-Y)>4.0
    def transform(p):
        x2,y2=rotate(dx,dy,p[2]);return x2+b['center'][0]+p[0],y2+b['center'][1]+p[1]
    # Validity mask accounts for bilinear footprints and avoids the saturated core in both mosaics.
    for sx in [-1.2,1.2]:
        for sy in [-1.2,1.2]:eligible &= np.isfinite(sample(b['image'],dx+b['center'][0]+sx,dy+b['center'][1]+sy))
    def difference(p):return v-p[3]*sample(b['image'],*transform(p))-p[4]
    def resid(p):return np.nan_to_num(difference(p)[eligible]/np.clip(abs(v[eligible])*.03,.08,2))
    # Mosaics have detector-oriented pixel axes. A small rigid registration absorbs the interpolation offset.
    fit=least_squares(resid,[0,0,0,1,0],bounds=([-.65,-.65,-.035,.82,-1],[.65,.65,.035,1.18,1]),loss='soft_l1')
    xr,yr=transform(fit.x)
    return xx,yy,xr,yr,difference(fit.x),fit.x


def forced(a,b,coordinates,loci,method,which='both',inject=False):
    xx,yy,xr,yr,res,p=coordinates; l0,l1=loci
    d0,e0=xx-l0[0],yy-l0[1];d1,e1=xr-l1[0],yr-l1[1]
    def template(key):
        return np.nan_to_num(sample(a[key],18+d0,18+e0))-p[3]*np.nan_to_num(sample(b[key],18+d1,18+e1))
    radius=3.7 if method=='plane' else 5.;m0=np.hypot(d0,e0)<radius;m1=np.hypot(d1,e1)<radius
    masks=[(m0,d0,e0),(m1,d1,e1)]
    if which=='2016':masks=masks[:1]
    if which=='2019':masks=masks[1:]
    keep=np.isfinite(res)&np.logical_or.reduce([t[0] for t in masks]); basis=[template('psf')]
    for region,x,y in masks:
        basis.extend([region,x*region,y*region])
        if method=='quadratic':basis.extend([x*x*region,x*y*region,y*y*region])
    design=np.column_stack([column[keep] for column in basis]);value=res[keep]
    unit=abs(np.linalg.det(a['wcs'].pixel_scale_matrix))*(np.pi/180)**2*1.e12
    if inject:value=value+25*template('held')[keep]/unit
    if len(value)<2*len(basis) or np.linalg.matrix_rank(design)<len(basis):return np.nan
    return np.linalg.lstsq(design,value,rcond=None)[0][0]*unit


def suite(a,b,coordinates,track,method):
    loci=[a[track],b[track]];centers=[a['center'],b['center']];rmin=min(np.linalg.norm(p-c) for p,c in zip(loci,centers))
    n=int(np.floor(2*np.pi*rmin/(1.72/.6)));angles=np.arange(n)*2*np.pi/n
    control=[]
    for ang in angles[1:]:
        test=[c+rotate(*(p-c),ang) for p,c in zip(loci,centers)]
        # Exclude regions overlapping either of the actual trial lobes in detector coordinates.
        if min(np.linalg.norm((t-c)-(l-d)) for t,c in zip(test,centers) for l,d in zip(loci,centers))<1.72/.6:continue
        val=forced(a,b,coordinates,test,method)
        if np.isfinite(val):control.append(val)
    val=forced(a,b,coordinates,loci,method); bg=np.median(control);sig=1.4826*np.median(np.abs(np.array(control)-bg));
    return {'hypothesis':track,'local_background':method,'signed_uJy':val,'control_median_uJy':bg,
            'control_MAD_sigma_uJy':sig,'n_controls':len(control),'controls_uJy':control,
            'control_centered_MAD_units':(val-bg)/sig,'lobe2016_uJy':forced(a,b,coordinates,loci,method,'2016'),
            'lobe2019_uJy':forced(a,b,coordinates,loci,method,'2019'),'25uJy_recovered_increment':forced(a,b,coordinates,loci,method,inject=True)-val}


def main(cache):
    rows=json.loads((HERE/'forced-photometry.json').read_text());a,b=[mosaic(cache,id,rows) for id in ('57886720','65885696')]; coords=subtracted(a,b)
    out={'registration':[float(x) for x in coords[-1]], 'gaia_ties':[{'median_shift_pix':q['tie_offset_pix'].tolist(),
        'rms_pix':q['tie_scatter_pix'].tolist(),'stars_used':sum(x[-1] for x in q['tie_rows'])} for q in [a,b]],
        'fits':[suite(a,b,coords,track,method) for track in ['fixed_sky','co_moving_fixed_offset'] for method in ['plane','quadratic']]}
    # How far does each individual residual peak lie from the predicted fixed position?
    # Search this 2.4-arcsec square only as a morphology diagnostic, without assigning a significance.
    grid=np.arange(-1.2,1.21,.15);m0,m1=a['fixed_sky'],b['fixed_sky'];peak=[]
    for which in ['2016','2019']:
        vals=[]
        for E in grid:
            for N in grid:
                test=[np.array(q['wcs'].world_to_pixel_values(*offset(JW_RA,JW_DEC,EAST+E,NORTH+N)))+q['tie_offset_pix'] for q in [a,b]]
                vals.append((forced(a,b,coords,test,'plane',which),E,N))
        f,e,n=max(vals,key=lambda t:t[0] if np.isfinite(t[0]) else -np.inf)
        peak.append({'year':which,'best_east_shift_arcsec':e,'best_north_shift_arcsec':n,'peak_flux_uJy':f})
    out['fixed_sky_peak_scan']=peak
    (HERE/'paired-results.json').write_text(json.dumps(out,indent=2,default=lambda x:x.item(),allow_nan=False)+'\n')
    for r in out['fits']: print(r['hypothesis'],r['local_background'],r['signed_uJy'],r['control_median_uJy'],r['control_MAD_sigma_uJy'],r['control_centered_MAD_units'],r['lobe2016_uJy'],r['lobe2019_uJy'])
    print(peak)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);main(ap.parse_args().cache)
