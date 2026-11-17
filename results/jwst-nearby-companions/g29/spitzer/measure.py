"""Astrometric and finite-aperture pixel checks on archival Spitzer images of G29-38."""
import csv
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
SOURCE = (352.193797, 5.246759)
HOST = (352.19670801911326, 5.247213099043008)  # Gaia DR3, reference epoch J2016.0
MOTION = np.array([-0.3982458354, -0.2667438097])  # east, north (arcsec/year)
JW_EPOCH = Time('2025-11-13T19:25:00').jyear
FWHM = {'I2': 1.72, 'I4': 1.98, 'M1': 5.9}  # arcsec, initial width for centroid fit

def offset(position, delta):
    ra, dec = position
    return ra + delta[0] / (3600 * np.cos(np.deg2rad(dec))), dec + delta[1] / 3600

def sample(image, xy):
    return map_coordinates(image, [xy[1], xy[0]], order=1, mode='constant', cval=np.nan)

def pixels(center, radius):
    cx, cy = center
    iy, ix = np.mgrid[int(np.floor(cy-radius)):int(np.ceil(cy+radius))+1,
                      int(np.floor(cx-radius)):int(np.ceil(cx+radius))+1]
    rel = np.array([ix-cx, iy-cy])
    return ix, iy, rel, np.hypot(*rel)

def centroid(img, pos, fw):
    ix, iy, rel, rad = pixels(pos, 1.7*fw)
    if not (ix.min()>=0 and iy.min()>=0 and ix.max()<img.shape[1] and iy.max()<img.shape[0]):
        return None
    vals = img[iy, ix]
    keep = np.isfinite(vals) & (rad <= 1.7*fw)
    if keep.sum() < 18:
        return None
    dx, dy = rel[:, keep]
    v = vals[keep]
    amp = max(v.max()-np.median(v), 1e-8)
    def residual(p):
        x, y, a, sig, bg, gx, gy = p
        return a*np.exp(-((dx-x)**2+(dy-y)**2)/(2*sig*sig)) + bg + gx*dx + gy*dy - v
    fit = least_squares(residual, [0,0,amp,fw/2.355,np.median(v),0,0],
                        bounds=([-fw,-fw,0,fw*.28,-np.inf,-np.inf,-np.inf],
                                [fw,fw,np.inf,fw*.8,np.inf,np.inf,np.inf]))
    noise = 1.4826*np.median(abs(fit.fun-np.median(fit.fun)))
    return np.array(pos)+fit.x[:2], fit.x[2]/max(noise,1e-10), fit.x[3]*2.355

def registration(img, wcs, epoch, fw, host):
    shifts = []
    for star in csv.DictReader((HERE/'gaia-field.csv').open()):
        sky = (float(star['ra']), float(star['dec']))
        pm = np.array([float(star['pmra']),float(star['pmdec'])])/1000
        position = np.array(wcs.world_to_pixel_values(*offset(sky, (epoch-2016)*pm)))
        if np.linalg.norm(position-host)<8*fw:
            continue
        fit = centroid(img, position, fw)
        if fit is not None and fit[1]>15 and np.linalg.norm(fit[0]-position)<.72*fw:
            shifts.append(fit[0]-position)
    if len(shifts)<4:
        return np.zeros(2), len(shifts), None
    shifts = np.array(shifts)
    med = np.median(shifts,axis=0)
    mad = np.maximum(1.4826*np.median(abs(shifts-med),axis=0),.10)
    shifts = shifts[np.all(abs(shifts-med)<3*mad,axis=1)]
    return np.median(shifts,axis=0), len(shifts), np.std(shifts,axis=0).tolist()

def finite_flux(img, host, pos, fw, ujypix):
    # Local-plane sky removal on the host-reflected image. Oversample the
    # circular 1.15-FWHM radius aperture for fractional edge pixel weights.
    ix,iy,rel,rad = pixels(pos,2.65*fw)
    vals=img[iy,ix] - sample(img, np.array([2*host[0]-ix, 2*host[1]-iy]))
    ann = (rad>=1.55*fw)&(rad<2.6*fw)&np.isfinite(vals)
    if ann.sum()<12:
        return np.nan
    dx,dy=rel
    basis=np.stack([np.ones_like(dx),dx,dy],axis=-1)
    for _ in range(2):
        coeff=np.linalg.lstsq(basis[ann],vals[ann],rcond=None)[0]
        resid=vals-basis@coeff
        mad=1.4826*np.median(np.abs(resid[ann]-np.median(resid[ann])))
        ann &= abs(resid)<max(3*mad,1e-12)
    weights=np.zeros_like(dx)
    subs=np.arange(-.4375,.5,.125)
    for offx in subs:
        for offy in subs:
            weights+=(np.hypot(dx+offx,dy+offy)<1.15*fw)/len(subs)**2
    if not np.all(np.isfinite(resid[weights>0])):
        return np.nan
    return np.sum(weights[weights>0]*resid[weights>0])*ujypix

def measure(path):
    img, h = fits.getdata(path,header=True)
    if h['BUNIT'] != 'MJy/sr':
        raise ValueError('Unexpected surface brightness units')
    w=WCS(h).celestial
    det=abs(np.linalg.det(w.pixel_scale_matrix))
    scale=det**.5*3600
    ujypix=det*(np.pi/180)**2*1e12
    date=Time(h['MJD_OBS'],format='mjd')
    epoch=date.jyear
    band=path.name.split('_')[1]
    fw=FWHM[band]/scale
    nominal=np.array(w.world_to_pixel_values(*offset(HOST,MOTION*(epoch-2016))))
    hostfit=centroid(img,nominal,fw)
    shift,n,scatter=registration(img,w,epoch,fw,nominal)
    host=hostfit[0]
    fixed=np.array(w.world_to_pixel_values(*SOURCE))+shift
    # Use the observed host center for the relative-position prediction.
    jw_host=offset(HOST,MOTION*(JW_EPOCH-2016))
    distance=np.array([(SOURCE[0]-jw_host[0])*3600*np.cos(np.deg2rad(SOURCE[1])),
                       (SOURCE[1]-jw_host[1])*3600])
    comoving=host+np.array(w.world_to_pixel_values(*offset(offset(HOST,MOTION*(epoch-2016)),distance)))-nominal
    found=centroid(img,fixed,fw)
    sky=None
    if found:
        sky=w.pixel_to_world_values(*(found[0]-shift))
    out=dict(filename=path.name,date=date.isot,band=band,scale_arcsec_pix=scale,gaia_ties=n,
             wcs_shift_pix=shift.tolist(),wcs_scatter_pix=scatter,host_xy_pix=host.tolist(),
             source_offset_fixed_vs_comoving_arcsec=np.linalg.norm(fixed-comoving)*scale,
             fixed_fit_ra_deg=None if sky is None else float(sky[0]),fixed_fit_dec_deg=None if sky is None else float(sky[1]),
             fixed_centroid_distance_arcsec=None if not found else np.linalg.norm(found[0]-fixed)*scale,
             fixed_peak_over_residual_mad=None if not found else found[1],
             fixed_fit_fwhm_arcsec=None if not found else found[2]*scale)
    for key,pos in [('fixed_sky',fixed),('comoving',comoving)]:
        flux=finite_flux(img,host,pos,fw,ujypix)
        dr=pos-host
        r=np.linalg.norm(dr)
        theta=np.arctan2(dr[1],dr[0])
        # Exclude fixed, moving and both reflected test apertures.
        angles=theta+2*np.pi*np.arange(max(8,int(2*np.pi*r/fw)))/max(8,int(2*np.pi*r/fw))
        controls=[]
        for a in angles:
            loc=host+r*np.array([np.cos(a),np.sin(a)])
            if min(np.linalg.norm(loc-p) for p in [fixed,comoving,2*host-fixed,2*host-comoving])<1.6*fw:
                continue
            controls.append(finite_flux(img,host,loc,fw,ujypix))
        ctr=np.array(controls);ctr=ctr[np.isfinite(ctr)]
        out[key]=dict(xy_pix=pos.tolist(),radius_arcsec=r*scale,flux_uJy=flux,aperture_radius_arcsec=1.15*fw*scale,
                      angular_control_median_uJy=np.median(ctr) if len(ctr) else None,
                      angular_control_mad_uJy=1.4826*np.median(abs(ctr-np.median(ctr))) if len(ctr)>3 else None,
                      angular_control_n=len(ctr))
    return out

if __name__=='__main__':
    outputs=[measure(p) for p in sorted((HERE/'cache').glob('*maic.fits'))]
    (HERE/'measurements.json').write_text(json.dumps(outputs,indent=2,allow_nan=False)+'\n')
    for x in outputs:
        print(x['band'],x['date'],x['gaia_ties'],x['fixed_centroid_distance_arcsec'],x['fixed_peak_over_residual_mad'],x['fixed_sky'],x['comoving'])
