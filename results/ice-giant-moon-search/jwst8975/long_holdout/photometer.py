"""F140M capsule-aperture photometry on four native 2026 Perdita trails."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from scipy.spatial import cKDTree
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time

root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root))
import geometry,ephemeris,search
from prepare import inverse_sky


def native_motion(h,ephem,center):
    """Propagate a near-circular path, using the measured orientation of the actual detector."""
    head=h[0].header;w=WCS(h['SCI'].header,naxis=2)
    dates=np.linspace(head['EXPSTART'],head['EXPEND'],65)
    _,_,ids,states=ephemeris.interpolate(ephem,dates)
    ra,dec,distance=states[0,:,0],states[0,:,1],states[0,:,4]*geometry.AU
    mjd=Time(dates,scale='utc',format='mjd').tdb.jd-distance/geometry.C/86400.
    ox=geometry.sky_xy(states[list(ids).index(725),:,0],states[list(ids).index(725),:,1],ra,dec,1.,[0,0])
    los,east,north=geometry.tangent_basis(ra,dec);pole=geometry.pole('Uranus',float(np.mean(mjd)))
    ray=los+(-ox[:,0,None]*east+ox[:,1,None]*north)/geometry.ARCSEC
    sat=(distance*(los@pole)/(ray@pole))[:,None]*ray-distance[:,None]*los
    angle=np.deg2rad(-19.8);rot=sat*np.cos(angle)+np.cross(pole,sat)*np.sin(angle)
    den=distance+(rot*los).sum(-1);e=(rot*east).sum(-1)/den;n=(rot*north).sum(-1)/den
    mid=len(dates)//2;cx,cy=center;rr,dd=w.all_pix2world(cx,cy,0)
    wr,wd=inverse_sky(rr,dd,e-e[mid],n-n[mid]);px,py=w.all_world2pix(wr,wd,0,quiet=True)
    # Local WCS Jacobian maps pixel offsets at the source to tangent-plane steradians.
    wr,wd=w.all_pix2world([cx,cx+1,cx],[cy,cy,cy+1],0)
    tangent=geometry.sky_xy(wr,wd,wr[0],wd[0],1.,[0,0])
    area=abs(np.linalg.det(tangent[1:]-tangent[0]))/geometry.ARCSEC**2
    return np.stack([px,py],axis=-1),area


def image_model(shape,track,scale,table):
    """An optical EE model averaged in native subpixels and over shutter-open time."""
    half=18;oversample=6
    edges=np.r_[0,table['aper_radius']];brightness=np.diff(np.r_[0,table['F140M']])/(np.pi*np.diff(edges**2))
    axis=(np.arange(-half,half+1)[:,None]+(np.arange(oversample)+.5)/oversample-.5).ravel()
    image=np.zeros(shape);weights=np.ones(len(track));weights[[0,-1]]=.5;weights/=weights.sum()
    for (x,y),weight in zip(track,weights,strict=True):
        left,top=int(np.floor(x)),int(np.floor(y))
        r=np.hypot(axis[:,None]-(y-top), axis[None,:]-(x-left))*scale
        psf=(brightness[np.searchsorted(edges[1:],r)]*(scale/oversample)**2).reshape(2*half+1,oversample,2*half+1,oversample).sum((1,3))
        placed=search.clipped_stamp(shape,top-half,left-half,weight*psf)
        if placed is not None:
            sl,s=placed;image[sl]+=s
    return image


def aperture(science,good,model,track,radius,order):
    xy=track[len(track)//2]
    x0,y0=np.maximum(np.floor(track.min(0)-15).astype(int),0)
    x1,y1=np.minimum(np.ceil(track.max(0)+16).astype(int),np.array(science.shape)[::-1])
    if x1<=x0 or y1<=y0:return None
    sy,sx=np.mgrid[y0:y1,x0:x1];tree=cKDTree(track)
    point=np.stack([sx,sy],axis=-1).reshape(-1,2)
    d=tree.query(point)[0];ap=np.zeros(len(point))
    for oy in (np.arange(4)+.5)/4-.5:
        for ox in (np.arange(4)+.5)/4-.5:
            ap+=tree.query(point+[ox,oy])[0]<=radius
    ap/=16;bg=(d>=7)&(d<=13)
    sl=np.s_[y0:y1,x0:x1];valid=good[sl].ravel();f=science[sl].ravel()
    support=(ap*valid).sum()/max(1.,ap.sum())
    if support<.8 or (bg&valid).sum()<70:return None
    weight=ap*valid;dx,dy=(point-xy).T
    basis=np.array([np.ones(len(dx)),dx,dy]+([dx**2,dx*dy,dy**2] if order==2 else [])).T
    # Fit a sigma-clipped background across both sides of the entire curved trail.
    def measurement(values):
        use=bg&valid
        for _ in range(3):
            par=np.linalg.lstsq(basis[use],values[use],rcond=None)[0]
            resid=values-basis@par
            s=1.4826*np.median(abs(resid[use]-np.median(resid[use])))
            if s<=0:break
            use=bg&valid&(abs(resid-np.median(resid[use]))<4*s)
        return np.nansum(weight*resid)
    value=measurement(f);unit=model[sl].ravel()
    # Keep native measurements and local sky estimates inside the throughput calibration.
    throughput=measurement(f+.50*unit)-measurement(f-.50*unit)
    if not .30<throughput<1.05:return None
    return value/throughput,throughput,support


def main(a):
    rows=[];table=np.genfromtxt(a.ee,names=True)
    inputs=[]
    locations=list(csv.DictReader(a.long_centers.open()))
    for r in locations:
        if r['accepted']=='True':inputs.append((a.long_files/(r['frame']+'.fits'),[float(r['native_x']),float(r['native_y'])],a.early_ephem))
    for r in csv.DictReader(a.short_centers.open()):
        v=int(r['filename'][7:10])
        inputs.append((a.short_files/r['filename'],[float(r['x']),float(r['y'])],a.early_ephem if v==6 else a.late_ephem))
    for path,center,ephem in inputs:
        with fits.open(path) as h:
            hdr=h['SCI'].header
            if hdr['BUNIT']!='MJy/sr':raise ValueError('native MJy/sr required')
            track,area=native_motion(h,ephem,center)
            science=h['SCI'].data.astype(float)*area*1e12
            good=np.isfinite(science)&np.isin(h['DQ'].data,[0,4])&np.isfinite(h['ERR'].data)&(h['ERR'].data>0)
            model=image_model(science.shape,track,geometry.ARCSEC*np.sqrt(area),table)
            length=np.linalg.norm(track[-1]-track[0])
            for radius in [2.5,3.5,5.0]:
                for order in [1,2]:
                    measured=aperture(science,good,model,track,radius,order)
                    if measured is None:continue
                    flux,through,support=measured
                    rows.append(dict(filename=path.name,visit=h[0].header['VISIT_ID'],detector=h[0].header['DETECTOR'],
                        source_sha256=search.file_hash(path),trail_length_px=length,radius_px=radius,background_order=order,
                        flux_uJy=flux,throughput=through,geometric_aperture_support=support,local_pixel_area_sr=area,
                        nominal_pixel_area_sr=hdr['PIXAR_SR'],effexptm_s=h[0].header['EFFEXPTM']))
    with (a.out/'photometry.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    summary={}
    for rad in (2.5,3.5,5.0):
        for order in (1,2):
            v=np.array([r['flux_uJy'] for r in rows if r['radius_px']==rad and r['background_order']==order])
            if len(v):
                summary[f'aperture_{rad}_sky_order_{order}']=dict(n=len(v),median_uJy=float(np.median(v)),
                    scaled_MAD_uJy=float(1.4826*np.median(abs(v-np.median(v)))),range_uJy=[float(v.min()),float(v.max())])
    primary=summary['aperture_3.5_sky_order_1'];flux=primary['median_uJy']
    out=dict(filter='F140M',primary=primary,AB_mag=float(-2.5*np.log10(flux*1e-6/3631)),alternates=summary,
        ee_sha256=search.file_hash(a.ee),EE_URL='https://jwst-docs.stsci.edu/files/216457060/216457083/1/1762453590926/Encircled_Energy_SW_ETCv2.txt',
        interpretation='Median across four exposures. Scaled MAD measures across-exposure scatter. Alternate apertures and backgrounds expose model sensitivity.')
    (a.out/'photometry.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['primary'],out['AB_mag'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['long_files','short_files','early_ephem','late_ephem','long_centers','short_centers','ee','out']:p.add_argument(k,type=Path)
    main(p.parse_args())
