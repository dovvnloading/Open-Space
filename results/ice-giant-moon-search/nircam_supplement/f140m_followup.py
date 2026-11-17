"""Measure the two 2026 F140M short-frame sources and compare the known moons."""
import argparse
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from screen import ROOT, frame, write, geometry, ephemeris, prepare

EPOCH=2461321.533281779


def track(f, rad, deg):
    theta=np.deg2rad(deg)+geometry.frequency('Uranus',rad)*(f['t']-EPOCH)*86400
    u,v=rad*np.cos(theta),rad*np.sin(theta);m=f['m']
    den=(f['delta']+m[2,0]*u+m[2,1]*v)/geometry.ARCSEC
    return np.stack([(m[0,0]*u+m[0,1]*v)/den,(m[1,0]*u+m[1,1]*v)/den],axis=-1)@f['M']+f['cen']


def run(cache,eph):
    data=[];catalog=[];astrometry=[]
    for visit,ref in zip((6,8),eph):
        file=next(cache.glob(f'jw08975{visit:03d}*nrc*?_cal.fits'))
        if 'long' in file.name:
            file=next(p for p in cache.glob(f'jw08975{visit:03d}*_cal.fits') if 'long' not in p.name)
        f=frame(file,[ref],'Uranus');reftrack=track(f,76000,355)
        with fits.open(file,memmap=False) as h:
            im=h['SCI'].data;dq=h['DQ'].data;start,end=h[0].header['EXPSTART'],h[0].header['EXPEND']
            y,x=np.indices(im.shape);x0,y0=reftrack.mean(0)
            window=(abs(x-x0)<18)&(abs(y-y0)<18)&np.isfinite(im)&((dq&np.uint32(0xffffffff-4))==0)
            yg,xg=y[window]-y0,x[window]-x0;observed=im[window]
        xp,yp=(reftrack-[x0,y0]).T
        def model(p):
            dx,dy,amplitude,width,background,gx,gy=p
            shape=np.exp(-((xg[:,None]-xp-dx)**2+(yg[:,None]-yp-dy)**2)/(2*width**2)).mean(1)
            return amplitude*shape+background+gx*xg+gy*yg
        med=np.median(observed);noise=max(1e-9,1.4826*np.median(abs(observed-med)))
        initial=[0,0,max(observed.max()-med,1e-6),1.05,med,0,0]
        fit=least_squares(lambda p:(model(p)-observed)/noise,initial,loss='soft_l1',
                bounds=([-7,-7,0,.55,-np.inf,-np.inf,-np.inf],[7,7,np.inf,3,np.inf,np.inf,np.inf]))
        dx,dy,amp,width,*_=fit.x
        source=np.array([x0+dx,y0+dy]);r,dec=np.linalg.solve(f['M'].T,source-f['cen'])
        plane=np.linalg.solve(f['m'][:2],np.array([r,dec])*f['delta']/geometry.ARCSEC)
        phase=np.rad2deg(np.arctan2(plane[1],plane[0]))
        data.append(dict(filename=file.name,start_mjd=start,end_mjd=end,x=source[0],y=source[1],dx_from_trial=dx,dy_from_trial=dy,
                    psf_sigma_px=width,peak_amplitude=amp,local_mad=noise,amplitude_over_local_mad=amp/noise,
                    equatorial_deprojection_radius_km=np.linalg.norm(plane),phase_deg=phase))
        _,_,ids,obs=ephemeris.interpolate(ref,np.array([(start+end)/2]));b=obs[0,0]
        source_ra,source_dec=prepare.inverse_sky(*b[:2],r/geometry.ARCSEC,dec/geometry.ARCSEC)
        data[-1].update(ICRF_ra_deg=float(source_ra),ICRF_dec_deg=float(source_dec),relative_east_arcsec=r,relative_north_arcsec=dec,
                        fit_jacobian_x_se_px=float(np.sqrt(np.linalg.pinv(fit.jac.T@fit.jac)[0,0])),
                        fit_jacobian_y_se_px=float(np.sqrt(np.linalg.pinv(fit.jac.T@fit.jac)[1,1])))
        astrometry.append((f,b[:2],source))
        for k,body in enumerate(ids[1:],1):
            off=geometry.sky_xy(*obs[k,0,:2],*b[:2],1.,[0,0]);off[0]*=-1
            xy=off@f['M']+f['cen'];bodyplane=np.linalg.solve(f['m'][:2],off*f['delta']/geometry.ARCSEC)
            trueph=np.rad2deg(np.arctan2(bodyplane[1],bodyplane[0]))
            catalog.append(dict(filename=file.name,id=int(body),nominal_separation_px=np.linalg.norm(xy-source),
                         projected_equatorial_radius_km=np.linalg.norm(bodyplane),phase_offset_deg=(phase-trueph+180)%360-180))
    for i,row in enumerate(data):
        other=data[1-i];f,sky,source=astrometry[i]
        projected=geometry.sky_xy(other['ICRF_ra_deg'],other['ICRF_dec_deg'],*sky,1.,[0.,0.]);projected[0]*=-1
        frozen_sky=projected@f['M']+f['cen']
        frozen_planet=np.array([other['relative_east_arcsec'],other['relative_north_arcsec']])@f['M']+f['cen']
        row['miss_if_sky_fixed_in_other_visit_px']=float(np.linalg.norm(frozen_sky-source))
        row['miss_if_planet_fixed_in_other_visit_px']=float(np.linalg.norm(frozen_planet-source))
    write('f140m_sources.csv',data);write('f140m_catalog.csv',catalog)
    for r in data:print(r['filename'],r['equatorial_deprojection_radius_km'],r['phase_deg'],r['psf_sigma_px'],r['amplitude_over_local_mad'],
                    sorted((round(c['nominal_separation_px'],1),c['id'],round(c['phase_offset_deg'],2)) for c in catalog if c['filename']==r['filename'])[:5])
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('ephem',nargs=2,type=Path)
    a=p.parse_args();run(a.cache,a.ephem)
