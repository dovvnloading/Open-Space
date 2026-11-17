"""Fit the February counterpart on individual native RATEINTS planes and compare the orbital aliases."""
import argparse,csv,json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import warnings
import numpy as np
from scipy.optimize import least_squares,brentq
from astropy.time import Time
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'jwst9482'));sys.path.insert(0,str(Path(__file__).resolve().parent))
from link_visits import geometry,search,trials
from centroids import fit_patch
from independent_vet import save
from displaced import measurements,model
from prepare import inverse_sky
import positive_control as pc


def run(a):
    c=search.load_npz(a.cube);mapc=search.load_npz(a.maps)
    records=json.loads(Path(str(a.cube)+'.inputs.json').read_text())
    a.out.mkdir(parents=True,exist_ok=False)
    with TemporaryDirectory() as td:
        t=Path(td)
        b=trials(t/'e.csv',a.maps,t/'b.npz',np.array([[75511.75,0,0,0,0,254.55348]]),2460955.765061979)
    fitted=[]
    for j,rec in enumerate(records):
        m=pc.metadata(a.meta,rec['exposure_id'].removesuffix('_cal'))
        im,x0,y0=pc.wide_images(a.pixels,m)
        ra,dec=m['wcs'].all_pix2world(rec['center_native_x'],rec['center_native_y'],0)
        take=mapc['frame']==j
        coords=[]
        for k in [0,2,4]:
            sx=b['x_path'][0,take,k];sy=b['y_path'][0,take,k]
            sc=c['scale_arcsec']/geometry.ARCSEC
            wra,wdec=inverse_sky(ra,dec,-(sx-c['center_xy'][0])*sc,(sy-c['center_xy'][1])*sc)
            xx,yy=m['wcs'].all_world2pix(wra,wdec,0,quiet=True)
            coords.append(np.stack([xx,yy],axis=1))
        v=coords[2]-coords[0]
        patches=np.array(pc.patches_for_positions(im,x0,y0,{0:coords[1]})[0])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning);stack=np.nanmedian(patches,axis=0)
        yy,xx=np.mgrid[-pc.RADIUS:pc.RADIUS+1,-pc.RADIUS:pc.RADIUS+1]
        ok=(np.isfinite(patches).sum(axis=0)>=4)&(xx*xx+yy*yy<=6**2)
        if not np.isfinite(stack[ok]).any():continue
        ix=np.unravel_index(np.nanargmax(np.where(ok,stack,np.nan)),stack.shape)
        shift=np.array([xx[ix],yy[ix]])
        for k,g in enumerate(coords[1]+shift):
            p=fit_patch(im[k],x0,y0,g,*v[k])
            if p is not None:
                fitted.append(dict(j=j,k=k,t=m['times'][k],**p))
    save(a.out/'native_fits.csv',fitted)
    accepted=[r for r in fitted if r['edge_parameter']==0 and r['model_peak_contrast']>=3]
    cube,xy,pp,t,near,positions=measurements(a.cube,accepted,records,a.files)
    for i,r in enumerate(positions):
        r['fit_contrast']=accepted[i]['model_peak_contrast'];r['psf_sigma']=accepted[i]['psf_sigma']
    save(a.out/'astrometry.csv',positions)
    late=list(csv.DictReader(a.late_ast.open()));lc=search.load_npz(a.late_cube)
    lt=np.array([float(r['mjd_utc']) for r in late]);lxy=np.array([[float(r['x_sky']),float(r['y_sky'])] for r in late])
    lp=[]
    for r in late:
        j=np.flatnonzero(lc['exposure_id']==r['frame'])[0];tm=float(r['mjd_utc'])
        lp.append([np.interp(tm,lc['path_mjd'][j],k[j]) for k in [lc['path_ra'],lc['path_dec'],lc['path_delta']]])
    lp=np.array(lp)
    e2=2460955.765061979
    def geom(times,p,epoch):
        dt=(Time(times,format='mjd',scale='utc').tdb.jd-p[:,2]*geometry.AU/geometry.C/86400-epoch)*86400
        return dt,*geometry.tangent_basis(p[:,0],p[:,1])
    e1=float(np.median(Time(t,format='mjd',scale='utc').tdb.jd-pp[:,2]*geometry.AU/geometry.C/86400))
    def residual(par,cube,xy,p,tm,epoch):
        dt,los,east,north=geom(tm,p,epoch)
        pred=model(par[0],par[1],'Uranus',dt,cube,geometry.basis('Uranus',epoch),los,east,north,p[:,2])
        return (pred-xy).ravel()
    def fit1(cc,xy,p,tm,epoch,guess):
        fit=least_squares(residual,guess,args=(cc,xy,p,tm,epoch),x_scale=[1e4,1],loss='soft_l1',f_scale=.6)
        dr=np.linalg.norm(fit.fun.reshape(-1,2),axis=1)
        return fit.x,dict(a_km=fit.x[0],phase_deg=float(np.rad2deg(fit.x[1])%360),rms_px=np.sqrt(np.mean(dr**2)),median_px=np.median(dr))
    phase0=(np.deg2rad(254.55348)+geometry.frequency('Uranus',75511.75)*(e1-e2)*86400)%(2*np.pi)
    q1,ind1=fit1(cube,xy,pp,t,e1,[75511.75,phase0]);q2,ind2=fit1(lc,lxy,lp,lt,e2,[76337,np.deg2rad(254.746)])
    dt=(e2-e1)*86400
    pd=(q2[1]-q1[1])/(2*np.pi)
    spins=lambda r:geometry.frequency('Uranus',r)*dt/(2*np.pi)
    records2=[]
    for cycle in range(int(np.ceil(spins(78000)-pd)),int(np.floor(spins(74000)-pd))+1):
        r0=brentq(lambda r:spins(r)-cycle-pd,74000,78000)
        def both(par):return np.r_[residual(par,cube,xy,pp,t,e2),residual(par,lc,lxy,lp,lt,e2)]
        fit=least_squares(both,[r0,q2[1]],x_scale=[.2,.005],loss='soft_l1',f_scale=.6)
        rr=both(fit.x).reshape(-1,2);d1=np.linalg.norm(rr[:len(t)],axis=1);d2=np.linalg.norm(rr[len(t):],axis=1)
        records2.append(dict(cycle=cycle,a_km=fit.x[0],phase_late_deg=np.rad2deg(fit.x[1])%360,
            feb_rms=np.sqrt(np.mean(d1**2)),oct_rms=np.sqrt(np.mean(d2**2)),feb_median=np.median(d1),oct_median=np.median(d2),
            robust_cost=fit.cost))
    records2.sort(key=lambda z:z['robust_cost'])
    save(a.out/'circular_aliases.csv',records2)
    (a.out/'summary.json').write_text(json.dumps(dict(fits=len(fitted),accepted=len(accepted),late_n=len(late),epoch1=e1,epoch2=e2,
        independent_feb=ind1,independent_oct=ind2,best_circular=records2[0],second_circular=records2[1],
        median_feb_contrast=float(np.median([r['model_peak_contrast'] for r in accepted]))),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','maps','meta','pixels','files','late_cube','late_ast','out']:p.add_argument(k,type=Path)
    run(p.parse_args())
