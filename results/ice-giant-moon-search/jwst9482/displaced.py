"""Fit astrometry of the compact feature displaced from the screened orbit."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy import optimize
from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS
from link_visits import geometry,search
from independent_vet import save
ROOT=Path(__file__).resolve().parent


def measurements(cube,records,cent,files):
    c=search.load_npz(cube);n=len(c['exposure_id']);xy=[];pp=[];times=[];nearest=[];positions=[]
    for row in records:
        j=int(row['j']);t=float(row['t']);xp,yp=float(row['x']),float(row['y'])
        fid=str(c['exposure_id'][j]);hd=fits.getheader(files/(fid+'.fits'),'SCI');w=WCS(hd,naxis=2);ref=cent[j]
        center=w.all_pix2world(ref['center_native_x'],ref['center_native_y'],0)
        sky=w.all_pix2world(xp,yp,0)
        x,y=geometry.sky_xy(*sky,*center,float(c['scale_arcsec']),c['center_xy'])
        pos=np.array([np.interp(t,c['path_mjd'][j],q[j]) for q in (c['path_ra'],c['path_dec'],c['path_delta'])])
        # Use the actual telescope-relative planet position to convert offsets to ICRF.
        los,east,north=geometry.tangent_basis(pos[0],pos[1])
        native=geometry.direction(*sky);l,e,north0=geometry.tangent_basis(*center)
        vec=los+east*(native@e)/(native@l)+north*(native@north0)/(native@l);vec/=np.linalg.norm(vec)
        sky_ra=np.rad2deg(np.arctan2(vec[1],vec[0]))%360;sky_dec=np.rad2deg(np.arcsin(vec[2]))
        positions.append(dict(frame=fid,integration=int(row['k'])+1,mjd_utc=t,ra_deg=round(sky_ra,9),dec_deg=round(sky_dec,9),x_sky=round(x,4),y_sky=round(y,4)))
        known=np.array([[np.interp(t,c['path_mjd'][j],v[:,k]) for k in range(2)] for v in c['known_path'][:,j]])
        near=np.linalg.norm(known-[x,y],axis=1);nearest.append(near)
        xy.append([x,y]);pp.append(pos);times.append(t)
    return c,np.array(xy),np.array(pp),np.array(times),np.array(nearest),positions


def model(a,phase,planet,dt,c,frame,los,east,north,delta):
    q=phase+geometry.frequency(planet,a)*dt
    xyz=a*(np.cos(q)[:,None]*frame[0]+np.sin(q)[:,None]*frame[1])
    den=delta*geometry.AU+np.sum(xyz*los,axis=1)
    return c['center_xy']+np.column_stack([-np.sum(xyz*east,axis=1),np.sum(xyz*north,axis=1)])*(geometry.ARCSEC/float(c['scale_arcsec']))/den[:,None]


def run(a):
    rows=list(csv.DictReader(a.centroids.open()));det=json.loads((ROOT/'inputs.json').read_text())
    c,xy,pp,t,near,pos=measurements(a.cube,rows,det,a.files)
    cfg=json.loads((ROOT/'independent_vet.json').read_text())['F150W2'];epoch=cfg['emission_epoch_tdb_jd']
    dt=(Time(t,format='mjd',scale='utc').tdb.jd-pp[:,2]*geometry.AU/geometry.C/86400-epoch)*86400
    los,east,north=geometry.tangent_basis(pp[:,0],pp[:,1]);frame=geometry.basis('Uranus',epoch)
    def residual(par):
        return (model(par[0],par[1],'Uranus',dt,c,frame,los,east,north,pp[:,2])-xy).ravel()
    best=optimize.least_squares(residual,[cfg['radius_km'],np.deg2rad(cfg['phase_deg'])],x_scale=[1e4,1],loss='soft_l1',f_scale=1.)
    offsets=residual(best.x).reshape(-1,2);dr=np.linalg.norm(offsets,axis=1)
    report=dict(epoch_tdb_jd=epoch,radius_km=float(best.x[0]),phase_deg=float(np.rad2deg(best.x[1])%360),
                median_fit_px=float(np.median(dr)),rms_fit_px=float(np.sqrt((dr**2).mean())),
                nearest_known=[dict(id=int(c['known_id'][k]),median_pix=float(np.median(near[:,k])),minimum_pix=float(np.min(near[:,k]))) for k in np.argsort(np.median(near,axis=0))[:6]])
    for r,dx,dy,ds in zip(pos,offsets[:,0],offsets[:,1],near):
        r['model_minus_x_px']=round(dx,3);r['model_minus_y_px']=round(dy,3)
        k=ds.argmin();r['known_id']=int(c['known_id'][k]);r['known_sep_px']=round(ds[k],3)
    save(ROOT/'integration_check/displaced_astrometry.csv',pos)
    (ROOT/'integration_check/displaced_circular.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','files','centroids']:ap.add_argument(k,type=Path)
    run(ap.parse_args())
