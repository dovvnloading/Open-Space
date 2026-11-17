"""Fit near-circular, inclined first-order epicycles to 2025 only, then predict 2024."""
import argparse,csv,json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares
from astropy.time import Time
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import ARCSEC,AU,C,FIELD,basis,sky_xy,tangent_basis
from ephemeris import interpolate


def data(files,eph):
    rows=[]
    for file in files:
        with open(file) as f:rows.extend(csv.DictReader(f))
    t=np.array([float(r['mjd_utc']) for r in rows])
    _,_,_,ep=interpolate(eph,t);planet=ep[0]
    ra,dec=planet[:,0],planet[:,1];d=planet[:,4]*AU
    o=sky_xy(np.array([float(r['ra_deg']) for r in rows]),np.array([float(r['dec_deg']) for r in rows]),ra,dec,.031,(0.,0.))
    emit=Time(t,format='mjd',scale='utc').tdb.jd-d/C/86400
    return emit,d,*tangent_basis(ra,dec),o


def model(par,d,epoch):
    a,lam,k,h,qi,pi=par
    gm,r0,j2,j4=FIELD['Uranus'];q=(r0/a)**2;n2=gm/a**3
    az=np.sqrt(n2*(1+1.5*j2*q-1.875*j4*q*q))
    radial=np.sqrt(n2*(1-1.5*j2*q+5.625*j4*q*q))
    vertical=np.sqrt(2*az**2-radial**2)
    t,dist,los,east,north,_=d;dt=(t-epoch)*86400
    long=lam+az*dt
    A=long-(az-radial)*dt
    # k=e*cos(peri_long); h=e*sin(peri_long), defined at the emission epoch.
    radius=a*(1-k*np.cos(A)-h*np.sin(A))
    theta=long+2*az/radial*(k*np.sin(A)-h*np.cos(A))
    B=long-(az-vertical)*dt
    z=a*(qi*np.sin(B)-pi*np.cos(B))
    u,v,pole=basis('Uranus',epoch)
    pos=radius[:,None]*(np.cos(theta)[:,None]*u+np.sin(theta)[:,None]*v)+z[:,None]*pole
    denom=dist+(pos*los).sum(-1)
    return np.stack((-(pos*east).sum(-1)/denom,(pos*north).sum(-1)/denom),axis=1)*ARCSEC/.031


def run(feb,october,holdout,eph1,eph2,eph3,aliases,out):
    early=data([feb],eph1);late=data([october],eph2);hold=data(holdout,eph3)
    with open(Path(aliases)/'circular_aliases.csv') as f:trained=list(csv.DictReader(f))
    epoch=json.loads((Path(aliases)/'summary.json').read_text())['epoch2']
    rec=[]
    for row in trained:
        anchor=float(row['a_km']);angle=np.deg2rad(float(row['phase_late_deg']));cycle=int(row['cycle'])
        start=np.array([anchor,angle,0.,0.,0.,0.])
        def residual(p):
            return np.concatenate([(model(p,d,epoch)-d[-1]).ravel() for d in (early,late)])
        # Loose first-order sensitivity family with component caps, enclosing
        # the tabulated e=0.005 and i=1.6 deg. No priors or 2024 data enter.
        ext=np.array([25.,.1,.010,.010,np.deg2rad(2),np.deg2rad(2)])
        fit=least_squares(residual,start,bounds=(start-ext,start+ext),loss='soft_l1',f_scale=.6,
                          x_scale=[1.,.01,.002,.002,.01,.01],max_nfev=600)
        def rms(p,d):return float(np.sqrt(np.mean(np.sum((model(p,d,epoch)-d[-1])**2,axis=1))))
        rec.append(dict(cycle=cycle,model='linear eccentric/inclined epicycle',a_km=fit.x[0],
            phase_epoch_deg=np.rad2deg(fit.x[1])%360,e=np.hypot(*fit.x[2:4]),i_deg=np.rad2deg(np.hypot(*fit.x[4:])),
            k=fit.x[2],h=fit.x[3],i_cos_node=fit.x[4],i_sin_node=fit.x[5],
            feb_rms_px=rms(fit.x,early),oct_rms_px=rms(fit.x,late),training_cost=fit.cost,
            heldout_2024_rms_px=rms(fit.x,hold),circular_same_alias_2024_rms_px=rms(start,hold)))
    with open(out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rec[0]);w.writeheader();w.writerows(rec)
    print('best 2025 loss',min(rec,key=lambda r:r['training_cost']))
    print('min held-out residual',min(rec,key=lambda r:r['heldout_2024_rms_px']))


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for k in ['feb','oct','feb_ephem','oct_ephem','holdout_ephem','aliases','out']:a.add_argument(k)
    a.add_argument('--holdout',nargs='+',required=True)
    p=a.parse_args();run(p.feb,p.oct,p.holdout,p.feb_ephem,p.oct_ephem,p.holdout_ephem,p.aliases,p.out)
