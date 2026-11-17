"""Backward orbit ensembles in the Gala v2 MW plus a moving Hernquist LMC, in the MW frame."""
import argparse
import numpy as np
import astropy.units as u
from astropy.constants import G
from infer import GALAXY,draws,load_inputs
from make_queue import number,save
from screen_desi import cartesian

K=(u.km/u.s).to(u.kpc/u.Myr)
GM=G.to_value(u.kpc**3/u.Myr**2/u.Msun)*1.5e11
SCALE=17.14
FOCUS={'2752545725694798336','6340918677332110720','2294414169958935296','1717410430709663616','4444095982221534720'}

def lmc_accel(offset):
    r=np.linalg.norm(offset,axis=0);return -GM*offset/(np.maximum(r,1e-8)*(r+SCALE)**2)

def mw_accel(x):
    return GALAXY.acceleration(x*u.kpc).to_value(u.kpc/u.Myr**2)

def closest(a,b):
    delta=b-a; t=np.clip(-np.sum(a*delta,axis=0)/np.maximum(np.sum(delta**2,axis=0),1e-30),0,1)
    return np.linalg.norm(a+t*delta,axis=0)

def integrate(s,g,model,dt=.25,n=256,duration=300.):
    sample=draws(g,s,model,n)
    if sample is None:return None
    _,x,v=sample;v*=K
    # Independent observational draws of the LMC COM phase space. West PM sign converted to east.
    rng=np.random.default_rng(7913)
    L,W=cartesian(*[rng.normal(mu,err,n) for mu,err in [(78.76,.52),(-69.19,.25),(49.59,.54),(1.910,.020),(.229,.047),(262.2,3.4)]]);W*=K
    def accelerations(x,L):
        a=mw_accel(np.hstack([x,L])); reflex=lmc_accel(-L)
        return a[:,:n]+lmc_accel(x-L)-reflex,a[:,n:]-reflex
    a,b=accelerations(x,L);h=-dt;min_gc=np.linalg.norm(x,axis=0);min_lmc=np.linalg.norm(x-L,axis=0)
    disk_radius=np.full(n,np.nan);disk_speed=disk_radius.copy();disk_time=disk_radius.copy()
    for step in range(round(duration/dt)):
        xn=x+h*v+.5*h*h*a;Ln=L+h*W+.5*h*h*b;an,bn=accelerations(xn,Ln)
        vn=v+.5*h*(a+an); W+=.5*h*(b+bn)
        min_gc=np.minimum(min_gc,closest(x,xn)); min_lmc=np.minimum(min_lmc,closest(x-L,xn-Ln));
        c=(x[2]*xn[2]<=0)&np.isnan(disk_radius)
        if c.any():
            f=x[2,c]/(x[2,c]-xn[2,c]); q=x[:,c]+f*(xn[:,c]-x[:,c]);w=(v[:,c]+f*(vn[:,c]-v[:,c]))/K;R=np.hypot(q[0],q[1]);q[2]=0
            vc=GALAXY.circular_velocity(q*u.kpc).to_value(u.km/u.s)
            disk_radius[c]=R;disk_time[c]=(step+f)*dt
            disk_speed[c]=np.linalg.norm(w-np.vstack([vc*q[1]/R,-vc*q[0]/R,np.zeros_like(R)]),axis=0)
        x,v,L,a,b=xn,vn,Ln,an,bn
    out=dict(source_id=s['source_id'],model=model,draws=n,dt_myr=dt,lookback_myr=duration,
        p_gc_1kpc=np.mean(min_gc<1.),p_lmc_1kpc=np.mean(min_lmc<1.),p_lmc_5kpc=np.mean(min_lmc<5.),
        p_plane=np.mean(np.isfinite(disk_radius)),p_plane_R_lt25=np.mean(disk_radius<25))
    for k,col in [('gc_min_kpc',min_gc),('lmc_min_kpc',min_lmc),('first_plane_R_kpc',disk_radius),('first_plane_t_myr',disk_time),('disk_rel_kms',disk_speed)]:
        finite=col[np.isfinite(col)]
        for pct,val in zip((16,50,84),np.percentile(finite,[16,50,84]) if len(finite) else [np.nan]*3):out[f'{k}_{pct}']=round(float(val),4)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dt',type=float,default=.25);p.add_argument('--n',type=int,default=256);p.add_argument('--output',default='orbits.csv');p.add_argument('--only',nargs='*');a=p.parse_args()
    stars,gs=load_inputs();res=[]
    for s in stars:
        if a.only and s['source_id'] not in a.only:continue
        if s['source_id'] not in FOCUS and s['channel'] not in ('sdss','known_control','RRc_control'):continue
        for model in ('spec_plx','bj_photogeo'):
            r=integrate(s,gs[s['source_id']],model,a.dt,a.n)
            if r is not None:res.append(r)
    save(a.output,res)
