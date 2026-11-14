"""Joint timing fit using individual OGLE and Gaia magnitudes and TESS timing estimates."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from gaia_clock import to_hjd_utc
from orbit import get, P, EP, Y, TC, SCALE, circ
HERE = Path(__file__).resolve().parent
DATA = HERE / 'inputs'
HOURS = 24*P
KINDS = ['constant','linear','quadratic','cubic','kepler','sinusoid_offset','single_trough']

def load():
    h,t,ph,m,e=get()
    dg=(h>=2455000).astype(int)
    b=DATA/'gaia-epochs.csv'
    assert hashlib.sha256(b.read_bytes()).hexdigest() == json.loads((DATA/'gaia_receipt.json').read_text())['inputs'][b.name]['sha256']
    with b.open() as f:
        rows=[r for r in csv.DictReader(f) if r['variability_flag_g_reject']!='true']
    gh=to_hjd_utc(np.array([float(r['g_transit_time']) for r in rows])+2455197.5)
    gm=[float(r['g_transit_mag']) for r in rows]
    ge=np.hypot(1.0857362/np.array([float(r['g_transit_flux_over_error']) for r in rows]),.012)
    optical=dict(h=np.r_[h,gh],mag=np.r_[m,gm],error=np.r_[e,ge],survey=np.r_[dg,np.full(len(gh),2)])
    return optical, json.loads((HERE/'tess_measurements.json').read_text())

def kepler_displacements(t,U,e,angle):
    """Exact edge-on conjunction longitudes, expressed as departures from circular timing."""
    omega=angle+2*np.pi*t/U
    def shift(longitude):
        nu=(longitude-omega+np.pi)%(2*np.pi)-np.pi
        E=2*np.arctan2(np.sqrt(1-e)*np.sin(nu/2), np.sqrt(1+e)*np.cos(nu/2))
        return (E-e*np.sin(E)-nu)/(2*np.pi)
    return np.array([shift(np.pi/2),shift(3*np.pi/2)])

def geometry(kind,p,t):
    if kind=='kepler':
        g=kepler_displacements(t,*p)
        # Central finite-difference derivatives of analytic centers only.
        deriv=[]
        for i,eps in enumerate([.0001,1e-6,1e-5]):
            a=p.copy();b=p.copy()
            a[i]+=eps;b[i]-=eps
            deriv.append((kepler_displacements(t,*a)-kepler_displacements(t,*b))/(2*eps))
        return g,np.array(deriv)
    if kind=='single_trough':
        plateau,amplitude,t0,width=p
        q=(t-t0)/width;expo=np.exp(-.5*q*q)
        d=plateau-amplitude*expo
        mat=np.array([np.ones(len(t)),-expo,-amplitude*expo*q/width,-amplitude*expo*q*q/width])
    elif kind=='sinusoid_offset':
        U,offset,c,s=p
        theta=2*np.pi*t/U;co=np.cos(theta);sn=np.sin(theta)
        d=offset+c*co+s*sn
        mat=np.array([(c*sn-s*co)*theta/U,np.ones(len(t)),co,sn])
    else:
        mat=np.array([(t/10)**i for i in range(len(p))]);d=p@mat
    return np.array([-.5*d,.5*d]),np.stack([-.5*mat,.5*mat],axis=1)

def initial(kind,U=28):
    if kind=='single_trough':
        return [.055,.1,-3.,6.],[-.2,.0,-27.,.75],[.23,.45,22.,50.]
    if kind=='kepler':
        return [U,.065,3.13],[2.5,.00001,-20],[200,.3,20]
    if kind=='sinusoid_offset':
        return [U,0,-.04,0],[2.5,-.2,-.23,-.23],[200,.2,.23,.23]
    n=['constant','linear','quadratic','cubic'].index(kind)+1
    return [0]*n,[-.4]*n,[.4]*n

class JointFit:
    def __init__(self,optical,tess,kind,keep_o=None,keep_t=None,which='phase',sys_floor=0):
        self.kind=kind
        self.tess=tess
        if keep_o is None: keep_o=np.ones(len(optical['h']),dtype=bool)
        self.h=optical['h'][keep_o]; self.y=optical['mag'][keep_o]; self.err=optical['error'][keep_o]
        self.t=(self.h-TC)/Y; self.ph=(self.h-EP)/P%1
        orig=optical['survey'][keep_o]
        self.surveys,self.group=np.unique(orig,return_inverse=True)
        self.tkeys=[k for k in tess if keep_t is None or k in keep_t]
        times=[np.array(tess[k]['hjd_by_orbit']) for k in self.tkeys]
        self.bounds=np.cumsum([0]+[len(t) for t in times])
        self.times=np.r_[self.t,(np.concatenate(times)-TC)/Y] if times else self.t
        self.ty=np.array([tess[k][which] for k in self.tkeys])
        self.te=np.array([np.hypot(tess[k]['sigma'],sys_floor) for k in self.tkeys])
        p,lo,hi=initial(kind);self.ng=len(p)
        s=[]; l=[]; u=[]
        for g in self.surveys:
            base=18.15 if g==2 else 18.23
            s.extend([base,.42,.023,.41,.023])
            l.extend([base-.3,.015,.008,.015,.008]);u.extend([base+.3,1.,.08,1.,.08])
        self.start=np.r_[.015,0,p,s];self.lo=np.r_[-.15,-.12,lo,l];self.hi=np.r_[.15,.12,hi,u]
        self.cache=None
    def evaluate(self,p):
        if self.cache is not None and np.array_equal(self.cache,p):return
        self.cache=p.copy()
        n=len(self.h);ng=self.ng
        g,dg=geometry(self.kind,p[2:2+ng],self.times)
        c=p[0]+p[1]*self.t/10
        r1=circ(self.ph-c-g[0,:n]);r2=circ(self.ph-.5-c-g[1,:n])
        base,d1,w1,d2,w2=p[2+ng:].reshape(-1,5)[self.group].T
        f1=np.exp(-.5*(r1/w1)**2);f2=np.exp(-.5*(r2/w2)**2)
        c1=d1*f1*r1/w1**2;c2=d2*f2*r2/w2**2
        self.optical_pred=base+d1*f1+d2*f2
        j=np.zeros((n+len(self.ty),len(p)))
        j[:n,0]=c1+c2;j[:n,1]=(c1+c2)*self.t/10
        j[:n,2:2+ng]=(c1*dg[:,0,:n]+c2*dg[:,1,:n]).T
        row=np.arange(n); ind=2+ng+5*self.group
        for k,v in enumerate([1.,f1,c1*r1/w1,f2,c2*r2/w2]): j[row,ind+k]=v
        j[:n]/=self.err[:,None]
        ts=.5+g[1,n:]-g[0,n:];tj=dg[:,1,n:]-dg[:,0,n:]
        self.tess_pred=[]
        for i,(a,b) in enumerate(zip(self.bounds[:-1],self.bounds[1:])):
            self.tess_pred.append(ts[a:b].mean())
            j[n+i,2:2+ng]=tj[:,a:b].mean(axis=1)/self.te[i]
        self.resid=np.r_[(self.optical_pred-self.y)/self.err,(self.tess_pred-self.ty)/self.te]
        self.jac=j
    def fun(self,p): self.evaluate(p);return self.resid
    def derivative(self,p): self.evaluate(p);return self.jac
    def rho(self,z):
        # Robust soft-L1 for magnitudes; quadratic loss for bootstrap-standardized TESS spacings.
        a=np.sqrt(1+z)
        out=np.array([2*(a-1),1/a,-.5/a**3])
        n=len(self.y)
        out[:,n:]=[z[n:],np.ones_like(z[n:]),np.zeros_like(z[n:])]
        return out
    def solve(self,start=None):
        sol=least_squares(self.fun,self.start if start is None else start,jac=self.derivative,
                          bounds=(self.lo,self.hi),loss=self.rho,f_scale=SCALE,x_scale='jac',
                          ftol=1e-9,gtol=1e-8,max_nfev=220)
        self.evaluate(sol.x)
        return sol
    def multi(self,grid):
        if self.kind=='single_trough':
            fits=[]
            for center in [-10.,-2.,5.]:
                for width in [2.,6.,14.]:
                    p=self.start.copy();p[4:6]=center,width
                    fits.append(self.solve(p))
            return min(fits,key=lambda sol:sol.cost)
        if self.kind not in ['kepler','sinusoid_offset']: return self.solve()
        fits=[]
        for U in grid:
            p=self.start.copy();p[2]=U
            fits.append(self.solve(p))
        return min(fits,key=lambda sol:sol.cost)

def timed_prediction(kind,pars,times):
    t=(np.atleast_1d(times)-TC)/Y
    g,_=geometry(kind,np.asarray(pars)[2:2+len(initial(kind)[0])],t)
    return .5+g[1]-g[0]

def stats(fit,sol,tess):
    p=sol.x
    gap={}
    for k,d in tess.items():
        pred=timed_prediction(fit.kind,p,d['hjd_by_orbit']).mean()
        gap[k]=dict(pred_h=HOURS*pred,observed_h=d['phase']*HOURS,resid_min=(d['phase']-pred)*HOURS*60,
                    standardized_residual=(d['phase']-pred)/d['sigma'])
    out=dict(loss=2*sol.cost,n_params=len(p),parameters=p.tolist(),tess=gap,
             predictions_hours={date:float(HOURS*timed_prediction(fit.kind,p,jd)[0])
                                for date,jd in [('2027-04-16',2461511.5),('2027-07-01',2461587.5),
                                                ('2028-07-01',2461953.5),('2030-07-01',2462683.5)]})
    if fit.kind in ['kepler','sinusoid_offset']:out['period_years']=p[2]
    if fit.kind=='kepler':out['edge_on_eccentricity']=p[3]
    return out

def study():
    optical,tess=load();s=optical['survey'];h=optical['h']
    full=np.ones(len(h),bool);tk=list(tess)
    schemes={
      'full':(full,tk,'phase',0.),
      'wider_shapes':(full,tk,'wide_phase',0.),
      'alternate_extraction':(full,tk,'alternate',0.),
      'extra_0.008_cycle_tess_error':(full,tk,'phase',.008),
      'ogle_i_only':(s!=2,[],'phase',0.),
      'all_optical_only':(full,[],'phase',0.),
      'before_2025':(full,tk[:-1],'phase',0.),
      'before_2022':(h<2459580.5,tk[:2],'phase',0.),
    }
    result={}
    for name,(sel,ts,key,floor) in schemes.items():
        result[name]={}
        for kind in KINDS:
            f=JointFit(optical,tess,kind,sel,ts,key,floor)
            grid=np.geomspace(3,160,24) if name=='full' else [8,17,28,65,140]
            sol=f.multi(grid)
            row=stats(f,sol,tess)
            if name=='before_2022':
                use=(h>=2459580.5)&(s!=2)
                # Extract shape parameters for the trained OGLE-IV survey only.
                g,_=geometry(kind,sol.x[2:2+f.ng],(h[use]-TC)/Y)
                c=sol.x[0]+sol.x[1]*(h[use]-TC)/Y/10
                phase=(h[use]-EP)/P%1
                ix=list(f.surveys).index(1)
                b,d1,w1,d2,w2=sol.x[2+f.ng:].reshape(-1,5)[ix]
                pred=b+d1*np.exp(-.5*(circ(phase-c-g[0])/w1)**2)+d2*np.exp(-.5*(circ(phase-.5-c-g[1])/w2)**2)
                r=(pred-optical['mag'][use])/optical['error'][use]
                row['held_late_ogle_loss']=float(sum(2*SCALE**2*(np.sqrt(1+(r/SCALE)**2)-1)))
                row['held_late_ogle_n']=int(sum(use))
            result[name][kind]=row
        print(name,{k:[round(v['loss'],1),round(v.get('period_years',0),2)] for k,v in result[name].items()},flush=True)
    (HERE/'fits.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': study()
