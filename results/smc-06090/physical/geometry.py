"""Fit precessing, spherical, mutually occulting stars to both OGLE I-band series."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
P, EPOCH, TREF = 1.1696719, 2460000.0369, 2456000.
HASH = ('d57614117ef140adf1413baab42f92182b9c39752eed39706b6fc040e203a897',
        '429caaf74f2086f3a8f53c71a91bfb3c07049dbc2e81122c06646387c8fa21a0')
KEYS = ('phase_offset','phase_drift_per_decade','apsidal_years','eccentricity','omega_ref_rad',
        'sum_fractional_radii','radius_ratio_b_a','cos_inclination','intensity_ratio_b_a',
        'baseline_iii','baseline_iv','cos2_amplitude_mag','cos1_amplitude_mag')
START = np.array([.0153,-.0031,28.3,.0678,-2.6632,.35,.85,.075,.95,18.23,18.23,0.,0.])
LO = np.array([-.12,-.1,6.,.00001,-7.,.14,.32,.00001,.35,18.0,18.0,-.15,-.15])
HI = np.array([.12,.1,145.,.32,7.,.62,2.6,.35,1.8,18.5,18.5,.15,.15])
SCALE=1.8

def load(folder):
    chunks=[]
    for g, label in enumerate(('ogle3.dat','ogle4.dat')):
        p=folder/label
        assert sha256(p.read_bytes()).hexdigest()==HASH[g]
        a=np.loadtxt(p)[:,:3]
        a=a[np.isfinite(a).all(axis=1)&(a[:,2]>0)&(a[:,2]<=.1)]
        if g==0: a[:,0]+=2450000.
        chunks.append(np.c_[a,np.full(len(a),g)])
    h,mag,err,survey=np.concatenate(chunks).T
    return h,mag,np.hypot(err,.012),survey.astype(int)

def overlap(r,k,d):
    """Area / pi shared by circles of radii r and k. Supports broadcast."""
    rr,kk,dd=np.broadcast_arrays(r,k,d)
    out=np.minimum(rr,kk)**2
    out=np.where(dd>=rr+kk,0.,out)
    cross=(dd>abs(rr-kk))&(dd<rr+kk)
    a,b,s=rr[cross],kk[cross],dd[cross]
    q=np.maximum((a+b+s)*(-a+b+s)*(a-b+s)*(a+b-s),0.)
    alpha=np.arccos(np.clip((s*s+a*a-b*b)/(2*s*a),-1,1))
    beta=np.arccos(np.clip((s*s+b*b-a*a)/(2*s*b),-1,1))
    out[cross]=(a*a*alpha+b*b*beta-.5*np.sqrt(q))/np.pi
    return out

def occulted_fraction(d, ratio, u, ring, weights):
    """Foreground disk radius and projected distance in units of eclipsed radius."""
    base=overlap(1.,ratio,d)
    if u==0:return base
    # sqrt(1-r^2) is the superposition of concentric, uniform disks with radii
    # sqrt(1-v^2), integrated uniformly from v=0 through v=1.
    limb=overlap(ring[None,:],ratio,d[:,None])@weights
    return ((1-u)*base+u*limb)/(1-u/3)

def sky(t,ph,q,sense):
    offset,drift,period,e,angle=q[:5]
    om=angle+sense*2*np.pi*t/period
    mean=2*np.pi*(ph-offset-drift*t/10)+np.pi/2-om
    E=mean+e*np.sin(mean)
    for _ in range(4):
        E-=(E-e*np.sin(E)-mean)/(1-e*np.cos(E))
    ax=np.cos(E)-e
    ay=np.sqrt(1-e*e)*np.sin(E)
    x=ax*np.cos(om)-ay*np.sin(om)
    y=ax*np.sin(om)+ay*np.cos(om)
    return np.hypot(x,q[7]*y),y>0

class LightCurve:
    def __init__(self,data,u,third,sense,order=16):
        h,self.mag,self.err,self.group=data
        self.t=(h-TREF)/365.25
        self.ph=((h-EPOCH)/P)%1
        self.u=u
        self.third=third
        self.sense=sense
        node,w=leggauss(order)
        self.ring=np.sqrt(1-((node+1)/2)**2)
        self.weights=w/2

    def predict_at(self,q,t,ph):
        d,front=sky(t,ph,q,self.sense)
        sumr,k,_,J=q[5:9]
        ra=sumr/(1+k)
        rb=k*ra
        in_eclipse=d<sumr
        dip=np.zeros_like(d)
        for where,back,ratio,flux in ((front,ra,k,1.),(~front,rb,1/k,J*k*k)):
            idx=in_eclipse&where
            dip[idx]=flux*occulted_fraction(d[idx]/back,ratio,self.u,self.ring,self.weights)
        # Identical fixed limb coefficient on both stars makes disk integrated flux ratio J*k².
        rel=1-(1-self.third)*dip/(1+J*k*k)
        angle=2*np.pi*(ph-q[0]-q[1]*t/10)
        return -2.5*np.log10(rel)+q[11]*np.cos(2*angle)+q[12]*np.cos(angle)

    def predict(self,q):
        return q[9:11][self.group]+self.predict_at(q,self.t,self.ph)

    def residuals(self,q):
        return (self.predict(q)-self.mag)/self.err

def scores(resid):
    return 2*SCALE**2*(np.hypot(1,resid/SCALE)-1)

def fit(curve,mask=None,start=None,fixed=None):
    if mask is None:mask=np.ones(len(curve.t),bool)
    q=START.copy() if start is None else np.array(start)
    if start is None and curve.sense<0:q[4]*=-1
    keep=np.ones(len(q),bool)
    for idx,v in (fixed or {}).items():
        keep[idx]=False
        q[idx]=v
    def fill(v):
        full=q.copy()
        full[keep]=v
        return full
    def fun(v):
        return curve.residuals(fill(v))[mask]
    sol=least_squares(fun,q[keep],bounds=(LO[keep],HI[keep]),loss='soft_l1',f_scale=SCALE,
                      x_scale='jac',ftol=5e-7,xtol=5e-7,max_nfev=140)
    pars=fill(sol.x)
    resid=curve.residuals(pars)
    loss=scores(resid)
    held=np.ones(len(loss),bool)
    held[mask]=False
    return {'loss':float(loss[mask].sum()),'held_loss':float(loss[held].sum()),
            'n_train':int(len(loss[mask])),'n_hold':int(held.sum()),'evals':sol.nfev,
            'success':bool(sol.success),'rms_mmag':1000*float(np.sqrt(np.mean((resid[mask]*curve.err[mask])**2))),
            'inclination_deg':float(np.degrees(np.arccos(pars[7]))),'pars':pars.tolist()}

def best(curve,mask=None,init=None,fixed=None,starts=True):
    q=START.copy() if init is None else np.array(init)
    if init is None and curve.sense<0:q[4]*=-1
    configs=[q.copy()]
    if starts:
        for k,c in ((.7,.05),(1.25,.08)):
            x=q.copy()
            x[6:8]=[k,c]
            configs.append(x)
    return min((fit(curve,mask,v,fixed) for v in configs),key=lambda f:f['loss'])

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=HERE/'inputs')
    ap.add_argument('--run',choices=['main','cv','sensitivity','profiles'],default='main')
    args=ap.parse_args()
    data=load(args.input)
    if args.run=='main':
        results={}
        for sense in (1,-1,0):
            curve=LightCurve(data,.4,0.,sense)
            result=best(curve,fixed={2:28.3} if sense==0 else None)
            results[str(sense)]=result
            print(sense,result['loss'],result['pars'],flush=True)
        (HERE/'full.json').write_text(json.dumps(results,indent=2)+'\n')
        return
    full=json.loads((HERE/'full.json').read_text())
    results=[]
    if args.run=='cv':
        grp=np.digitize(data[0],[2453500,2455000,2457000,2458300,2459700])
        for sense in (1,-1,0):
            for i in range(6):
                c=LightCurve(data,.4,0,sense)
                # Initial guesses are fixed across folds; held-out values never set the initialization.
                res=best(c,grp!=i,fixed={2:28.3} if sense==0 else None,starts=False)
                results.append(dict(sense=sense,fold=i,**res))
                print('cv',sense,i,res['held_loss'],flush=True)
    if args.run=='sensitivity':
        for u,l in [(0.,0.),(.7,0.),(.4,.10),(.4,.20)]:
            for sense in (1,-1):
                c=LightCurve(data,u,l,sense)
                res=best(c,init=full[str(sense)]['pars'],starts=False)
                results.append(dict(sense=sense,u=u,third_light=l,**res))
                print('sensitivity',u,l,sense,res['loss'],flush=True)
    if args.run=='profiles':
        for sense in (1,-1):
            c=LightCurve(data,.4,0.,sense)
            for i in (77.,79.,81.,83.,85.,87.,89.):
                res=best(c,init=full[str(sense)]['pars'],fixed={7:np.cos(np.radians(i))},starts=False)
                results.append(dict(sense=sense,inclination_fixed=i,**res))
                print('profile',sense,i,res['loss'],flush=True)
    (HERE/(args.run+'.json')).write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
