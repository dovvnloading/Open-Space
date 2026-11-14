"""Compare eclipse separation histories of OGLE-SMC-ECL-06090 using individual I measurements.

Run with OPENBLAS_NUM_THREADS=1 and Python / NumPy / SciPy. No network access is required.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
EP=2460000.0369
P=1.1696719
TC=2456000
Y=365.25
SCALE=1.8
HASHES={'ogle3.dat':'d57614117ef140adf1413baab42f92182b9c39752eed39706b6fc040e203a897','ogle4.dat':'429caaf74f2086f3a8f53c71a91bfb3c07049dbc2e81122c06646387c8fa21a0'}
def get():
    data=[]
    for i,n in enumerate(HASHES):
        f=HERE.parent/'joint'/'inputs'/n;assert hashlib.sha256(f.read_bytes()).hexdigest()==HASHES[n]
        a=np.loadtxt(f)[:,:3]
        if i==0:a[:,0]+=2450000
        a=a[np.isfinite(a).all(axis=1)&(a[:,2]>0)&(a[:,2]<=.1)]
        data.append(a)
    h,m,e=np.concatenate(data).T
    return h,(h-TC)/Y,(h-EP)/P%1,m,np.hypot(e,.012)
def basis(t,model,period):
    if model=='constant':return np.ones((len(t),1))
    if model=='linear':return np.array([np.ones(len(t)),t/10]).T
    if model=='quadratic':return np.array([np.ones(len(t)),t/10,(t/10)**2]).T
    if model=='cubic':return np.array([np.ones(len(t)),t/10,(t/10)**2,(t/10)**3]).T
    wave=np.array([np.cos(2*np.pi*t/period),np.sin(2*np.pi*t/period)]).T
    return wave if model=='apsidal' else np.c_[np.ones(len(t)),wave]
def circ(q):return (q+.5)%1-.5
class Model:
    def __init__(self,h,t,ph,m,e,design,group,shape):
        self.t=t; self.ph=ph;self.m=m;self.e=e;self.group=group;self.B=design
        self.k=design.shape[1];self.shape=shape;self.groups=int(max(group)+1);self.ell=7 if shape=='asymmetric' else 5
        self.start=np.r_[.015,0.,np.zeros(self.k),np.tile([18.23,.41,.023,.40,.023]+([0.,0.] if self.ell==7 else []),self.groups)]
        # Bounds allow much greater eclipse timing variation than is measured.
        self.lo=np.r_[-.15,-.12,[-.23]*self.k,np.tile([18.00,.12,.010,.12,.010]+([-.45,-.45] if self.ell==7 else []),self.groups)]
        self.hi=np.r_[.15,.12,[.23]*self.k,np.tile([18.45,.9,.065,.9,.065]+([.45,.45] if self.ell==7 else []),self.groups)]
    def eval(self,a):
        q=a[2+self.k:].reshape(-1,self.ell)[self.group]
        level,d1,w1,d2,w2=q[:,:5].T
        shared=a[0]+a[1]*self.t/10
        half=.5*(self.B@a[2:self.k+2])
        r1,r2=circ(self.ph-(shared-half)),circ(self.ph-(.5+shared+half))
        def one(r,d,w,skew):
            rs=np.sign(r)
            fac=np.exp(skew*rs)
            z=r/(w*fac)
            if self.shape=='triangular':
                term=np.maximum(1.-abs(z)/2.5,0)
                dc=np.where(term>0,d*np.sign(z)/(2.5*w*fac),0)
                dw=np.where(term>0,d*abs(z)/(2.5*w),0)
            else:
                term=np.exp(-.5*z*z)
                dc=d*term*z/(w*fac)
                dw=d*term*z*z/w
            return term,dc,dw,dw*w*rs
        f1,c1,v1,s1=one(r1,d1,w1,q[:,5] if self.ell==7 else 0.)
        f2,c2,v2,s2=one(r2,d2,w2,q[:,6] if self.ell==7 else 0.)
        self.pred=level+d1*f1+d2*f2
        self.jac=np.zeros((len(q),len(a)))
        self.jac[:,0]=c1+c2
        self.jac[:,1]=self.jac[:,0]*self.t/10
        self.jac[:,2:2+self.k]=.5*(c2-c1)[:,None]*self.B
        ids=np.arange(len(q));ind=2+self.k+self.ell*self.group
        for j,v in enumerate([1.,f1,v1,f2,v2]+([s1,s2] if self.ell==7 else [])): self.jac[ids,ind+j]=v
        self.cache=a.copy()
    def fun(self,a):
        self.eval(a)
        return (self.pred-self.m)/self.e
    def derivative(self,a):
        if not np.array_equal(a,self.cache):self.eval(a)
        return self.jac/self.e[:,None]
def estimate(data,model,period=25.,mask=None,shape='gaussian',group=None,initial=None):
    h,t,ph,m,e=data
    if group is None:group=(h>=2455000).astype(int)
    if mask is None:mask=np.ones(len(h),bool)
    b=basis(t,model,period)
    o=Model(*(v[mask] for v in data),b[mask],group[mask],shape)
    # A leave-out fold can omit one survey entirely; such folds are never used here.
    sol=least_squares(o.fun,o.start if initial is None else initial,jac=o.derivative,
        bounds=(o.lo,o.hi),loss='soft_l1',f_scale=SCALE,max_nfev=130,ftol=1e-8)
    oo=Model(*data,b,group,shape); resid=oo.fun(sol.x)
    loss=2*SCALE**2*(np.sqrt(1+(resid/SCALE)**2)-1)
    return {'deviance':loss.sum(),'loss':loss,'x':sol.x,'pred':oo.pred}
def profile(data,model,periods,mask=None,group=None,shape='gaussian'):
    # Uniform frequency samples cover short aliases. Tighten every local minimum.
    y=np.array([estimate(data,model,T,mask,shape,group)['loss'][mask].sum() for T in periods])
    local=np.r_[0,find_peaks(-y)[0],len(y)-1].astype(int);opts=[]
    for i in local:
        if i in (0,len(y)-1): pp=periods[i];loss=y[i]
        else:
            r=minimize_scalar(lambda p:estimate(data,model,p,mask,shape,group)['loss'][mask].sum(),
                              bounds=(periods[i-1],periods[i+1]),method='bounded',options={'xatol':.002})
            pp,loss=r.x,r.fun
        opts.append([float(pp),float(loss)])
    opts.sort(key=lambda z:z[1])
    return opts,np.c_[periods,y]
def foldmask(data):
    h=data[0]
    return np.digitize(h,[2453500,2455000,2457000,2458300,2459700])
def simple(x):return {'deviance':float(x['deviance']),'theta':x['x'].tolist()}
def run():
    data=get();mask=np.ones(len(data[0]),bool); groups=foldmask(data); periods=1/np.linspace(1/200.,1/2.5,165)[::-1]
    out={'revision':'1', 'hashes':HASHES,'n':len(mask),'loss':'sum of twice the soft-L1 objective on standardized individual magnitudes, scale 1.8, photometric error floor 0.012 mag','models':{}}
    for m in ('constant','linear','quadratic','cubic'):out['models'][m]=simple(estimate(data,m))
    for m in ('apsidal','free_offset'):
        opt,grid=profile(data,m,periods,mask)
        out['models'][m]={**simple(estimate(data,m,opt[0][0])),'best_period_yr':opt[0][0],'local_minima':opt}
        np.savetxt(HERE/(m+'-profile.csv'),grid,delimiter=',',header='period_years,robust_deviance',comments='')
    u=out['models']['apsidal']['best_period_yr']
    for shape in ['triangular','asymmetric']:
        opt,_=profile(data,'apsidal',np.linspace(.75*u,1.25*u,14),mask,shape=shape)
        out[shape]={'local_period_and_deviance':opt[0]}
    # Independent survey amplitudes and widths instead allowed to change in six epochs.
    opt,_=profile(data,'apsidal',np.linspace(.75*u,1.25*u,14),mask,group=groups)
    out['six_shape_groups']={'local_period_and_deviance':opt[0]}
    return out
if __name__=='__main__':
    out=run();(HERE/'orbit.json').write_text(json.dumps(out,indent=2)+'\n')
