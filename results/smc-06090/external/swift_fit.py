"""Fit Swift UVW1 eclipse depths using phases predicted in advance from OGLE."""
from common import D, np, PARAMS,PERIOD,shape
from astropy.table import Table
from scipy.optimize import least_squares
import csv,json
def gauss(dx,width):
    return np.exp(-.5*(dx/width)**2)
def design(h,off=0):
    _,d1,d2=shape(h+off*PERIOD)
    return gauss(d1,PARAMS[4]),gauss(d2,PARAMS[5])
def fit(y,e,h,year,rate,rerr,ref=False,off=0):
    # Global normalized ratio. Per-year linear nuisance offsets allow changes in sensitivity.
    epochs,idx=np.unique(year,return_inverse=True)
    Z=np.eye(len(epochs))[idx];a,b=design(h,off)
    X=np.column_stack((Z,-a,-b))
    if ref:
        e=np.hypot(e/rate,y*rerr/rate**2);y=y/rate
    # Constant base estimate is used only to express fractional amplitudes and error floor.
    scale=np.median(y[(a<.005)&(b<.005)])
    y,e=y/scale,e/scale
    # Include an additional 6% scatter on the normalized measurements.
    e=np.hypot(e,.06)
    return np.linalg.lstsq(X/e[:,None],y/e,rcond=None)[0],X,y,e
def main():
    rows=list(csv.DictReader((D/'swift-photometry.csv').open()));out={}
    for rad in [2.7,3.5,4.5]:
      for kind in ['target','target_ratio','reference']:
        who='ref2' if kind=='reference' else 'target'
        def fetch(key):return np.array([float(x.get(key) or np.nan) for x in rows])
        h=fetch('hjd_mid');ye=np.array([int(r['start_utc'][:4]) for r in rows])
        x=fetch(f'{who}_{rad}_rate');e=fetch(f'{who}_{rad}_err');ref=fetch(f'ref1_{rad}_rate');re=fetch(f'ref1_{rad}_err')
        texp=fetch(f'{who}_{rad}_exposure');r0exp=fetch(f'ref1_{rad}_exposure')
        # Restrict coverage by exposure time relative to the nominal header integration later if needed; local maps in photometry.
        ok=np.isfinite(x+e+ref+re)&(texp>20)&(r0exp>20)&(ref>10)&(ref<70)&(e>0)&(re>0)
        h,ye,x,e,ref,re=[v[ok] for v in [h,ye,x,e,ref,re]]
        # Include 2012 long pointing with distinct per-year normalization only where it samples continuum and eclipse.
        p,X,y,s=fit(x,e,h,ye,ref,re,kind!='target');E=np.unique(ye)
        # Cluster sandwich covariance, grouping observations on the same UT date (including multiple snapshot images).
        nights=np.unique(np.floor(h-.5));J=X/s[:,None];res=(y-X@p)/s
        inv=np.linalg.pinv(J.T@J);groups=np.floor(h-.5)
        score=np.array([J[groups==g].T@res[groups==g] for g in nights]);cov=inv@score.T@score@inv*len(nights)/(len(nights)-1)*max(1,(len(x)-1)/(len(x)-len(p)));se=np.sqrt(np.maximum(np.diag(cov),0))
        test=[]
        for off in np.arange(1,80)/80:
            p0,x0,y0,e0=fit(x,e,h,ye,ref,re,kind!='target',off)
            test.append(sum(((y0-x0@p0)/e0)**2))
        a,b=design(h)
        rr=dict(n=len(x),nights=len(nights),years=E.tolist(),norm_by_year=p[:-2].tolist(),depths=p[-2:].tolist(),se_cluster=se[-2:].tolist(),counts_core=[int((a>.5).sum()),int((b>.5).sum())],nphase_shifts_scored_better=sum(bool(q<sum(((y-X@p)/s)**2)) for q in test),n_shifts=len(test),latest_mjd=float(h.max()-2400000.5))
        if kind=='target_ratio' and rad==3.5:
            # Secondary/primary independent center corrections common to all years, profiled across both free depths.
            def fn(z):
                _,d1,d2=shape(h)
                return (y-(X[:,:-2]@z[:-4]-z[-4]*gauss(d1-z[-2],PARAMS[4])-z[-3]*gauss(d2-z[-1],PARAMS[5])))/s
            bound=np.r_[np.full(len(p),np.inf),.06,.06]
            f=least_squares(fn,np.r_[p,0.,0.],bounds=(-bound,bound))
            Ji=f.jac;rsi=f.fun;Ai=np.linalg.pinv(Ji.T@Ji)
            sc=np.array([Ji[groups==g].T@rsi[groups==g] for g in nights]);cv=Ai@sc.T@sc@Ai*len(nights)/(len(nights)-1)*(len(x)-1)/(len(x)-len(f.x))
            rr['free_offsets_min']=(f.x[-2:]*PERIOD*1440).tolist();rr['free_offsets_SE_min']=(np.sqrt(np.diag(cv))[-2:]*PERIOD*1440).tolist()
            rr['separation_residual_min']=float((f.x[-1]-f.x[-2])*PERIOD*1440);rr['separation_residual_SE_min']=float(np.sqrt(cv[-2,-2]+cv[-1,-1]-2*cv[-1,-2])*PERIOD*1440)
        out[f'{kind}_{rad}']=rr
    (D/'swift.json').write_text(json.dumps(out,indent=2));print(json.dumps(out['target_ratio_3.5']))
if __name__=='__main__':main()
