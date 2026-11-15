"""Fit independent XMM imaging data against phase windows from Chandra."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2
EPOCH=56180.216165685095
INGRESS=56191.61960343367
LATE_EGRESS=59373.11583257826
# Offset of ULX-2 relative to ULX-1 from published Chandra coordinates, arcsec.
NEIGHBOR=np.array([-5.09483324,6.7])
HERE=Path(__file__).resolve().parent

def components(x,y,rc,alpha,half_size):
    x2,y2=NEIGHBOR
    def psf(dx,dy):
        return (1+(dx**2+dy**2)/rc**2)**(-alpha)
    a=(x-x2/2)/half_size
    b=(y-y2/2)/half_size
    return np.vstack([psf(x,y),psf(x-x2,y-y2),(1+a)*(1+b),(1+a)*(1-b),(1-a)*(1+b),(1-a)*(1-b)])

def fit(x,y,rc=4.5,alpha=1.5,half_size=19):
    counts=len(x)
    if counts<25:
        return None
    # Midpoint integration normalizes each component over the actual fit box.
    grid=(np.arange(90)+.5)*2*half_size/90-half_size
    gx,gy=np.meshgrid(grid+NEIGHBOR[0]/2,grid+NEIGHBOR[1]/2)
    norm=components(gx.ravel(),gy.ravel(),rc,alpha,half_size).mean(axis=1)*(2*half_size)**2
    design=components(x,y,rc,alpha,half_size)/norm[:,None]
    def solve(active):
        design_a=design[active]
        def objective(params):
            model=np.maximum(params@design_a,1e-12)
            return params.sum()-np.log(model).sum(),1-(design_a/model).sum(axis=1)
        result=minimize(objective,np.full(len(active),counts/len(active)),jac=True,method='L-BFGS-B',bounds=[(0,None)]*len(active),options={'ftol':1e-12,'maxiter':600})
        if not result.success:
            raise ValueError(result.message)
        return result.x,result.fun
    pars,logl=solve(np.arange(6))
    _,null=solve(np.arange(1,6))
    h=design/np.maximum(pars@design,1e-12)
    err=np.sqrt(np.linalg.pinv(h@h.T)[0,0])
    return pars[0],err,pars[1],2*(null-logl)

def windows(gti,period,duration,margin,in_eclipse):
    low=gti[:,0].min()
    high=gti[:,1].max()
    if in_eclipse:
        begin,end=-duration+margin,-margin
    else:
        begin,end=margin,period-duration-margin
    rows=[]
    for k in range(int((low-EPOCH)//period)-1,int((high-EPOCH)//period)+3):
        epoch=EPOCH+k*period
        l=np.maximum(gti[:,0],epoch+begin)
        r=np.minimum(gti[:,1],epoch+end)
        rows.extend(zip(l[r>l],r[r>l]))
    return np.array(rows).reshape(-1,2)

def read_visit(file):
    prefix=file.name.split('-')[0]
    photons=pd.read_csv(file)
    gti=pd.read_csv(file.parent/(prefix+'-gti.csv')).values
    return prefix,photons,gti

def measure(visit,k,rc=4.5,alpha=1.5,dx=0.,dy=0.,lo=900,hi=7000,half_size=19,margin=.08):
    prefix,photons,gti=visit
    period=(LATE_EGRESS-EPOCH)/k
    family=1 if k<300 else 2
    dur=family*period-(INGRESS-EPOCH)
    xy=photons[['x_arcsec','y_arcsec']].values-np.array([dx,dy])
    x,y=xy.T
    m=(abs(x-NEIGHBOR[0]/2)<half_size)&(abs(y-NEIGHBOR[1]/2)<half_size)&photons.energy_ev.between(lo,hi)
    t=photons.mjd_tdb_geo.to_numpy()
    result=[]
    for label in ['in','out']:
        w=windows(gti,period,dur,margin,label=='in')
        exposure=(w[:,1]-w[:,0]).sum()*86.4
        if exposure<1:continue
        ii=np.clip(np.searchsorted(w[:,0],t)-1,0,len(w)-1)
        mask=m&(t>=w[ii,0])&(t<=w[ii,1])
        model=fit(x[mask],y[mask],rc,alpha,half_size)
        if model is None:continue
        source,error,neighbor,dc=model
        result.append(dict(k=k,period_days=period,duration_days=dur,visit=prefix[1:11],detector=prefix[11:13],state=label,exposure_ks=exposure,box_events=int(mask.sum()),ulx1_counts=source,ulx1_count_error=error,ulx2_counts=neighbor,delta_c=dc,ulx1_rate_per_ks=source/exposure,rate_error_per_ks=error/exposure))
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,default=HERE/'inputs')
    parser.add_argument('--output',type=Path,default=HERE/'xmm-alias-rates.csv')
    parser.add_argument('--cycles',nargs='+',type=int,default=[*range(243,258),*range(484,515)])
    args=parser.parse_args()
    visits=[read_visit(f) for f in sorted(args.inputs.glob('*-photons.csv.gz'))]
    rows=[r for v in visits for k in args.cycles for r in measure(v,k)]
    pd.DataFrame(rows).to_csv(args.output,index=False,float_format='%.7g')
if __name__=='__main__':main()
