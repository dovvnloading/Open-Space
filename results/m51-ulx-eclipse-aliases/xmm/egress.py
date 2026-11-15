"""Localize the candidate 2018 transition with independent XMM detectors."""
import numpy as np,pandas as pd
from analyze import *
from scipy.stats import chi2

def profile(v,grid):
    prefix,p,g=v
    t=p.mjd_tdb_geo.to_numpy();x=p.x_arcsec.to_numpy();y=p.y_arcsec.to_numpy()
    good=(p.energy_ev.between(900,7000)).to_numpy()
    # Exposure interval fixed to first ~0.5 d to avoid later intrinsic variability.
    low,high=58251.92,58252.43
    gt=np.column_stack([np.maximum(g[:,0],low),np.minimum(g[:,1],high)])
    gt=gt[gt[:,1]>gt[:,0]]
    good&=(t>=low)&(t<=high)&(abs(x-NEIGHBOR[0]/2)<19)&(abs(y-NEIGHBOR[1]/2)<19)
    t,x,y=t[good],x[good],y[good]
    rc,alpha,half=4.5,1.5,19
    gx=np.linspace(-half,half,101);xx,yy=np.meshgrid(gx+NEIGHBOR[0]/2,gx+NEIGHBOR[1]/2)
    norm=components(xx.ravel(),yy.ravel(),rc,alpha,half).mean(axis=1)*(2*half)**2
    X=components(x,y,rc,alpha,half)/norm[:,None]
    total_exp=np.sum(gt[:,1]-gt[:,0])*86.4
    res=[]
    for b in grid:
        pre=t<b
        dur=np.array([np.maximum(0,np.minimum(gt[:,1],b)-gt[:,0]).sum()*86.4,0]);dur[1]=total_exp-dur[0]
        # Neighbor and background free on each side, so the alleged change need not
        # be carried by ULX-1. Extended Poisson likelihood includes exposure.
        rows=[]
        for j,mask in enumerate([pre,~pre]):
            def fun(q):
                model=np.maximum(q@X[:,mask],1e-14)
                return dur[j]*sum(q)-np.log(model).sum(),dur[j]-(X[:,mask]/model).sum(axis=1)
            result=minimize(fun,np.ones(6)*mask.sum()/dur[j]/6,jac=True,method='L-BFGS-B',bounds=[(0,None)]*6,options={'ftol':1e-12})
            rows.append((result.fun,*result.x[:2]))
        res.append(dict(detector=prefix[11:13],mjd=b,minuslogL=rows[0][0]+rows[1][0],pre_ks=dur[0],post_ks=dur[1],ulx1_pre_rate=rows[0][1],ulx1_post_rate=rows[1][1],ulx2_pre_rate=rows[0][2],ulx2_post_rate=rows[1][2]))
    return res
if __name__=='__main__':
    grid=np.arange(58252.045,58252.199,.00025)
    rows=[]
    for ff in sorted((HERE/'inputs').glob('P0824450901*-photons.csv.gz')):
        rows+=profile(read_visit(ff),grid)
    f=pd.DataFrame(rows)
    f.to_csv(HERE/'xmm-egress-profile.csv.gz',index=False,float_format='%.8f')
    joint=f.groupby('mjd').minuslogL.sum()
    j=joint.idxmin();s=joint[joint<joint.min()+1.92]
    print('combined',j,s.index.min(),s.index.max())
    for name,part in f.groupby('detector'):
        ix=part.minuslogL.idxmin();z=part.loc[ix];print(name,z[['mjd','ulx1_pre_rate','ulx1_post_rate','ulx2_pre_rate','ulx2_post_rate']].round(4).tolist())
