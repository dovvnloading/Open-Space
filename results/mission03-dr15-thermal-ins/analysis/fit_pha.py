"""Fit a PPS source/background PHA using its ARF and RMF and joint Poisson counts."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.special import xlogy
from astropy.io import fits

# Morrison & McCammon (1983), ApJ 270, 119, table 2. sigma in cm2 per H atom.
BREAKS=np.array([.03,.1,.284,.4,.532,.707,.867,1.303,1.840,2.471,3.210,4.038,7.111,8.331])
C0=np.array([17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,202.7,342.7,352.2,433.9,629.0,701.2])
C1=np.array([608.1,267.9,18.8,66.8,145.8,-380.6,169.3,146.8,104.7,18.7,18.7,-2.4,30.9,25.2])
C2=np.array([-2150.,-476.1,4.3,-51.4,-61.1,294.,-47.7,-31.5,-17.,0.,0.,0.,0.,0.])

def sigma(e):
    i=np.maximum(0,np.searchsorted(BREAKS,e,side='right')-1)
    return 1e-24*(C0[i]+e*C1[i]+e*e*C2[i])/e**3

def photon_shape(e,shape,param):
    if shape=='blackbody': return e**2/np.expm1(np.minimum(e/param,700))
    if shape=='brems_continuum': return np.exp(-e/param)/e
    return e**(-param)

def profile_stat(n,b,alpha,s):
    # Profile nuisance expectation beta in background aperture; n~Pois(s+a*beta), b~Pois(beta).
    linear=(1+alpha)*s-alpha*(n+b)
    beta=np.maximum((-linear+np.sqrt(linear**2+4*alpha*(1+alpha)*b*s))/(2*alpha*(1+alpha)),1.e-200)
    mu=np.maximum(s+alpha*beta,1.e-200)
    return 2*np.sum(mu-n+xlogy(n,n/mu)+beta-b+xlogy(b,b/beta))

def read(path,minimum):
    with fits.open(path) as f:
        h=f['SPECTRUM']
        meta=h.header.copy(); chan=h.data['CHANNEL'].astype(int)
        n=h.data['COUNTS'].astype(float)
        good=h.data['QUALITY']==0
    with fits.open(path.parent/meta['BACKFILE']) as f:
        h=f['SPECTRUM']; b=h.data['COUNTS'].astype(float)
        if not np.array_equal(chan,h.data['CHANNEL']):raise ValueError('background channel grid mismatch')
        alpha=meta['EXPOSURE']*meta['BACKSCAL']/(h.header['EXPOSURE']*h.header['BACKSCAL'])
    p=path.parent/meta['RESPFILE']
    if not p.exists():p=p.with_name(p.name+'.gz')
    with fits.open(p) as f:
        rmf=f['MATRIX'].data; eb=f['EBOUNDS'].data
        true_lo, true_hi = np.asarray(rmf['ENERG_LO'],float),np.asarray(rmf['ENERG_HI'],float)
        e=(true_hi+true_lo)/2
        if not np.array_equal(eb['CHANNEL'],chan):raise ValueError('PHA/RMF channels differ')
        eb_mid=(eb['E_MIN']+eb['E_MAX'])/2
        bounds=np.unique(np.round(np.r_[np.arange(minimum,1.001,.05),np.arange(1.1,2.001,.1),2.5,3,4,5],6))
        bins=np.searchsorted(bounds,eb_mid,side='right')-1
        nb=len(bounds)-1
        legit=good&(bins>=0)&(bins<nb)
        sn=np.bincount(bins[legit],weights=n[legit],minlength=nb)
        bn=np.bincount(bins[legit],weights=b[legit],minlength=nb)
        m=np.zeros((nb,len(e)))
        for i,row in enumerate(rmf):
            k=int(row['N_GRP'])
            if not k: continue
            starts=np.atleast_1d(row['F_CHAN'])[:k].astype(int)
            lengths=np.atleast_1d(row['N_CHAN'])[:k].astype(int)
            c=np.concatenate([np.arange(s,s+z) for s,z in zip(starts,lengths)])
            channel_idx=c-chan[0]; mask=legit[channel_idx]
            vals=np.asarray(row['MATRIX'])[:len(c)]
            m[:,i]=np.bincount(bins[channel_idx[mask]],weights=vals[mask],minlength=nb)
    with fits.open(path.parent/meta['ANCRFILE']) as f:
        a=f['SPECRESP'].data
        if not (np.allclose(true_lo,a['ENERG_LO']) and np.allclose(true_hi,a['ENERG_HI'])):raise ValueError('ARF and RMF energy grids differ')
        m*=a['SPECRESP']*(true_hi-true_lo)*meta['EXPOSURE']
    # Drop bins with no good channels. Ignore out-of-calibration true energies below 0.03 keV.
    keep=np.bincount(bins[legit],minlength=nb)>0; te=(e>=.03)&(e<=10.)
    return e[te],m[keep][:,te],sn[keep],bn[keep],float(alpha),bounds[:-1][keep],bounds[1:][keep],true_lo[te],true_hi[te]

def run(path,minimum,fixed_nh=None):
    energy,m,n,b,alpha,lo,hi,_,_=read(path,minimum)
    nh_scale=sigma(energy)*1.e21
    result={'file':path.name,'energy_edges_keV':[[float(l),float(h)] for l,h in zip(lo,hi)],'source_counts':n.astype(int).tolist(),'background_counts':b.astype(int).tolist(),'background_scale':alpha,'models':{}}
    for model,shape_bounds in [('blackbody',[.02,.8]),('powerlaw',[.5,8.]),('brems_continuum',[.025,10.])]:
        limits=[(-8,16),(0,10) if fixed_nh is None else (fixed_nh,fixed_nh),tuple(np.log(shape_bounds))]
        def stat(x):
            total,nh,theta=np.exp(x[0]),x[1],np.exp(x[2])
            shape=photon_shape(energy,model,theta)
            unit=m@(np.exp(-nh*nh_scale)*shape)
            counts=total*unit/max(np.sum(unit),1e-200)
            return profile_stat(n,b,alpha,counts)
        starts=np.geomspace(*shape_bounds,5)[1:-1]
        fits_=[minimize(stat,[np.log(max(1.,n.sum()-alpha*b.sum())),nh,np.log(st)],bounds=limits,method='L-BFGS-B') for st in starts for nh in [0.,1.]]
        best=min(fits_,key=lambda z:z.fun)
        total,nh,theta=np.exp(best.x[0]),best.x[1],np.exp(best.x[2])
        norm=total/np.sum(m@(np.exp(-nh*nh_scale)*photon_shape(energy,model,theta)))
        def profile(theta_log):
            def objective(remaining):return stat([remaining[0],remaining[1],theta_log])
            pair=minimize(objective,best.x[:2],bounds=limits[:2],method='L-BFGS-B')
            return pair.fun-best.fun-2.706
        ci=[]
        for end in limits[-1]:
            ci.append(float(np.exp(brentq(profile,min(best.x[2],end),max(best.x[2],end)))) if abs(end-best.x[2])>1e-6 and profile(end)>0 else None)
        egrid=np.linspace(.2,2,4001)
        energy_flux=1.602176634e-9*np.trapezoid(norm*photon_shape(egrid,model,theta)*np.exp(-nh*1e21*sigma(egrid))*egrid,x=egrid)
        result['models'][model]={'deviance':float(best.fun),'normalization':norm,'NH_1e21':nh,'shape':theta,'shape_profile90':ci,'absorbed_flux_0.2_2':float(energy_flux),'optimizer_converged':bool(best.success)}
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectrum',type=Path)
    p.add_argument('--minimum',type=float,choices=(.2,.3),default=.2)
    p.add_argument('--fixed-nh',type=float,help='Freeze column in units of 1e21 cm^-2')
    a=p.parse_args()
    print(json.dumps(run(a.spectrum,a.minimum,a.fixed_nh),indent=2))
if __name__=='__main__':main()
