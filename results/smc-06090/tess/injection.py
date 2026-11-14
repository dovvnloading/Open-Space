"""Check timing recovery by injecting a known two-eclipse signal in the real pixel cubes."""
import hashlib,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from photometry import lightcurves,wrap,PER,T0,RA,DEC
from tess_timing import GROUPS,fit
HERE=Path(__file__).resolve().parent

def folds(t,ph,aps,orb,cadence):
    S=[];N=[];b=np.floor((ph%1)*120).astype(int)
    for c in np.unique(orb):
        ok=orb==c;n=np.bincount(b[ok],minlength=120)
        if n.sum()*cadence<=.65*PER or (n>0).sum()<min(55,int(.65*PER/cadence)):continue
        S.append(np.bincount(b[ok],weights=aps[ok,0],minlength=120))
        N.append(n)
    return np.array(S),np.array(N)

def prepare(path):
    with fits.open(path) as h:
        t=h[1].data['TIME'].astype(float);cube=h[1].data['FLUX'].astype(float);flags=h[1].data['QUALITY']
        x,y=WCS(h[2].header,naxis=2).world_to_pixel_values(RA,DEC)
    yy,xx=np.indices(cube.shape[1:]);sp=np.exp(-.5*((xx-x)**2+(yy-y)**2)/(.7**2));sp/=sp.sum()
    meta={'x':float(x),'y':float(y),'cadence':float(np.median(np.diff(t[np.isfinite(t)])))}
    return t,cube,flags,sp,meta

def run():
    groups={}
    for name,sectors in [('2018_2019',[1,2,13]),('2025',[94,95,96])]:
        sets=[[],[]];nns=[[],[]]
        for s in sectors:
            p=HERE/'raw'/f's{s:03}.fits'
            t,cube,fl,sp,m=prepare(p)
            norig=lightcurves(p,return_native=True)
            # Inject at the same two arbitrary phases in both groups, outside the
            # true dips. Use identical 2 e-/s depth and a Gaussian PSF sigma 0.7 pix.
            centers=[.205,.735];sep=centers[1]-centers[0]
            ph=wrap((t+2457000.-70/86400.-T0)/PER)
            signal=np.exp(-.5*(wrap(ph-centers[0])/.025)**2)+np.exp(-.5*(wrap(ph-centers[1])/.025)**2)
            for arm in [0,1]:
                if arm==0:native=norig
                else:
                    native=lightcurves(p,return_native=True,inject=(signal,sp,2.))
                t2,phase,A,orb=native
                # Project the injected phase origin onto the timing fitter's axes.
                iz=np.floor((t2+2457000.-70/86400.-T0)/PER-.205+.2).astype(int)
                ss,nn=folds(t2,wrap(phase-.205),A,iz,m['cadence'])
                sets[arm].append(ss);nns[arm].append(nn)
        # Subtract the baseline folded flux on the same phase bins; the same reference
        # pixels and independent clipping criterion give identical counts in both arms
        # except when the faint PSF tail alters one clipping boundary.
        remainders=[];counts=[]
        for a,b,nn,mm in zip(sets[0],sets[1],nns[0],nns[1]):
            if a.shape!=b.shape:raise ValueError('Coverage changed')
            summed_a=a.sum(axis=0);summed_b=b.sum(axis=0);n=nn.sum(axis=0);m=mm.sum(axis=0);w=np.minimum(n,m)
            residual=(summed_b/m-summed_a/n)*w
            # Preserve the sector weight (number of usable orbits).
            remainders.append(np.tile(residual/len(a),(len(a),1)));counts.append(np.tile(w/len(a),(len(a),1)))
        q=fit(remainders,counts);groups[name]={'injected_sep':sep,'recovered_sep':float(q[2]),'bias_minutes':float((q[2]-sep)*PER*1440),'rec_depths_e_s':q[3:5].tolist()}
        print(name,groups[name],flush=True)
    (HERE/'injection.json').write_text(json.dumps(groups,indent=2)+'\n')

if __name__=='__main__':run()
