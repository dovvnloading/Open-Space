"""Compare Gaia's photocentre to the photocentric-orbit shape scaled freely by the RVs."""
import csv
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lamost'))
from check_orbits import predictions, unpack
from scipy.special import logsumexp

SID = '3378588057203660160'


def run():
    row = next(r for r in csv.DictReader((ROOT/'orbit/lamost_candidate_nss.csv').open()) if r['source_id']==SID)
    star = next(r for r in csv.DictReader((ROOT/'gaia/mass_screen.csv').open()) if r['source_id']==SID)
    obs = [r for r in csv.DictReader((ROOT/'lamost/exposures.csv').open()) if r['source_id']==SID]
    times = np.array([float(r['mjd_utc_start']) for r in obs])
    nights = np.array([r['obsdate'] for r in obs])
    err = np.array([float(r['err_kms']) for r in obs])
    vel = np.array([float(r['rv_kms']) for r in obs])
    dates = np.unique(nights)
    chunks=[]
    for j in range(8):
        seed=230021+j
        shape=predictions(row,times,np.random.default_rng(seed))[1:]
        draw=unpack(row,4096,np.random.default_rng(seed))
        A,B,F,G = [draw[c+'_thiele_innes'][1:] for c in 'abfg']
        a0 = (np.hypot(A+G,B-F)+np.hypot(A-G,B+F))/2
        p=draw['period'][1:]
        e=draw['eccentricity'][1:]
        plx=draw['parallax'][1:]
        assert len(shape)==len(p) and np.all((e>=0)&(e<1)&(plx>0)&(p>0))
        f0=(a0/plx)**3/(p/365.25)**2
        means=[np.average(shape[:,nights==date],weights=err[nights==date]**-2,axis=1) for date in dates]
        chunks.append((f0,np.stack(means,axis=1)))
    f0=np.concatenate([q[0] for q in chunks]); x=np.concatenate([q[1] for q in chunks])
    y=np.array([np.average(vel[nights==d],weights=err[nights==d]**-2) for d in dates])
    var=np.array([1/np.sum(err[nights==d]**-2) for d in dates])
    rng=np.random.default_rng(42301)
    out={}
    for floor in (1,2):
        w=1/(var+floor**2)
        xc=x-np.average(x,axis=1,weights=w)[:,None]
        yc=y-np.average(y,weights=w)
        scaling_info=np.sum(xc**2*w,axis=1)
        scale=np.sum(xc*yc*w,axis=1)/scaling_info
        chi2=np.sum((yc-scale[:,None]*xc)**2*w,axis=1)
        # Marginalize over the signed RV amplitude scale and intercept with flat priors.
        lw=-0.5*(chi2+np.log(scaling_info))
        weights=np.exp(lw-logsumexp(lw))
        i=rng.choice(len(f0),50000,p=weights)
        s=np.abs(scale[i]+rng.normal(size=len(i))/np.sqrt(scaling_info[i]))
        f=f0[i]*s**3
        z=rng.normal(size=len(f))
        m1=float(star['m1'])+z*np.where(z<0,float(star['m1'])-float(star['m1_lower']),
                                           float(star['m1_upper'])-float(star['m1']))
        low,high=f,f+2*m1
        for _ in range(50):
            m=(low+high)/2; large=m**3>f*(m1+m)**2
            low,high=np.where(large,low,m),np.where(large,m,high)
        Bmass=m/(m1+m)
        beta_same=Bmass*(1-1/s); beta_reverse=Bmass*(1+1/s)
        quant=lambda a:np.quantile(a,[.025,.16,.5,.84,.975]).tolist()
        out[str(floor)]=dict(importance_ess=1/sum(weights**2),m2=quant(m),
                a_phot_over_a1=quant(1/s),beta_same_signed=quant(beta_same),
                beta_reversed=quant(beta_reverse),reversed_fraction_phys=np.mean(beta_reverse<=1))
    return dict(source_id=SID,draws=len(f0),nights=len(dates),quantiles=[.025,.16,.5,.84,.975],
                by_nightly_floor_kms=out)


if __name__ == '__main__':
    (Path(__file__).parent/'photocentre.json').write_text(json.dumps(run(),indent=2)+'\n')
