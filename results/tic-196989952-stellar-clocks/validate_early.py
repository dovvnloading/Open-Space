"""Project TESS-defined frequencies onto earlier surveys with night-block sign randomization."""
from pathlib import Path
import json
import numpy as np
from analyze import inputs, outside

HERE=Path(__file__).resolve().parent

def block_test(t, y, freq, subtract_night, seed=196989952):
    """Preserve each nightly residual vector and observation times in a wild sign test."""
    rng=np.random.default_rng(seed)
    day=np.floor(t+.5)
    groups=np.unique(day)
    y=y-np.mean(y)
    t=t-np.mean(t)
    X=np.column_stack([fn(2*np.pi*t*f) for f in freq for fn in (np.sin,np.cos)])
    if subtract_night:
        for key in groups:
            use=day==key
            y[use]-=np.mean(y[use])
            X[use]-=np.mean(X[use],axis=0)
    else:
        X-=X.mean(axis=0)
    Q=np.linalg.qr(X)[0]
    scores=np.stack([y[day==k]@Q[day==k] for k in groups])
    observed=np.sum(scores.sum(axis=0)**2)
    exact=len(groups)<20
    total=2**(len(groups)-1) if exact else 30000
    exceed=0
    for start in range(0,total,1000):
        size=min(1000,total-start)
        if exact:
            # Global sign reversal leaves the quadratic statistic unchanged.
            signs=2*((np.arange(start,start+size)[:,None]>>np.arange(len(groups)))%2)-1
        else:
            signs=rng.choice([-1,1],size=(size,len(groups)))
        exceed+=np.count_nonzero(np.sum((signs@scores)**2,axis=1)>=observed-1e-13)
    return dict(n_points=len(t),n_groups=len(groups),partial_r2=observed/np.sum(y*y),
                exceedances=int(exceed),realizations=total,p=(exceed/total if exact else (1+exceed)/(1+total)))

def main():
    _,w,g=inputs()
    freq=json.loads((HERE/'measurements.json').read_text())['tess_frequencies_cpd']
    t=w.btjd.to_numpy();y=w.flux.to_numpy();k=outside(t)&(abs(y-1)<.09)
    tg=g.TimeG.to_numpy()-1802.5;yg=g.Gmag.to_numpy()
    kg=(g.GrVFlag.to_numpy()==0)&outside(tg)&np.isfinite(yg)
    result={}
    for label,time,flux in [('SuperWASP',t[k],y[k]),('Gaia G',tg[kg],yg[kg])]:
        result[label]={}
        for name,frequencies in [('first',[freq[0]]),('second',[freq[1]]),('joint',freq)]:
            result[label][name]={variant:block_test(time,flux,frequencies,night)
                                for variant,night in [('centered',False),('night_offsets',True)]}
    (HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    main()
