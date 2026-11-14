"""Test whether different observing times and changing eclipse shapes mimic the gap."""
import json
import numpy as np
from analyze import BASE, PERIOD, EPOCH, ogle, wrap, measure

def model(data, x):
    p=wrap((data[:,0]-EPOCH)/PERIOD)
    base, h, a, c, w, b, d, z=x
    return (base+h*np.cos(4*np.pi*p)+a*np.exp(-.5*(wrap(p-c)/w)**2)
            +b*np.exp(-.5*(wrap(p-d)/z)**2))

def experiment():
    old,new=ogle("ii"),ogle("iv")
    new=new[new[:,0]>=2459700]
    samples=[old,new]
    fits=[measure(sample) for sample in samples]
    centers=(fits[0][[3,6]]+fits[1][[3,6]])/2
    rng=np.random.default_rng(4602)
    changes=[]
    models=[]
    for sample,fit in zip(samples,fits):
        res=sample[:,1]-model(sample,fit)
        common=fit.copy()
        common[[3,6]]=centers
        _,night=np.unique(np.floor(sample[:,0]+.5), return_inverse=True)
        models.append((model(sample,common),res,night,common))
    for _ in range(300):
        gaps=[]
        for sample,(mean,res,night,start) in zip(samples,models):
            boot=sample.copy()
            # Preserve the sign relationships among residuals on the same night.
            signs=rng.choice([-1.,1.],night.max()+1)
            boot[:,1]=mean+res*signs[night]
            x=measure(boot,start=start)
            gaps.append((x[6]-x[3])*PERIOD*24)
        changes.append(gaps[1]-gaps[0])
    d2=np.loadtxt(BASE/"data/ii-04602-I.dat",usecols=(0,1,2,3))
    d2=d2[(d2[:,3]==0)&(d2[:,2]<=.10)&(d2[:,2]>0),:3]
    fit2=measure(d2)
    out={"fixed_phase_centers":centers.tolist(),"seed":4602,"trials":len(changes),
         "null_delta_hours_95_range":np.percentile(changes,[2.5,97.5]).tolist(),
         "null_delta_hours_all_range":[min(changes),max(changes)],
         "OGLE_II_zero_flags":{"n":len(d2),"hours":(fit2[6]-fit2[3])*PERIOD*24}}
    (BASE/"results/cadence.json").write_text(json.dumps(out,indent=2)+"\n")
    print(out)

if __name__=="__main__":
    experiment()
