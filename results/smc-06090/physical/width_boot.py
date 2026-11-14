"""Propagate paired observing-night resampling into the 2027 signed width forecast."""
import json
import numpy as np
from geometry import HERE,LightCurve,load,fit
from predictions import evaluate

def main(draws=120):
    data=load(HERE/'inputs')
    d=json.loads((HERE/'dilution.json').read_text())
    curves=[(LightCurve(data,.4,third,sense,order=16),
             np.array(next(r for r in d[str(sense)] if r['third']==third)['pars']))
            for sense in (1,-1) for third in (0.,.3)]
    rng=np.random.default_rng(523031)
    h=data[0];night=np.floor(h-.5)
    block=np.digitize(h,[2453500,2455000,2457000,2458300,2459700])
    groups=[[np.where((block==g)&(night==n))[0] for n in np.unique(night[block==g])] for g in range(6)]
    out=[[] for _ in curves]
    for i in range(draws):
        idx=np.concatenate([g[j] for g in groups for j in rng.integers(len(g),size=len(g))])
        for j,(c,p) in enumerate(curves):
            res=fit(c,idx,p)
            forecast=evaluate(c,np.array(res['pars']),2461587.5,True)
            out[j].append([forecast['difference_min'],forecast['separation_hr']])
        if i%20==0:print(i,flush=True)
    entries=[]
    for (c,p),samples in zip(curves,out):
        entries.append(dict(sense=c.sense,third=c.third,u=.4,
                            columns=['secondary_minus_primary_FWHM_minutes','separation_hours'],
                            quantiles=np.percentile(samples,[2.5,50,97.5],axis=0).tolist(),draws=samples))
    (HERE/'width_boot.json').write_text(json.dumps(entries,indent=2)+'\n')
if __name__=='__main__':main()
