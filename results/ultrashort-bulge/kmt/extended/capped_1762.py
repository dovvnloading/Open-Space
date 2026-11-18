"""Profile 2024-only DIAPL night bins for 1762 with an empirical total-light bound."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution,minimize

from tier2 import GEOMETRY
from recurrence import binned
from support import HERE,engine


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True);args=p.parse_args()
    event='KMT-2024-BLG-1762'
    raw=pd.read_csv(args.cache/event/'BLG31M0105.058473.cache')
    sites=binned(raw,GEOMETRY[event])[2024]
    cal=json.loads((HERE/'diapl-scale.json').read_text())
    caps={x['site']:x['largest_calibrator_p84_counts_I20'] for x in cal['transfer']}
    order=list(sites)
    t=np.concatenate([sites[s][0] for s in order])-60512.47902
    y=np.concatenate([sites[s][1] for s in order]);w=np.concatenate([1/sites[s][2]**2 for s in order])
    g=np.repeat(np.arange(len(order)),[len(sites[s][0]) for s in order])
    wg=np.bincount(g,weights=w);ym=np.bincount(g,weights=w*y)/wg
    vb=engine();fits=[];previous=[]
    for multiplier in (1.,2.,float('inf')):
        hi=np.array([caps[s]*multiplier for s in order])
        for model,rb in [('ps',None),('fs',(np.log(.02),np.log(15.))),('large',(0.,np.log(15.))),('small',(np.log(.02),0.))]:
            domain=[(-50,50),(np.log(.15),np.log(200)),(np.log(.002),np.log(5.))]+([rb] if rb else [])
            def evaluate(x,params=False):
                center,lt,lu=x[:3];te,u=np.exp([lt,lu])
                if model=='ps':
                    q=((t-center)/te)**2+u*u;a=(q+2)/np.sqrt(q*(q+4))-1
                else:
                    a=np.asarray(vb.ESPLLightCurve([lu,lt,center,x[3]],t.tolist())[0])-1
                am=np.bincount(g,weights=w*a)/wg
                ah=a-am[g]
                den=np.bincount(g,weights=w*ah**2)
                fs=np.clip(np.bincount(g,weights=w*ah*y)/np.maximum(den,1e-30),0,hi)
                res=y-ym[g]-fs[g]*ah
                return dict(zip(order,map(float,fs))) if params else float(np.sum(w*res**2))
            opt=differential_evolution(evaluate,domain,seed=527,popsize=7,maxiter=60,polish=False,tol=.0003)
            starts=[opt.x,[0.,np.log(2.596),np.log(.002)]+([np.clip(np.log(15.),*rb)] if rb else [])]
            if rb:
                for T,U,R in [(34.,.15,1.55),(46.,.15,.96),(40.,1.05,1.2),(24.,.05,2.)]:
                    starts.append([0.,np.log(T),np.log(U),np.clip(np.log(R),*rb)])
                starts.extend([np.clip(x,[a for a,b in domain],[b for a,b in domain]) for x in previous])
            fit=min([minimize(evaluate,x,method='Nelder-Mead',bounds=domain,
                      options=dict(maxiter=1200,xatol=.00005,fatol=.001)) for x in starts],key=lambda y:y.fun)
            x=fit.x
            if rb:previous.append(x)
            fits.append(dict(cap_factor=multiplier if np.isfinite(multiplier) else None,model=model,
                chi2=round(float(fit.fun),3),t0=round(float(x[0])+60512.47902,4),te=round(float(np.exp(x[1])),4),u0=round(float(np.exp(x[2])),4),
                rho=round(float(np.exp(x[3])),4) if rb else 0.,source_counts=evaluate(x,True)))
    (HERE/'capped-1762.json').write_text(json.dumps(dict(n=len(t),sites=order,empirical_total_I20_counts=caps,fits=fits),indent=2)+'\n')


if __name__=='__main__':main()
