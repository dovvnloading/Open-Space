"""Expand four-dip blend optimization with broad frequency and random-phase initializations."""
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import least_squares
from compare import HERE, get_data, fold_id, setup, solve_model, INPUTS

def run(task):
    fold,block,previous=task
    D=get_data()
    train=np.ones(len(D),bool) if fold=='all' else fold_id(D[:,0])!=fold
    T=D[train]
    rng=np.random.default_rng(30000*block+(787 if fold=='all' else 145*(fold+1)))
    base,lo,hi=setup('blend')
    lo[2]=-0.24;hi[2]=.25;lo[3]=.22;hi[3]=.78
    lo[8]=-.11;hi[8]=.11
    begin=np.array(previous)
    def fun(p):return (solve_model(p,T,'blend')-T[:,4])/T[:,5]
    def jac(p):return solve_model(p,T,'blend',True)/T[:,5,None]
    best=None
    for j in range(40):
        if j%4<2:
            p=begin.copy()
            p[[2,3,10,11]]+=rng.normal(0,.028,4)
            p[[8,16]]+=rng.normal(0,.003 if j%2 else .02,2)
            if j%11==0: p[[2,3,4,5,6,7,8,10,11,12,13,14,15,16]]=p[[10,11,12,13,14,15,16,2,3,4,5,6,7,8]]
        else:
            p=base.copy();p[[2,3,10,11]]+=rng.uniform(-.21,.21,4)
            p[[8,16]]=rng.uniform(-.10,.10,2) if j%2 else rng.uniform(-.018,.018,2)
            p[[6,7,14,15]]=rng.uniform(.05,.27,4)
        p=np.clip(p,lo+1e-6,hi-1e-6)
        r=least_squares(fun,p,jac=jac,bounds=(lo,hi),loss='soft_l1',f_scale=1.8,
                        ftol=3e-7,xtol=3e-7,gtol=3e-7,max_nfev=350,x_scale='jac')
        if best is None or r.cost<best.cost: best=r
    r=(solve_model(best.x,D[~train],'blend')-D[~train,4])/D[~train,5]
    return {'fold':fold,'block':block,'train_loss':2*best.cost,'params':best.x.tolist(),
            'test_loss':float(2*1.8**2*np.sum(np.sqrt(1+(r/1.8)**2)-1)) if len(r) else None}
if __name__=='__main__':
    prior={q['fold']:q for q in json.loads((HERE/'blend.json').read_text())['fitted']}
    tasks=[(fold,j,prior[fold]['params']) for fold in ('all',5) for j in range(4)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        rec=list(pool.map(run,tasks))
    print([(q['fold'],round(q['train_loss'],2),q['test_loss']) for q in rec])
    (HERE/'stress.json').write_text(json.dumps(rec,indent=2)+'\n')
