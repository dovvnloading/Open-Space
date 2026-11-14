"""Challenge the separate-clock model at +/-1 sidereal-day and +/-1 year aliases.

A linear local correction to each alias frequency is optimized from scratch per
holdout. The orbital-phase catalog supplies the reference only. Initial alias
phases are measured from dim training samples folded on the ordinary period;
these are starting points, not labels used by the full-photometry loss.
"""
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import least_squares
from compare import get_data, fold_id, wrapped, PERIOD, INPUTS, EPOCH, TIME_ORIGIN
HERE=Path(__file__).resolve().parent
DAY=1.00273790935
ALIASES={'orbital':0.,'annual_plus':1/365.25,'annual_minus':-1/365.25,
         'daily_plus':DAY,'daily_minus':-DAY}
def job(task):
    fold,la,lb=task
    D=get_data()
    tr=np.ones(len(D),bool) if fold=='all' else fold_id(D[:,0])!=fold
    ni=np.flatnonzero(tr)
    f1,f2=1/PERIOD+ALIASES[la],1/PERIOD+ALIASES[lb]
    x1=((D[:,0]-TIME_ORIGIN)*f1)%1; x2=((D[:,0]-TIME_ORIGIN)*f2)%1
    y=D[:,4];e=D[:,5];yr=D[:,1];sur=D[:,3];h=np.cos(4*np.pi*D[:,2])
    dim=(y<.84)&tr
    bases=[]
    for xx,center in [(x1,.02),(x2,.51)]:
        zz=np.exp(2j*np.pi*xx[dim&(abs(wrapped(D[:,2]-center))<.12)])
        bases.append(float(np.angle(zz.mean())/(2*np.pi)))
    def calc(p,idx,jac=False):
        z=np.where(sur[idx]==0,p[0],p[1])+p[10]*h[idx]
        if jac:
            J=np.zeros((len(idx),11));J[:,0]=sur[idx]==0;J[:,1]=sur[idx]==1;J[:,10]=h[idx]
        for xx,ic,iw,ia,ir in [(x1,2,4,6,8),(x2,3,5,7,9)]:
            w=p[iw]
            u=wrapped(xx[idx]-p[ic]-yr[idx]*p[ir])/w
            k=np.exp(-.5*u*u)
            if jac:
                dc=-p[ia]*k*u/w
                J[:,ic]=dc;J[:,ir]=dc*yr[idx];J[:,ia]=-k;J[:,iw]=-p[ia]*k*u*u/w
            z-=p[ia]*k
        return J if jac else z
    guesses=[]
    for j in range(6):
        phases=np.array(bases)+([0,0] if j<4 else [(-1)**j*.18,(-1)**(j+1)*.18])
        rates=([0,0],[.003,-.003],[-.010,.010],[.019,-.019],[0,0],[0,0])[j]
        p0=[1,1,*phases,.025,.025,.30,.30,*rates,0]
        res=least_squares(lambda p:(calc(p,ni)-y[ni])/e[ni],p0,jac=lambda p:calc(p,ni,True)/e[ni,None],
                          loss='soft_l1',f_scale=1.8,bounds=([.9,.9,-2,-2,.004,.004,0,0,-.065,-.065,-.06],
                                                              [1.1,1.1,2,2,.10,.10,.80,.80,.065,.065,.06]),
                          max_nfev=220,x_scale='jac',xtol=5e-7,ftol=5e-7)
        guesses.append(res)
    b=min(guesses,key=lambda r:r.cost);test=np.flatnonzero(~tr)
    out={'fold':fold,'first':la,'second':lb,'train_loss':2*b.cost,'params':b.x.tolist()}
    if len(test):
        r=(calc(b.x,test)-y[test])/e[test]
        out['test_loss']=2*1.8**2*np.sum(np.sqrt(1+(r/1.8)**2)-1)
    return out
if __name__=='__main__':
    jobs=[(f,a,b) for f in ['all',5,0,2,1,4,3] for a in ALIASES for b in ALIASES]
    with ProcessPoolExecutor(max_workers=3) as p:res=list(p.map(job,jobs,chunksize=3))
    (HERE/'aliases.json').write_text(json.dumps({'alias_frequency_offsets_c_per_day':ALIASES,'fits':res},indent=2)+'\n')
    for f in ['all',5,0,2,1,4,3]:
        rows=sorted([k for k in res if k['fold']==f],key=lambda k:k['train_loss'])
        print(f,[(k['first'],k['second'],round(k['train_loss']),round(k.get('test_loss',0))) for k in rows[:3]])
