"""Validate spatial attribution and extrapolation using compressed TESS measurements."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from tess_timing import GROUPS, fit

HERE=Path(__file__).resolve().parent
ap='3x3'

def image_center(image,target):
    yy,xx=np.indices(image.shape)
    r=np.hypot(xx-target[0],yy-target[1]);ok=r<=2.7
    def estimate(p):
        a,x,y,lsx,lsy,level=p
        return level+a*np.exp(-.5*((xx-x)**2/np.exp(2*lsx)+(yy-y)**2/np.exp(2*lsy)))
    iy,ix=np.unravel_index(np.argmax(np.where(ok,image,-np.inf)),image.shape)
    fit=least_squares(lambda p:(estimate(p)-image)[ok],[max(image[iy,ix],.01),ix,iy,-.3,-.3,0],
          bounds=([0,target[0]-2.6,target[1]-2.6,-1.05,-1.05,-5],[20,target[0]+2.6,target[1]+2.6,.75,.75,5]),loss='soft_l1',f_scale=.05)
    return {'xy':fit.x[1:3].tolist(),'offset_pixel':float(np.hypot(fit.x[1]-target[0],fit.x[2]-target[1])),
            'xy_sigma':np.exp(fit.x[3:5]).tolist(),'amplitude':float(fit.x[0])}

def spatial():
    ph=(np.arange(120)+.5)/120
    c=np.array([-.0061,.5337])
    G=np.exp(-.5*(((ph[:,None]-c+.5)%1-.5)/.025)**2)
    meta=json.loads((HERE/'tess_inputs.json').read_text());out={}
    D=np.column_stack([np.ones(120),-G])
    def co(y,n):
        w=np.sqrt(n)[:,None];return np.linalg.lstsq(D*w,y*w,rcond=None)[0][1:,:].reshape(2,11,11)
    with np.load(HERE/'folded_pixels.npz') as data:
        for s in meta:
            n=data[f's{s}_counts'];y=data[f's{s}_images'].reshape(120,121)
            imgs=co(y,n);target=meta[s]['target_pixel'];locs=[image_center(im,target) for im in imgs]
            out[s]={'primary':locs[0],'secondary':locs[1],'combined':image_center(np.mean(imgs,axis=0),target)}
    return out

def contrast():
    rng=np.random.default_rng(99754)
    meta=json.loads((HERE/'tess_inputs.json').read_text());out={}
    with np.load(HERE.parent/'joint'/'inputs'/'orbital_folds.npz') as data:
        pool={}
        for s in [1,2,13,94,95,96]:
            n=data[f's{s}_counts'];ss=data[f's{s}_sums'][:,:,0];dt=meta[str(s)]['cadence_days']
            ok=(n.sum(axis=1)*dt>.65*1.1696719)&((n>0).sum(axis=1)>=min(55,int(.65*1.1696719/dt)))
            pool[s]=(ss[ok],n[ok])
        for label,sectors in [('2018_2019',[1,2,13]),('2025',[94,95,96])]:
            S,N=map(list,zip(*[pool[s] for s in sectors]))
            vals=[]
            # Nonoverlapping contiguous segments of up to three accepted orbits,
            # independently resampled within sectors, to allow short-run correlation.
            blocks=[[np.arange(j,min(j+3,len(a))) for j in range(0,len(a),3)] for a in S]
            for _ in range(650):
                inds=[np.concatenate([b[k] for k in rng.integers(0,len(b),len(b))]) for b in blocks]
                vals.append(fit([a[i] for a,i in zip(S,inds)],[n[i] for n,i in zip(N,inds)])[2])
            p=fit(S,N)[2];out[label]={'point':p,'draws':vals}
        v=out['2025']['point']-out['2018_2019']['point']
        d=np.array(out['2025']['draws'])-out['2018_2019']['draws']
        return {'point_delta_cycles':v,'ci95_delta_cycles':np.percentile(d,[2.5,97.5]).tolist(),
          'period_hours':24*1.1696719,'epochs':{k:{'separation':r['point'],'ci95':np.percentile(r['draws'],[2.5,97.5]).tolist()} for k,r in out.items()},'cluster_size':3,'bootstrap_draws':650}
if __name__=='__main__':
    maps=spatial();print([(s,round(maps[s]['combined']['offset_pixel'],3)) for s in maps],flush=True)
    diff=contrast();print(diff)
    (HERE/'independent_checks.json').write_text(json.dumps({'centroids':maps,'earliest_to_latest':diff},indent=2)+'\n')
