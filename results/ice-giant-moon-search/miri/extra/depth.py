"""Off-grid native-image injection and neighboring-bank recovery for the Uranus MIRI arcs."""
import argparse,csv,math
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.stats import beta
from coverage import OUT,frame,get_states,pixels
from fetch import products
from native_search import maps,visit_scores
from miri_audit import write_csv


def get_xy(frames,rr,ph,origin):
    xy=np.full((len(rr),len(frames),2),np.nan)
    for j,g in enumerate(frames):
        i,x,y=pixels(g,rr,ph,origin);xy[i,j,:]=np.column_stack((x,y))
    return xy


def sample(frames,scoremaps,positions):
    r=np.zeros(positions.shape[:2],'f4');ok=np.zeros(r.shape,bool)
    for j,((scores,mask),g) in enumerate(zip(scoremaps,frames)):
        good=np.isfinite(positions[:,j,0]);xy=positions[good,j,::-1].T
        r[good,j]=ndi.map_coordinates(scores,xy,order=1,mode='constant')
        ok[good,j]=ndi.map_coordinates(mask.astype('f4'),xy,order=1,mode='constant')>.995
    return r,ok


def qualify(raw,support,frames):
    _,_,z,cover,vs=visit_scores(raw,support,frames,'Uranus')
    worst=np.full(z.shape,np.inf)
    for j,v in enumerate(vs):
        ids=np.flatnonzero([g['f']['visit']==v for g in frames])
        blocks=np.array([frames[k]['f']['file'].split('_')[1] for k in ids])
        for b in np.unique(blocks):
            subset=ids[blocks!=b];n=support[:,subset].sum(1)
            score=np.where(support[:,subset],np.clip(raw[:,subset],-7,7),0).sum(1)/np.sqrt(np.maximum(n,1))
            worst[:,j]=np.minimum(worst[:,j],np.where(n>=3,score,-99))
    return (np.minimum(z,worst/3.2)>=1).sum(1)>=2,(cover&(worst>-99)).sum(1)>=2


def run(cache,head):
    state=get_states();frames=[];original=[]
    for f,p in products():
        if f['planet']=='Uranus':
            g=frame(f,p,head,state[f['file']],{});g['product']=p;frames.append(g)
    frames.sort(key=lambda f:f['t']);origin=frames[0]['t']
    for g in frames:original.append(maps(g,g['product'],cache))
    rng=np.random.default_rng(1248362026)
    # Uniform continuous radius/phase draws conditioned on unplanted imaging support.
    radii=rng.uniform(400000,600000,5000);phase=rng.uniform(0,2*np.pi,len(radii))
    xy=get_xy(frames,radii,phase,origin);raw,sup=sample(frames,original,xy);prior,covered=qualify(raw,sup,frames)
    pool=list(np.flatnonzero(covered&~prior));plants=[];trials=[]
    for batch in range(2):
        chosen=[]
        for i in pool:
            if i in plants:continue
            if chosen:
                dif=np.linalg.norm(xy[chosen]-xy[i],axis=-1);finite=np.isfinite(dif)
                if np.any(finite&(dif<34.)):continue
            chosen.append(i)
            if len(chosen)==16:break
        if len(chosen)<12:raise RuntimeError(f'separated batch too small: {len(chosen)}')
        plants+=chosen
        bank=[];indices=[]
        for j,i in enumerate(chosen):
            a=round(radii[i]/1000)*1000;phi=round(np.rad2deg(phase[i])/.15)*.15
            for da in [-1000,0,1000]:
                for dp in [-.15,0,.15]:bank.append((a+da,np.deg2rad((phi+dp)%360)));indices.append(j)
        bank=np.array(bank);indices=np.array(indices);grid_xy=get_xy(frames,bank[:,0],bank[:,1],origin)
        target=xy[np.array(chosen)[indices]]
        delta=np.linalg.norm(target-grid_xy,axis=2);match=np.sqrt(np.nanmean(delta**2,axis=1))<2
        unplanted,oks=sample(frames,original,grid_xy);initial=qualify(unplanted,oks,frames)[0]
        for flux in [0.1,0.4,1.6,6.4,25.6]:
            planted=[]
            for j,g in enumerate(frames):
                ij=[(*xy[i,j],flux) for i in chosen if np.isfinite(xy[i,j]).all()]
                planted.append(maps(g,g['product'],cache,ij))
            raw,mask=sample(frames,planted,grid_xy);passing=qualify(raw,mask,frames)[0]
            for j,i in enumerate(chosen):
                success=np.any((indices==j)&match&~initial&passing)
                trials.append(dict(batch=batch,draw=i,a_km=radii[i],phase_deg=np.rad2deg(phase[i]),flat_fnu_ujy=flux,AB=23.9-2.5*np.log10(flux),recovered=int(success)))
            print(batch,flux,sum(r['recovered'] for r in trials if r['batch']==batch and r['flat_fnu_ujy']==flux),flush=True)
            del planted
    levels=[]
    for flux in sorted({r['flat_fnu_ujy'] for r in trials}):
        a=[r for r in trials if r['flat_fnu_ujy']==flux];k=sum(r['recovered'] for r in a);n=len(a)
        levels.append(dict(flat_fnu_ujy=flux,AB=23.9-2.5*np.log10(flux),recovered=k,plants=n,
                           binomial_95_lo=0 if k==0 else beta.ppf(.025,k,n-k+1),binomial_95_hi=1 if k==n else beta.ppf(.975,k+1,n-k)))
    write_csv(OUT/'uranus_plants.csv',list(trials[0]),trials);write_csv(OUT/'uranus_depth.csv',list(levels[0]),levels)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('headers',type=Path);a=p.parse_args();run(a.cache,a.headers)
