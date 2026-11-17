"""Supplementary two-exposure sensitivity on cells excluded by the repeatability test."""
import csv,json
from itertools import combinations
import numpy as np
from scipy.ndimage import map_coordinates,maximum_filter,uniform_filter
from common import CACHE,HERE,MOTION,HALF,N,world
from search_patch import MAG,SEED,sources,psf_field,wilson,injected_psf
from dense import R,kernel,score_tile

def main():
    coverage=json.loads((HERE/'coverage.json').read_bytes())
    D=np.array([np.load(CACHE/f'data_{i}.npy') for i in range(4)])
    W=np.array([np.load(CACHE/f'weight_{i}.npy') for i in range(4)])
    P,_,_=psf_field(D,W,coverage)
    counts=np.array([np.median([q['exp_e_per_vis25'] for q in v['quadrants']]) for v in coverage])
    gy,gx=np.indices((N,N))
    patch=(abs(gx-HALF)<=600)&(abs(gy-HALF)<=600)
    cat=sources()
    static=cat.query(np.column_stack((gx.ravel(),gy.ravel())))[0].reshape(gx.shape)<=15
    local=np.array([uniform_filter((w>0).astype('f4'),size=7)>.7 for w in W])
    draws=np.random.default_rng(SEED).uniform(-600,600,(2000,2))+HALF
    xx,yy=np.rint(draws.T).astype(int);out=np.zeros((len(draws),len(MAG)),bool)
    total=0;per_pair=[];peaks=[];n_eligible=0
    grid=np.arange(-R,R+1);dx,dy=np.meshgrid(grid,grid)
    for i,j in combinations(range(4),2):
        pair=[i,j];m=MOTION[pair]
        sel=patch&(local.sum(axis=0)==2)&local[pair].all(axis=0)&~static
        field=score_tile(D[pair],W[pair],m[None])[0]
        controls=field[::15,::15][sel[::15,::15]]
        if not len(controls):continue
        threshold=max(8.,float(np.percentile(controls,99.9)))
        positive=(field==maximum_filter(field,size=25))&(field>=threshold)&sel
        for y,x in zip(*np.nonzero(positive)):
            ra,de=world(x-HALF,y-HALF)
            peaks.append(dict(exposures=pair,ra=ra,dec=de,score=float(field[y,x])))
        # Reserve all 1.2-arcsec neighbourhoods of strong positive or negative residuals.
        ok=sel&(maximum_filter(abs(field),size=25)<threshold)
        ids=np.flatnonzero(ok[yy,xx]);n_eligible+=len(ids)
        per_pair.append(dict(exposures=pair,displacement_arcsec=float(np.linalg.norm(m[1]-m[0])*.1),
                             possible_cells=int(sel.sum()),searched_cells=int(ok.sum()),draws=len(ids),
                             tail_threshold=threshold,controls=len(controls),positive_peaks=int(positive.sum())))
        total+=int(ok.sum())
        K=np.stack([np.outer(kernel(oy),kernel(ox)) for ox,oy in m])
        rng=np.random.default_rng(54000+10*i+j)
        for t in range(0,len(ids),40):
            k=ids[t:t+40]
            region=[yy[k,None,None]+dy,xx[k,None,None]+dx]
            d=D[pair][:,region[0],region[1]].transpose(1,0,2,3)
            w=W[pair][:,region[0],region[1]].transpose(1,0,2,3)
            sw=np.maximum(1e-20,w.sum(axis=1))
            centered=K[None]-(w*K[None]).sum(axis=1)[:,None]/sw[:,None]
            coeff=(w*centered)/np.sqrt(np.maximum(1e-20,np.sum(w*centered**2,axis=(1,2,3))))[:,None,None,None]
            src=[]
            for h,q in enumerate(pair):
                ox,oy=m[h];x0=draws[k,0]-xx[k]+ox;y0=draws[k,1]-yy[k]+oy
                src.append(injected_psf(P[q],dx[None],dy[None],x0[:,None,None],y0[:,None,None]))
            src=np.array(src).transpose(1,0,2,3)
            ne=counts[pair][None,:,None,None]
            for p,mag in enumerate(MAG):
                fake=rng.poisson(np.maximum(0.,src*10**(.4*(25-mag))*ne))
                z=np.sum(coeff*(d+fake/ne),axis=(1,2,3))
                out[k,p]=z>=threshold
    (HERE/'two_frame_patch.json').write_text(json.dumps(dict(pairs=per_pair,positive_peaks=peaks,
          eligible_cells=total,eligible_draws=n_eligible),indent=2)+'\n')
    base=list(csv.DictReader((HERE/'injections.csv').open()))
    with (HERE/'supplementary_injections.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['VIS_AB','all_draws','three_or_four_frames_recovered','two_frames_recovered','either_recovered','either_fraction','wilson95_lo','wilson95_hi'])
        for idx,m in enumerate(MAG):
            num=int(out[:,idx].sum());m0=int(base[idx]['recovered_draws']);tot=m0+num
            w.writerow([m,len(draws),m0,num,tot,tot/len(draws),*[float(np.clip(q,0,1)) for q in wilson(tot,len(draws))]])
    print('pairs',per_pair,'peaks',peaks)

if __name__=='__main__':main()
