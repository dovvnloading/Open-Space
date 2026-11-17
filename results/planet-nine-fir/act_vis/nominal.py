"""Two-frame differential-flux test at the exact tabulated ACT track."""
import csv,json
import numpy as np
from scipy.ndimage import map_coordinates,uniform_filter
from common import CACHE,HERE,MOTION,HALF,N
from search_patch import psf_field,sources,injected_psf
from dense import R,kernel,score_tile

def main():
    D=np.array([np.load(CACHE/f'data_{i}.npy') for i in range(4)])
    W=np.array([np.load(CACHE/f'weight_{i}.npy') for i in range(4)])
    coverage=json.loads((HERE/'coverage.json').read_bytes())
    pp,ns,ap=psf_field(D,W,coverage)
    runs=np.array([0,3]);bank=MOTION[runs][None]
    fs=D[runs];ws=W[runs];psf=pp[runs]
    counts=np.array([np.median([q['exp_e_per_vis25'] for q in coverage[i]['quadrants']]) for i in runs])
    z=score_tile(fs,ws,bank)[0]
    # Compare with an equally weighted sample of every fifteenth valid pixel cell.
    gy,gx=np.mgrid[30:N-30:15,30:N-30:15]
    good=np.all([uniform_filter((w>0).astype('f4'),size=7)[gy,gx]>.7 for w in ws],axis=0)
    dist=sources().query(np.column_stack((gx.ravel(),gy.ravel())))[0].reshape(gx.shape)
    good&=dist>15
    controls=z[gy[good],gx[good]];threshold=max(8.,float(np.percentile(controls,99.9)))
    sl=np.s_[HALF-R:HALF+R+1,HALF-R:HALF+R+1]
    ds=fs[:,sl[0],sl[1]];wt=ws[:,sl[0],sl[1]];dw=np.maximum(1e-20,wt.sum(axis=0))
    ks=np.array([np.outer(kernel(py),kernel(px)) for px,py in bank[0]])
    centered=ks-(wt*ks).sum(axis=0)/dw
    coeff=wt*centered/np.sqrt((wt*centered**2).sum())
    score=float((coeff*ds).sum())
    assert abs(float(z[HALF,HALF])-score)<.0001
    off=np.arange(-R,R+1);dx,dy=np.meshgrid(off,off)
    src=np.array([injected_psf(p,dx,dy,ox,oy) for (ox,oy),p in zip(bank[0],psf)])
    # A gray 5 M_earth object with 2.71 R_earth at 300 au and p=0.1.
    vis=-26.986100831142572-2.5*np.log10(.1*(2.71*6371/149597870.7)**2/(300**4))
    mags=np.array([vis,22.,23.,24.,24.5,25.,25.25,25.5,25.75,26.,26.5,27.])
    rng=np.random.default_rng(6742);res=[]
    for m in mags:
        sample=rng.poisson(np.maximum(0.,src*10**(.4*(25-m))*counts[:,None,None]),size=(300,*src.shape))
        zi=np.sum(coeff*(ds+sample/counts[:,None,None]),axis=(1,2,3))
        res.append(dict(VIS_AB=m,shot_noise_repeats=len(zi),passes=int((zi>=threshold).sum()),p05=np.percentile(zi,5),median=np.median(zi),p95=np.percentile(zi,95)))
    out=dict(exposure_indices=runs.tolist(),actual_score=score,control_count=len(controls),control_999=float(np.percentile(controls,99.9)),threshold=threshold,
        control_max=float(np.max(controls)),static_separation_arcsec=float(sources().query([HALF,HALF])[0]*.1),model_VIS_AB=vis,results=res)
    (HERE/'two_frame_nominal.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
