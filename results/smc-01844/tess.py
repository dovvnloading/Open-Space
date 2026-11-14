"""Recompute blind TESS peak rankings, fixed-phase tests and spatial difference images."""
import argparse
import json
from pathlib import Path
from scipy.optimize import least_squares
from astropy.timeseries import LombScargle
import numpy as np

from replay import EPOCH, BINARY_PERIOD, DATA, HERE

SECTORS = (1, 27, 28, 67, 68, 94, 95)
FREQUENCIES = np.arange(.1, 15, .0033)
FREQ = 4.912246653669684
# Unit Fourier vector for flux, opposite the magnitude vector.
UNIT = np.array([0.00946839790662613, -0.0026546456546182352])
UNIT /= np.linalg.norm(UNIT)


def clock_matrix(t, n):
    arg = 2*np.pi*(t-EPOCH)/BINARY_PERIOD
    return np.column_stack([v for k in range(1,n+1) for v in (np.cos(k*arg),np.sin(k*arg))])


def detrend_blocks(t, v):
    """Project each column away from a constant and slope in each 0.8-day block."""
    blocks = np.floor((t-t[0])/.8).astype(int)
    out=v.copy()
    for k in np.unique(blocks):
        ix=blocks==k
        x=np.column_stack((np.ones(sum(ix)),t[ix]-np.mean(t[ix])))
        out[ix]-=x@np.linalg.lstsq(x,v[ix],rcond=None)[0]
    return out,blocks


def scan_sector(t, flux):
    clocks=clock_matrix(t,8)
    sf=2*np.pi*(t-EPOCH)*FREQ
    terms=np.column_stack((flux,clocks,np.cos(sf),np.sin(sf)))
    values,blocks=detrend_blocks(t,terms)
    y=values[:,0]
    x=values[:,1:-2]
    res=y-x@np.linalg.lstsq(x,y,rcond=None)[0]
    limit=5*np.std(res)
    power=LombScargle(t,np.clip(res,-limit,limit)).power(FREQUENCIES)
    peak=int(np.argmax(power))
    x=values[:,1:]
    b=np.linalg.lstsq(x,y,rcond=None)[0]
    r=y-x@b
    scores=x*r[:,None]
    meat=np.stack([np.bincount(blocks,weights=c) for c in scores.T])
    inv=np.linalg.pinv(x.T@x)
    cov=inv@meat@meat.T@inv
    return dict(cadences=len(t),start_hjd_approx=float(t[0]+2450000),
                end_hjd_approx=float(t[-1]+2450000),blocks=len(np.unique(blocks)),
                top_frequency=float(FREQUENCIES[peak]),top_power=float(power[peak]),
                signed_electrons_per_s=float(UNIT@b[-2:]),
                block_se=float(np.sqrt(UNIT@cov[-2:,-2:]@UNIT)))


def centroid(image,xy):
    y,x=np.indices(image.shape)
    xx,yy=xy
    use=(x-xx)**2+(y-yy)**2<3.1**2
    def residual(v):
        amp,cx,cy,sx,sy,bg=v
        return amp*np.exp(-.5*(((x[use]-cx)/sx)**2+((y[use]-cy)/sy)**2))+bg-image[use]
    result=least_squares(residual,[max(.01,image.max()),xx,yy,.8,.8,0],
               bounds=([0,xx-2,yy-2,.35,.35,-np.inf],[np.inf,xx+2,yy+2,2,2,np.inf]))
    return result.x[1:3]


def maps(x,y):
    coef=np.linalg.lstsq(np.sum(x,axis=0),np.sum(y,axis=0),rcond=1.e-10)[0]
    h=np.arange(1,7)
    difference=np.ravel(np.column_stack((np.cos(h*np.pi/2)-1.,np.sin(h*np.pi/2))))
    return (UNIT@coef[-2:]).reshape(7,7),(difference@coef[:-2]).reshape(7,7)


def spatial(blocks, metadata, sector, nboot):
    x,y=blocks[f'x{sector}'],blocks[f'y{sector}']
    xy=metadata['x_y_target_zero_based']
    transform=np.array(metadata['pixel_scale_matrix_deg'])*3600
    def offset(x,y):
        short,orb=maps(x,y)
        c_short,c_orb=centroid(short,xy),centroid(orb,xy)
        # local WCS plane has axes east/north in degrees
        return transform@(c_short-c_orb),c_short,c_orb
    disp,c1,c2=offset(x,y)
    rng=np.random.default_rng(1844+int(sector))
    draws=[]
    for _ in range(nboot):
        idx=rng.integers(0,len(x),size=len(x))
        draws.append(offset(x[idx],y[idx])[0])
    return dict(used_blocks=len(x),short_xy=c1.tolist(),orbital_xy=c2.tolist(),
                east_north_arcsec=disp.tolist(),percentiles_2_5_50_97_5=np.percentile(draws,[2.5,50,97.5],axis=0).tolist())


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--draws',type=int,default=150)
    args=ap.parse_args()
    metadata=json.loads((HERE/'sources.json').read_text())['TESS']
    with np.load(DATA/'pixel-blocks.npz') as blocks:
        out={}
        for sector in SECTORS:
            with np.load(DATA/f's{sector}.npz') as l:
                out[sector]=scan_sector(l['t'],l['flux'])
            out[sector]['position']=spatial(blocks,metadata[str(sector)],sector,args.draws)
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
