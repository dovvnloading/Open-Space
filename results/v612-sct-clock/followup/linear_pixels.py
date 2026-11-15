"""Independent *linear* reduction of original TESS FFIs, not the PR #24 PCA pipeline.

Per exposure: 4 aperture pixels minus 4 times the mean of eleven outer pixels.
Per UTC day: fit nuisance intercept and slope to BOTH data and design columns.
Then fit frozen f/2, f, 2f and 3f; resample whole days and flip day signs for nulls.
Downloaded raw cutouts are tracked by SHA256 and observation FFI names.
Requires numpy, pandas, scipy, astropy.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.special import ndtr

F=1/.15956200154633035
EPOCH=2459499.9473282346
RA,DEC=277.941193,-14.315539
RNG=np.random.default_rng(26449)


def basis(t):
    angle=2*np.pi*(t+2457000-EPOCH)
    return np.column_stack([x for f in [F/2,F,2*F,3*F]
                            for x in [np.cos(angle*f),np.sin(angle*f)]])


def daily_project(t,Y,X):
    # Independent, linear, daily polynomial nuisance; no nonlinear median and no PCA.
    day=np.floor(t).astype(int)
    Y=Y.copy()
    X=X.copy()
    for d in np.unique(day):
        m=day==d
        u=t[m]-t[m].mean()
        D=np.column_stack((np.ones(sum(m)),u))
        Y[m]-=D@np.linalg.lstsq(D,Y[m],rcond=None)[0]
        X[m]-=D@np.linalg.lstsq(D,X[m],rcond=None)[0]
    return Y,X


def angle(b):
    return float((np.arctan2(b[1],b[0])/(2*np.pi)+.5)%1-.5)


def bootstrap(X,y,t,sector,repeats=2500):
    groups=np.floor(t).astype(int)
    unique=np.unique(groups)
    grams=np.array([X[groups==g].T@X[groups==g] for g in unique])
    vectors=np.array([X[groups==g].T@y[groups==g] for g in unique])
    rng=np.random.default_rng(510+sector)
    index=rng.integers(len(unique),size=(repeats,len(unique)))
    return np.linalg.solve(grams[index].sum(axis=1),
                           vectors[index].sum(axis=1)[...,None])[...,0]


def block_sign_test(t,y,X,pair,reference=None,reps=10000):
    others=[j for j in range(X.shape[1]) if j not in pair]
    nullX=X[:,others]
    Z=X[:,pair]
    residual=y-nullX@np.linalg.lstsq(nullX,y,rcond=None)[0]
    projected=Z-nullX@np.linalg.lstsq(nullX,Z,rcond=None)[0]
    inv=np.linalg.inv(projected.T@projected)
    groups=np.floor(t).astype(int)
    score=np.array([projected[groups==d].T@residual[groups==d]
                    for d in np.unique(groups)])
    original=score.sum(axis=0)
    coef=inv@original
    flipped=RNG.choice([-1,1],size=(reps,len(score)))@score
    if reference is None:
        def stat(z):return np.sqrt(np.sum((z@inv)*z,axis=-1))
        value=stat(original)
        simulated=stat(flipped)
    else:
        u=np.array([np.cos(2*np.pi*reference),np.sin(2*np.pi*reference)])
        value=float(coef@u)
        simulated=(flipped@inv)@u
    k=int(np.sum(simulated>=value))
    return dict(amp_e_per_s=float(np.linalg.norm(coef)),phase_cycles=angle(coef),
                signflip_exceedances=k,signflip_trials=reps,
                p=(k+1)/(reps+1),test_statistic=float(value))


def gaussian(p):
    x,y,sigma,amplitude,bg=p
    grid=np.arange(5)
    vx=ndtr((grid+.5-x)/sigma)-ndtr((grid-.5-x)/sigma)
    vy=ndtr((grid+.5-y)/sigma)-ndtr((grid-.5-y)/sigma)
    return amplitude*np.outer(vy,vx).ravel()+bg


def fit_gaussian(image,sky):
    # NOT a calibrated TESS PRF. Do not interpret these errors as absolute astrometry.
    y=image.ravel()
    init=[float(sky[0]),float(sky[1]),.75,
          float(np.sum(y-np.median(y))),float(np.median(y))]
    result=least_squares(lambda p:gaussian(p)-y,init,
                         bounds=([sky[0]-.7,sky[1]-.7,.35,-100,-100],
                                 [sky[0]+.7,sky[1]+.7,1.4,100,100]),max_nfev=45)
    return result.x


def study(sector,path,spatial_boot=150):
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    with fits.open(path,memmap=True) as fits_file:
        rows=fits_file[1].data
        xy=WCS(fits_file[2].header).all_world2pix([[RA,DEC]],0)[0]
        good=((rows["QUALITY"]==0)&np.isfinite(rows["TIME"])&
              np.isfinite(rows["FLUX"]).all(axis=(1,2)))
        t=np.asarray(rows["TIME"][good],dtype=float)
        pix=np.asarray(rows["FLUX"][good],dtype=float).reshape((-1,25))
        first_ffi=str(rows["FFI_FILE"][good][0])
        last_ffi=str(rows["FFI_FILE"][good][-1])
    yy,xx=np.indices((5,5))
    r=np.hypot(xx-xy[0],yy-xy[1]).ravel()
    ap=r<1.15
    bg=r>=2.1
    y=pix[:,ap].sum(axis=1)-int(ap.sum())*pix[:,bg].mean(axis=1)
    A=basis(t)
    yres,X=daily_project(t,y,A)
    beta=np.linalg.lstsq(X,yres,rcond=None)[0]
    draws=bootstrap(X,yres,t,sector)
    row=dict(sector=sector,raw_fits_sha256=digest,
             first_ffi=first_ffi,last_ffi=last_ffi,
             n_good=int(len(t)),n_days=int(len(np.unique(np.floor(t)))),
             gaia_xy=xy.tolist(),aperture_xy=np.flatnonzero(ap).tolist(),
             outer_background_xy=np.flatnonzero(bg).tolist(),
             extraction="per-cadence mean of 11 outer pixels; fit per-day intercept+slope; no median or PCA")
    for i,key in enumerate(("f/2","f","2f","3f")):
        v=beta[2*i:2*i+2]
        z=draws[:,2*i:2*i+2]
        ph=angle(v)
        phases=(np.arctan2(z[:,1],z[:,0])/(2*np.pi)-ph+.5)%1-.5
        row[key]=dict(amp_e_per_s=float(np.linalg.norm(v)),
            phase_cycles=ph,
            bootstrap_amp_95=np.quantile(np.linalg.norm(z,axis=1),[.025,.975]).tolist(),
            bootstrap_phase_shift_95=np.quantile(phases,[.025,.975]).tolist())
    row["2f_null"]=block_sign_test(t,yres,X,[4,5])
    row["half_frequency_null"]=block_sign_test(t,yres,X,[0,1])
    # In a crowded field, only amplitude maps, not this Gaussian model, are
    # direct image measurements. A symmetric Gaussian may be biased by real PRF.
    raw_res=np.empty_like(pix)
    for j in range(25):
        raw_res[:,j],_=daily_project(t,pix[:,j],A)
    per_pixel=np.linalg.lstsq(X,raw_res,rcond=None)[0]
    maps={"f":per_pixel[2].reshape(5,5),
          "2f":-per_pixel[4].reshape(5,5)}
    location={}
    for key,mp in maps.items():
        g=fit_gaussian(mp,xy)
        location[key]=dict(gaussian_centroid_xy=g[:2].tolist(),
                           offset_xy=(g[:2]-xy).tolist(),
                           sigma_pixel=float(g[2]),
                           signed_cos_map_e_per_s=mp.round(6).tolist())
    row["gaussian_spatial_diagnostic"]=location
    return row,t,yres,X


def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--fits-dir",type=Path,required=True)
    arg.add_argument("--output",type=Path)
    arg.add_argument("--check",action="store_true")
    args=arg.parse_args()
    o={}
    first=None
    for s in (80,92):
        path=args.fits_dir/f"tess{s}.fits"
        r,t,y,X=study(s,path)
        if first is None:first=r["2f"]["phase_cycles"]
        else:
            r["2f_sector80_phase_locked_null"]=block_sign_test(t,y,X,[4,5],reference=first)
            r["2f_phase_separation_cycles"]=float((r["2f"]["phase_cycles"]-first+.5)%1-.5)
        o[str(s)]=r
    if args.check:
        assert .42<o["80"]["2f"]["amp_e_per_s"]<.9
        assert .32<o["92"]["2f"]["amp_e_per_s"]<.8
        assert abs(o["92"]["2f_phase_separation_cycles"])<.10
        assert o["92"]["2f_sector80_phase_locked_null"]["p"]<.01
        print("Pixel-level harmonic replay and phase-locked holdout passed.")
    message=json.dumps(o,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(exist_ok=True,parents=True)
        args.output.write_text(message)
    else:print(message)


if __name__=="__main__":
    main()
