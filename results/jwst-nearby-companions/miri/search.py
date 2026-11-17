"""Two-epoch, host-centered F1500W peak search with pre-filter image injections."""
import argparse
import csv
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from fetch_inputs import HERE, digest

STEP = 0.11
HALF = 320
INNER, OUTER = 4., 12.
MATCH, MARGIN = .20, .15  # arcsec
SIG = 6.
RATIO = (.001, .0025, .005, .01, .02, .04, .08, .16)
RADII = (5., 8., 11.)
PHASES = 40
warnings.simplefilter("ignore", FITSFixedWarning)


def smooth_response(image):
    return ndi.gaussian_filter(image, 1.7) - ndi.gaussian_filter(image, 5.5)


def refine(image, x, y):
    """Subpixel stationary point of a local quadratic, or the integer peak."""
    z = image[y-1:y+2, x-1:x+2]
    dx, dy = (z[1,2]-z[1,0])/2, (z[2,1]-z[0,1])/2
    xy = (z[2,2]-z[2,0]-z[0,2]+z[0,0])/4
    hess = np.array([[z[1,2]+z[1,0]-2*z[1,1], xy],
                     [xy, z[2,1]+z[0,1]-2*z[1,1]]])
    if np.max(np.linalg.eigvalsh(hess)) >= 0:
        return np.array([x,y], dtype=float)
    delta = -np.linalg.solve(hess,[dx,dy])
    return np.array([x,y])+np.clip(delta,-.8,.8)


def sky_grid(ra, dec):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.crval = [ra,dec]
    w.wcs.crpix = [HALF+1,HALF+1]
    w.wcs.cdelt = [STEP/3600,STEP/3600]
    return w


def sample(image, weight, native, grid):
    yy,xx = np.mgrid[-HALF:HALF+1,-HALF:HALF+1]+HALF
    world = grid.pixel_to_world_values(xx, yy)
    sx,sy = native.world_to_pixel_values(*world)
    pix = ndi.map_coordinates(image,[sy,sx],order=1,mode="constant",cval=np.nan)
    wt = ndi.map_coordinates(weight,[sy,sx],order=1,mode="constant",cval=0)
    return pix, wt


def load(path):
    with fits.open(path) as hd:
        p=hd[0].header
        h=hd["SCI"].header
        if h["BUNIT"]!="MJy/sr" or p["FILTER"]!="F1500W":
            raise ValueError(path)
        w = WCS(h, naxis=2)
        img,wt=hd["SCI"].data.astype(float), hd["WHT"].data.astype(float)
        x,y=w.world_to_pixel_values(p["TARG_RA"],p["TARG_DEC"])
        ix,iy = int(round(float(x))),int(round(float(y)))
        z=img[iy-8:iy+9,ix-8:ix+9]
        py,px = np.unravel_index(np.nanargmax(z),z.shape)
        x,y = refine(ndi.gaussian_filter(z,.8),px,py)+[ix-8,iy-8]
        ra, dec = w.pixel_to_world_values(x,y)
        grid=sky_grid(ra,dec)
        pixels,weights=sample(img,wt,w,grid)
        valid = np.isfinite(pixels) & (weights>.3*np.nanmedian(weights[weights>0]))
        good = ndi.minimum_filter(valid, size=31)
        pixels[~valid]=np.nanmedian(pixels[valid])
        filt = smooth_response(pixels)
        yy,xx = np.indices(pixels.shape)
        rad=np.hypot(xx-HALF,yy-HALF)*STEP
        # Independent local scatter for each field/epoch and three radial bins.
        edges = np.array([0,6.5,9.5,14,40])
        bins=np.digitize(rad,edges)-1
        med,noise = [],[]
        for i in range(len(edges)-1):
            mask=(bins==i)&good&(rad>INNER)
            values=filt[mask]
            m=np.median(values)
            med.append(m)
            noise.append(1.4826*np.median(abs(values-m)))
        idx = np.clip(bins,0,len(noise)-1)
        score=(filt-np.array(med)[idx])/np.array(noise)[idx]
        # Empirical host core. The taper suppresses nearby background sources.
        half=10
        stamp=pixels[HALF-half:HALF+half+1,HALF-half:HALF+half+1].copy()
        srad=rad[HALF-half:HALF+half+1,HALF-half:HALF+half+1]
        stamp-=np.median(stamp[(srad>.8)&(srad<1.1)])
        stamp*=np.clip((.9-srad)/.25,0,1)
        meta=dict(file=path.name,ra=float(ra),dec=float(dec),date=p["DATE-OBS"],
                  pipeline=p["CAL_VER"],crds=p.get("CRDS_CTX"),centroid_native_xy=[float(x),float(y)],
                  response_median=med,response_mad_sigma=noise,coverage_fraction=float(good[(rad>=INNER)&(rad<=OUTER)].mean()))
        return dict(score=score,noise=np.array(noise)[idx],good=good,rad=rad,grid=grid,psf=stamp,meta=meta)


def peaks(score, good, threshold=SIG):
    loc=(score==ndi.maximum_filter(score,size=5))&(score>=threshold)&good
    loc[:2]=loc[-2:]=False
    loc[:,:2]=loc[:,-2:]=False
    ys,xs=np.where(loc)
    positions=np.array([refine(score,x,y) for x,y in zip(xs,ys)]).reshape(-1,2)
    return positions,score[ys,xs]


def classify(p, second, sky_shift):
    """Nearest counterpart under each hypothesis; all distances in arcseconds."""
    db=np.linalg.norm(second-(p+sky_shift),axis=1)*STEP
    dc=np.linalg.norm(second-p,axis=1)*STEP
    if len(second)==0:
        return "unmatched",np.nan,np.nan,None,None
    ib,ic=int(np.argmin(db)),int(np.argmin(dc))
    b,c=db[ib],dc[ic]
    if c<=MATCH and b-c>=MARGIN:
        status="common-motion trigger"
    elif b<=MATCH and c-b>=MARGIN:
        status="stationary match"
    elif min(b,c)<=MATCH:
        status="ambiguous"
    else:
        status="unmatched"
    return status,b,c,ib,ic


def background_track(a,b):
    catalogs=[peaks(v["score"],v["good"]&(v["rad"]>14)&(v["rad"]<30),8.)[0] for v in (a,b)]
    p,q=catalogs
    # Project absolute coordinates into the second host frame for initial pairing.
    px,py=b["grid"].world_to_pixel_values(*a["grid"].pixel_to_world_values(p[:,0],p[:,1]))
    hint=np.transpose([px,py])
    d,idx=cKDTree(q).query(hint)
    keep=np.isfinite(d)&(d*STEP<.55)
    delta=q[idx]-p
    center=np.median(delta[keep],axis=0)
    resid=np.linalg.norm(delta-center,axis=1)*STEP
    keep &= resid < .18
    if keep.sum()<4:
        raise ValueError("Too few stationary tie sources")
    shift=np.median(delta[keep],axis=0)
    return shift,dict(n_ties=int(keep.sum()),offset_east_north_arcsec=(shift*STEP).tolist(),
                      radial_rms_arcsec=float(np.sqrt(np.mean(np.sum(((delta[keep]-shift)*STEP)**2,axis=1)))))


def injected_positions(obs, pos, contrast):
    # Local calculation of the linear filter is identical to adding a source
    # before filtering the full mosaic. Pad by the Gaussian's full support.
    x,y=np.rint(pos).astype(int)
    half=36
    base=obs["score"][y-half:y+half+1,x-half:x+half+1]
    sig=obs["noise"][y-half:y+half+1,x-half:x+half+1]
    valid=obs["good"][y-half:y+half+1,x-half:x+half+1]
    arr=np.zeros_like(base)
    stamp=obs["psf"]; h=stamp.shape[0]//2
    arr[half-h:half+h+1,half-h:half+h+1]=stamp
    arr=ndi.shift(arr,(pos[1]-y,pos[0]-x),order=1,mode="constant",prefilter=False)
    response=smooth_response(arr)/sig
    # Include all neighbors needed to compete under the stationary hypothesis.
    area=np.zeros_like(valid); area[half-17:half+18,half-17:half+18]=True
    for ratio in contrast:
        loc,_=peaks(base+ratio*response,valid&area)
        yield loc+[x-half,y-half]


def simulate(a,b,shift,name):
    rows=[]
    for radius in RADII:
        for angle in np.arange(PHASES)*2*np.pi/PHASES+.037*radius:
            pos=HALF+radius/STEP*np.array([np.cos(angle),np.sin(angle)])
            had=False
            for c,p,q in zip((0.,)+RATIO,injected_positions(a,pos,(0.,)+RATIO),injected_positions(b,pos,(0.,)+RATIO)):
                near=p[np.linalg.norm(p-pos,axis=1)*STEP<=MATCH]
                found=False
                for pt in near:
                    status,_,_,_,ic=classify(pt,q,shift)
                    if status=="common-motion trigger" and np.linalg.norm(q[ic]-pos)*STEP<=MATCH:
                        found=True
                if c==0:
                    had=found
                    continue
                rows.append(dict(target=name,radius_arcsec=radius,pa_deg_east_of_north=float((90-np.rad2deg(angle))%360),
                                 companion_host_ratio=c,baseline_match=int(had),recovered=int(found and not had)))
    return rows


def write_csv(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    par=argparse.ArgumentParser(description=__doc__);par.add_argument("cache",type=Path);par.add_argument("out",type=Path)
    args=par.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((HERE/"inputs.json").read_text())["inputs"]
    info={};sources=[];trials=[]
    for target,obs in [("WD1202-232","o002"),("WD2105-82","o003")]:
        inputs=sorted((r for r in manifest if "-"+obs+"_" in r["obs_id"]),key=lambda r:r["t_min"])
        if len(inputs)!=2:raise ValueError(inputs)
        for r in inputs:
            if digest(args.cache/r["filename"])!=r["sha256"]:raise ValueError(r["filename"])
        a,b=[load(args.cache/r["filename"]) for r in inputs]
        shift,reg=background_track(a,b)
        p,s1=peaks(a["score"],a["good"]&(a["rad"]>=INNER)&(a["rad"]<=OUTER))
        q,s2=peaks(b["score"],b["good"]&(b["rad"]>INNER-2)&(b["rad"]<OUTER+2))
        info[target]=dict(epochs=[o["meta"] for o in (a,b)],registration=reg,early_triggers=len(p))
        for i,pt in enumerate(p,1):
            status,db,dc,ib,ic=classify(pt,q,shift)
            ra,dec=a["grid"].pixel_to_world_values(*pt)
            east,north=(pt-HALF)*STEP
            row=dict(target=target,id=i,ra_early_deg=float(ra),dec_early_deg=float(dec),east_arcsec=east,
                     north_arcsec=north,separation_arcsec=np.hypot(east,north),early_score=s1[i-1],status=status,
                     stationary_residual_arcsec=db,comoving_residual_arcsec=dc)
            for label,index in (("static",ib),("moving",ic)):
                if index is not None:
                    offsets=(q[index]-HALF)*STEP
                    row.update({label+"_late_east_arcsec":offsets[0],label+"_late_north_arcsec":offsets[1],label+"_late_score":s2[index]})
            sources.append(row)
        trials.extend(simulate(a,b,shift,target))
    write_csv(args.out/"sources.csv",sources)
    write_csv(args.out/"injections.csv",trials)
    curves=[]
    for target in info:
        for r in RADII:
            for c in RATIO:
                part=[v for v in trials if (v["target"],v["radius_arcsec"],v["companion_host_ratio"])==(target,r,c)]
                n=sum(v["recovered"] for v in part)
                curves.append(dict(target=target,radius_arcsec=r,companion_host_ratio=c,n_injected=len(part),n_recovered=n,fraction=n/len(part)))
    write_csv(args.out/"contrast.csv",curves)
    (args.out/"measurement.json").write_text(json.dumps(info,indent=2)+"\n")

if __name__=="__main__":
    main()
