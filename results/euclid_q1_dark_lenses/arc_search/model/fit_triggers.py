"""Compare a one-source SIS+shear with independent emitting knots on held-out native VIS blocks."""
import argparse
import csv
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from urllib.request import urlopen

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
from scipy.ndimage import binary_dilation, gaussian_filter
from scipy.optimize import least_squares, nnls

HERE = Path(__file__).resolve().parent
ARC = HERE.parent
sys.path.insert(0,str(ARC))
from photometry import pix, read_primary

# Common native pixel grid; convolve a circular 0.17" VIS PSF after 2x integration.
PSF = 0.17/2.355/0.1


def targets():
    selected = [r for f in sorted((ARC/"profiles").glob("*.json"))
                for r in json.loads(f.read_text())["profiles"] if r["pass_vis"]]
    extra = json.loads((ARC/"unmodelled.json").read_text())
    for obj in extra:
        for j, p in enumerate(obj.get("positions",[])):
            if p["pass_vis"]:
                selected.append({"id":obj["id_str"].split("_")[0]+"_target_"+obj["id_str"]+"_"+str(j),
                                 "ra":p["ra"],"dec":p["dec"],"radius_arcsec":p["r"]})
    with (HERE/"vis_inputs.csv").open() as stream:
        cohort={r["id"] for r in csv.DictReader(stream)}
    result={r["id"]:r for r in selected if r["id"] in cohort}
    if set(result)!=cohort:
        raise ValueError("Search centre missing for a frozen follow-up entry")
    return result


def load_target(t, cache):
    with (HERE/"vis_inputs.csv").open() as stream:
        files = [r for r in csv.DictReader(stream) if r["id"]==t["id"] and "/VIS/" in r["product"]]
    maps={}
    for f in files:
        filename=f["product"].split("/")[-1]
        path=cache/(t["id"]+"_"+filename)
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()!=f["sha256"]:
            # Earlier validation can cache a 24" stamp under the same basename used for a 40" stamp.
            path=cache/(f["sha256"][:16]+"_"+filename)
        if not path.exists():
            with urlopen(f["url"],timeout=100) as response:
                path.write_bytes(response.read())
        raw=path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=f["sha256"]:
            raise ValueError("Bad SHA256: "+str(path))
        h,arr=read_primary(raw)
        key="rms" if "-RMS_" in filename else "flag" if "-FLAG_" in filename else "image"
        maps[key]=(h,arr)
    with (HERE/"initial_knots.csv").open() as stream:
        knots=np.array([[float(r["dx_arcsec"]),float(r["dy_arcsec"])] for r in csv.DictReader(stream)
                         if r["id"]==t["id"] and r["band"]=="VIS" and r["method"]=="matched_035" and r["point"]!="center"])
    return maps,knots


def gaussian(x,y,s,q,angle):
    ca,sa=np.cos(angle),np.sin(angle)
    return np.exp(-0.5*(((ca*x+sa*y)/s)**2+((-sa*x+ca*y)/(s*q))**2))


def fit_one(job):
    t,cache,out,fold,center_bound=job
    maps,knots=load_target(t,cache)
    h,im=maps["image"]
    cx,cy=pix(h,t["ra"],t["dec"])
    R=t["radius_arcsec"]
    radius=max(1.65*R,R+0.55)
    padding=int(np.ceil(radius*10))+5
    xx,yy=int(round(cx)),int(round(cy))
    crop=np.s_[yy-padding:yy+padding+1,xx-padding:xx+padding+1]
    data=im[crop].astype(float); rms=maps["rms"][1][crop].astype(float); flags=maps["flag"][1][crop]
    yp,xp=np.indices(data.shape)
    x=(xp-padding+xx-cx)*0.1; y=(yp-padding+yy-cy)*0.1
    aperture=np.hypot(x,y)<=radius
    valid=aperture&np.isfinite(data+rms)&(rms>0)&(rms<1e8)
    # Q1 flag OR may be nonzero for valid coadds; log fractions and use finite positive native RMS.
    med=np.median(rms[valid])
    values=data/med; errors=rms/med
    blocksize=6
    bx=xp//blocksize; by=yp//blocksize
    test=valid & ((bx+2*by)%5==fold)
    train=valid & ~binary_dilation(test,iterations=1)
    # Keep the trial centre and native orientation for all positions. Average four subpixels.
    X=x[None,:,:]+np.array([-.025,-.025,.025,.025])[:,None,None]
    Y=y[None,:,:]+np.array([-.025,.025,-.025,.025])[:,None,None]

    central=[gaussian_filter(np.exp(-(x*x+y*y)/(2*s*s)),PSF) for s in [.14,.45]]
    def profiles(p,kind):
        if kind=="lens":
            b,dx,dy,g1,g2,sx,sy,s,q,ph=p
            tx,ty=X-dx,Y-dy
            rad=np.sqrt(tx*tx+ty*ty+1e-12)
            ux=tx*(1-b/rad-g1)-g2*ty-sx
            uy=ty*(1-b/rad+g1)-g2*tx-sy
            comps=[gaussian_filter(gaussian(ux,uy,s,q,ph).mean(axis=0),PSF)]
        else:
            comps=[]
            for dx,dy,s,q,ph in np.asarray(p).reshape(-1,5):
                comps.append(gaussian_filter(gaussian(X-dx,Y-dy,s,q,ph).mean(axis=0),PSF))
        return [*comps,*central]
    def solve(p,kind,full=False):
        comps=profiles(p,kind)
        A=np.column_stack([z[train]/errors[train] for z in comps]); sky=1/errors[train]
        obs=values[train]/errors[train]
        # Project out a free constant, then solve positive amplitudes, keeping both models comparable.
        hsky=sky@sky
        Ap=A-np.outer(sky,sky@A/hsky)
        bp=obs-sky*(sky@obs/hsky)
        weights=nnls(Ap,bp)[0]
        level=sky@(obs-A@weights)/hsky
        if full:
            return level+sum(w*z for w,z in zip(weights,comps)),weights,level
        return A@weights+sky*level-obs
    b_bounds=(max(.17,.55*R),min(5,1.55*R))
    lo=np.array([b_bounds[0],-center_bound,-center_bound,-.32,-.32,-R,-R,.035,.30,-np.pi])
    hi=np.array([b_bounds[1],center_bound,center_bound,.32,.32,R,R,.70,1.,np.pi])
    sols=[]
    for i,(f,alpha) in enumerate([(.83,0),(1.,.5),(1.16,1.1),(1.,1.8)]):
        ux,uy=knots[i%len(knots)]
        vec=np.array([ux,uy]) / np.hypot(ux,uy)
        start=np.array([f*R,0,0,.13*np.cos(2*alpha),.13*np.sin(2*alpha),
                        *(vec * (.18*R)),.14,.72,alpha])
        s=least_squares(solve,np.clip(start,lo+1e-4,hi-1e-4),args=("lens",),
                        bounds=(lo,hi),max_nfev=90,ftol=.00015,xtol=.00015)
        sols.append(s)
    fitted=min(sols,key=lambda fit:fit.cost)
    # Independent positions may move +/-0.45" from the observed peak. Widths are intrinsic.
    init=np.array([[a,b,.18,.75,0] for a,b in knots]).flatten()
    lower=np.array([[a-.45,b-.45,.035,.30,-np.pi] for a,b in knots]).flatten()
    upper=np.array([[a+.45,b+.45,.80,1.,np.pi] for a,b in knots]).flatten()
    free=least_squares(solve,init,args=("knots",),bounds=(lower,upper),max_nfev=70,ftol=.00025,xtol=.00025)
    report={"fold":fold,"center_bound_arcsec":center_bound,"id":t["id"],"ra":t["ra"],"dec":t["dec"],"proposal_radius_arcsec":R,
            "fit_aperture_arcsec":radius,"training_pixels":int(train.sum()),"heldout_pixels":int(test.sum()),
            "vis_assumed_psf_fwhm_arcsec":.17,"nonzero_flag_fraction":float((flags[valid]!=0).mean()),
            "detected_knots":knots.tolist(),"heldout_blocks":int(len(set(zip(bx[test],by[test]))))}
    render=[]
    scores={}
    for kind,fit in [("lens",fitted),("knots",free)]:
        pred,weights,sky=solve(fit.x,kind,True)
        res=(values-pred)/errors
        render.append((pred*med,res));scores[kind]=res**2
        report[kind]={"parameters":fit.x.tolist(),"amplitudes_native":(weights*med).tolist(),
                      "sky_native":float(sky*med),"train_mean_sq":float(np.mean(res[train]**2)),
                      "test_mean_sq":float(np.mean(res[test]**2)),"optimizer_status":int(fit.status)}
    diff=scores["lens"]-scores["knots"]
    report["test_delta_lens_minus_knots"]=float(diff[test].mean())
    blockdiff=[float(np.mean(diff[test&(bx==a)&(by==b)])) for a,b in sorted(set(zip(bx[test],by[test])))]
    report["heldout_blocks_prefer_knots"]=sum(d>0 for d in blockdiff)
    report["block_mean_differences"]=blockdiff
    report["lens_multistart_train_means"]=[float(2*s.cost/train.sum()) for s in sols]
    if fold!=0 or center_bound!=.35:
        tag=f"fold{fold}"+("_wide" if center_bound!=.35 else "")
        (out/(t["id"]+"_"+tag+".json")).write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        return (t["id"]+"_"+tag,report["test_delta_lens_minus_knots"])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,5,figsize=(11.2,2.8),layout="constrained")
    limits=np.percentile(data[valid],[7,99])
    ims=[data,render[0][0],render[0][1],render[1][0],render[1][1]]
    titles=["Observed VIS","One-source lens","Lens residual / RMS","Independent knots","Knots residual / RMS"]
    for j,(ax,z,title) in enumerate(zip(axes,ims,titles)):
        res=j in [2,4]
        ax.imshow(np.where(valid,z,np.nan),origin="lower",cmap="RdBu_r" if res else "gray",
                  vmin=-5 if res else limits[0],vmax=5 if res else limits[1])
        ax.set_title(title,fontsize=9);ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle(t["id"],fontsize=9)
    fig.savefig(out/(t["id"]+".png"),dpi=125);plt.close(fig)
    (out/(t["id"]+".json")).write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    return (t["id"],report["test_delta_lens_minus_knots"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--output",type=Path,default=HERE)
    parser.add_argument("--ids",nargs="*")
    parser.add_argument("--folds",type=int,nargs="+",default=[0,1,2,3,4])
    parser.add_argument("--center-bound",type=float,default=.35)
    args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True);args.cache.mkdir(exist_ok=True,parents=True)
    ts=targets()
    jobs=[(ts[k],args.cache,args.output,fold,args.center_bound) for fold in args.folds for k in (args.ids or list(ts))]
    with ProcessPoolExecutor(args.jobs) as pool:
        for target,delta in pool.map(fit_one,jobs):
            print(target,round(delta,4),flush=True)
