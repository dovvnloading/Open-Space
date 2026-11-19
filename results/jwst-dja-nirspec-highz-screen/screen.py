#!/usr/bin/env python3
"""Exhaustive exploratory DJA-v4.4 1D redshift template screen.

Run AFTER retrieve.py download --all and retrieve.py freeze:
  python screen.py --cache /data/dja --out /data/dja/fit
  python screen.py --self-test

Fits ALL nine bands: 80,367 source spectra when every frozen input is present.
Uses 1D fluxes, full_err and valid masks, not DJA redshift as a fitting prior.
Outputs are deliberately named "proxy": simplified stellar templates, diagonal
errors, uncalibrated model prior and missing 2D/photometry cannot support any
discovery or scientifically secure probability of the interloper solution.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import log_ndtr, logsumexp

HERE=Path(__file__).resolve().parent
SOURCES=json.loads((HERE/"frozen_sources.json").read_text())
HI=np.round(np.arange(10.0,18.001,0.10),2)
LOW=np.round(np.arange(2.6,5.501,0.15),2)
EM=np.round(np.arange(0.5,9.501,0.35),2)
PRIORS={"high":.02,"dusty":.55,"emission":.28,"dwarf":.10,"artifact":.05}
OUTCOL=("band index root file srcid objid ra dec grade dja_z_best valid_bins "
        "two_sided q16_high q50_high q84_high p_gt14_proxy p_gt12_proxy "
        "p_lowz_gal_proxy p_dwarf_proxy p_artifact_proxy "
        "p_gt12_prior_0p002 p_gt12_prior_0p20 status").split()


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as inp:
        for chunk in iter(lambda:inp.read(2**20),b""):
            h.update(chunk)
    return h.hexdigest()


def gaussian(w,center,res):
    width=max(center/(2.355*res),.0003)
    return np.exp(-.5*((w-center)/width)**2)


def bank(w,band):
    w=np.asarray(w,float)
    R=100 if band=="prism" else 1000
    lib={}
    lib["high"]=np.asarray([
        (w/2)**s * np.where(w<.121567*(1+z),.008,
          np.exp(-tau*.025/np.maximum(w-.121567*(1+z),.025)))
        for z in HI for s in (-.5,.8) for tau in (0.,.9)
    ])
    models=[]
    for z in LOW:
        balmer=.4*(1+z)
        for slope in (-.8,0.,1.3):
            for blue_ratio in (.10,.35,.70):
                base=(w/2)**slope*np.where(w<balmer,blue_ratio,1.)
                redden=np.exp(-.6*(np.maximum(.35,2/w)-1.))
                lines=(.4*gaussian(w,.3727*(1+z),R)
                     + .22*gaussian(w,.4861*(1+z),R)
                     + .25*gaussian(w,.4959*(1+z),R)
                     + .75*gaussian(w,.5007*(1+z),R)
                     + .4*gaussian(w,.6563*(1+z),R))
                for a in (0.,2.):
                    models.append(base*redden+a*lines)
    lib["dusty"]=np.asarray(models)
    models=[]
    for z in EM:
        lines=(.6*gaussian(w,.3727*(1+z),R)
             + .2*gaussian(w,.4861*(1+z),R)
             + .3*gaussian(w,.4959*(1+z),R)
             + 1.0*gaussian(w,.5007*(1+z),R)
             + .5*gaussian(w,.6563*(1+z),R))
        for a in (1.,6.):
            models.append(np.ones_like(w)+a*lines)
    lib["emission"]=np.asarray(models)
    dips=np.ones_like(w)
    for x,depth,width in ((1.4,.65,.085),(1.9,.72,.09),
                           (2.7,.5,.18),(3.3,.65,.12)):
        dips*=1-depth*np.exp(-.5*((w-x)/width)**2)
    lib["dwarf"]=np.asarray([(w/2)**s*dips for s in (-1.5,-.5,.5,1.5)])
    lib["artifact"]=np.asarray([np.ones_like(w), w/2, 2/w,
                                np.where(w<2.45,.2,1.)])
    return lib


def amplitude_integral(templates,f,e,valid):
    # Whitened normal model with analytically marginalized A>=0:
    # log int exp[-(f-A*t)^T W(f-A*t)/2] dA
    # Data-only factor omitted; broad common amplitude scale is a proxy.
    w=np.where(valid,1/np.maximum(e,1e-12)**2,0)
    tm=templates@valid.astype(float)/np.maximum(valid.sum(axis=0),1)
    tm=np.maximum(tm,1e-12)
    precision=np.maximum(((templates**2)@w)/(tm**2),1e-20)
    projection=(templates@(f*w))/tm
    sn=projection/np.sqrt(precision)
    return .5*sn**2-.5*np.log(precision)+log_ndtr(sn)


def fit(w,f,e,valid,models):
    n=f.shape[1]
    ll={name:amplitude_integral(arr,f,e,valid) for name,arr in models.items()}
    raw={k:logsumexp(v,axis=0)-np.log(v.shape[0]) for k,v in ll.items()}
    def probabilities(prior_high):
        norm=(1-prior_high)/(1-PRIORS["high"])
        logs=np.stack([raw[k]+np.log(prior_high if k=="high" else PRIORS[k]*norm)
                       for k in PRIORS])
        den=logsumexp(logs,axis=0)
        return {k:np.exp(logs[i]-den) for i,k in enumerate(PRIORS)},den
    classes,den=probabilities(.02)
    priorlo,_=probabilities(.002)
    priorhi,_=probabilities(.20)
    nv=ll["high"].shape[0]//len(HI)
    per_z=logsumexp(ll["high"].reshape(len(HI),nv,n),axis=1)
    per_z+=np.log(PRIORS["high"])-np.log(ll["high"].shape[0])
    p_z_global=np.exp(per_z-den)
    p_z_given_high=p_z_global/np.maximum(p_z_global.sum(axis=0),1e-300)
    qs=np.array([np.interp([.16,.5,.84],np.cumsum(p_z_given_high[:,i]),HI)
                 for i in range(n)])
    edge=.121567*(1+HI[np.argmax(p_z_global,axis=0)])
    wave=np.asarray(w)[:,None]
    blue=((wave>edge-.18)&(wave<edge-.04)&valid).sum(axis=0)
    red=((wave>edge+.04)&(wave<edge+.18)&valid).sum(axis=0)
    high_gt12=p_z_given_high[HI>12].sum(axis=0)
    return dict(
        two_sided=(blue>=3)&(red>=3),
        q16_high=qs[:,0],q50_high=qs[:,1],q84_high=qs[:,2],
        p_gt14_proxy=p_z_global[HI>14].sum(axis=0),
        p_gt12_proxy=p_z_global[HI>12].sum(axis=0),
        p_lowz_gal_proxy=classes["dusty"]+classes["emission"],
        p_dwarf_proxy=classes["dwarf"],
        p_artifact_proxy=classes["artifact"],
        p_gt12_prior_0p002=priorlo["high"]*high_gt12,
        p_gt12_prior_0p20=priorhi["high"]*high_gt12,
        total_mass=sum(classes.values()),
    )


def self_test():
    w=np.linspace(.70,5.20,180)
    B=bank(w,"prism")
    assert B["high"].shape[0]==len(HI)*4
    assert all(np.isfinite(t).all() and (t>=0).all() for t in B.values())
    # Only an internal numerical invariant test; no injection/recovery
    # calibration or astronomical spectrum was used here.
    rng=np.random.default_rng(17)
    f=np.column_stack([B["high"][0]*.3,B["dusty"][0]*.5,
                       rng.normal(0,.05,len(w))])
    e=np.full(f.shape,.06); valid=np.ones(f.shape,bool)
    r=fit(w,f,e,valid,B)
    assert all(np.isfinite(r[k]).all() for k in r)
    assert np.allclose(r["total_mass"],1.,atol=1e-8)
    assert np.all(r["p_gt14_proxy"]<=r["p_gt12_proxy"]+1e-12)
    assert np.all(r["p_gt12_prior_0p002"]<=r["p_gt12_prior_0p20"]+1e-12)
    print("PASS internal numerical invariants (synthetic inputs only)")


def read_frozen(cache):
    receipt=json.loads((cache/"freeze_receipt.json").read_text())
    manifest=cache/"frozen_row_manifest.csv"
    if sha256(manifest)!=receipt["manifest_sha256"]:
        raise ValueError("Frozen source list changed")
    grouped={}
    with manifest.open(newline="") as inp:
        for row in csv.DictReader(inp):
            g=row["grating"].lower()
            b="prism" if g=="prism" else g+"-"+row["filter"].lower()
            grouped.setdefault(b,[]).append(row)
    if sum(map(len,grouped.values()))!=80367:
        raise ValueError("Expected complete 80,367-row freeze")
    return grouped,receipt


def run(cache,out,batch):
    from astropy.io import fits
    groups,freeze=read_frozen(cache)
    out.mkdir(parents=True,exist_ok=True)
    source_receipts=[]
    complete=0
    short=0
    mainfile=out/"screen_rows.csv"
    shortlist=out/"unvetted_ranked_1d.csv"
    with mainfile.open("w",newline="") as sf,shortlist.open("w",newline="") as cf:
        writer=csv.DictWriter(sf,fieldnames=OUTCOL)
        candidate=csv.DictWriter(cf,fieldnames=OUTCOL)
        writer.writeheader();candidate.writeheader()
        for spec in SOURCES["dja"]["spectral_shards"]:
            b=spec["band"]; path=cache/spec["filename"]
            if not path.exists():raise FileNotFoundError(path)
            with fits.open(path,memmap=True) as hdus:
                h=hdus[1].data
                fields=set(hdus[1].columns.names)
                required={"wave","flux","full_err","valid"}
                if not required.issubset(fields):
                    raise ValueError(f"{b}: missing fields {required-fields}")
                wave=np.asarray(h["wave"],float)
                expected=spec["expected_spectra"]
                if h["flux"].shape!=(len(wave),expected) or len(groups[b])!=expected:
                    raise ValueError(f"{b}: FITS/summary spectrum order mismatch")
                models=bank(wave,b)
                for i in range(0,expected,batch):
                    j=min(i+batch,expected)
                    f=np.asarray(h["flux"][:,i:j],float)
                    e=np.asarray(h["full_err"][:,i:j],float)
                    v=np.asarray(h["valid"][:,i:j],bool)
                    good=v&np.isfinite(f)&np.isfinite(e)&(e>0)
                    if "npix" in fields:
                        good &= np.asarray(h["npix"][:,i:j],float)>0
                    fitres=fit(wave,np.nan_to_num(f),np.nan_to_num(e),good,models)
                    for k in range(j-i):
                        meta=groups[b][i+k]
                        row=dict(band=b,index=int(meta["index"]),root=meta["root"],
                                 file=meta["file"],srcid=meta["srcid"],objid=meta["objid"],
                                 ra=meta["ra"],dec=meta["dec"],grade=meta["grade"],
                                 dja_z_best=meta["z_best"],valid_bins=int(good[:,k].sum()))
                        for key in OUTCOL:
                            if key in fitres:
                                val=fitres[key][k]
                                row[key]=bool(val) if key=="two_sided" else float(val)
                        row["status"]="UNVETTED_1D" if row["valid_bins"]>=20 else "INSUFFICIENT_1D"
                        writer.writerow(row)
                        if (row["status"]=="UNVETTED_1D" and row["two_sided"] and
                            row["p_gt12_proxy"]>=.5):
                            candidate.writerow(row);short+=1
                        complete+=1
            source_receipts.append(dict(band=b,sha256=sha256(path),rows=expected))
            print(f"{b}: fit {expected}; unvetted ranked={short}",flush=True)
    if complete!=80367:raise ValueError(f"Incomplete scan: {complete}/80367")
    (out/"analysis_receipt.json").write_text(json.dumps({
        "stage":"1D proxy only; not an astrophysical detection or complete null",
        "frozen_source_sha256":freeze["manifest_sha256"],
        "catalog_sha256":freeze["catalog_sha256"],
        "spectra_fit":complete,"unvetted_ranked":short,
        "results_sha256":sha256(mainfile),"priors":PRIORS,
        "highz_grid":HI.tolist(),"lowz_grid":LOW.tolist(),
        "line_grid":EM.tolist(),"inputs":source_receipts,
        "pending":["2D native slit + disjoint exposure checks",
                   "contamination/bad-pixel checks","NIRCam filter-convolved SEDs",
                   "calibrated brown-dwarf models","paper/ADS/arXiv crossmatch",
                   "held-out injection and interloper posterior calibration"]
    },indent=2)+"\n")
    print(json.dumps({"fit":complete,"unvetted_ranked":short}))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--batch-size",type=int,default=64)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:self_test()
    else:
        if not args.cache or not args.out:ap.error("--cache and --out required")
        if not 1<=args.batch_size<=512:ap.error("batch size must be 1..512")
        run(args.cache,args.out,args.batch_size)
