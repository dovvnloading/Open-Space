"""Additive, offline ZTF exposure-paired holdout and TESS waveform replay.

Reuses raw PR #24 records READ-ONLY; creates no files without --output.
Optional --asassn-dir is the directory of skypatrol CSVs made by harvest.py.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0,str(ROOT))
from replay import bjd  # Original archive's exact HJD(UTC) -> BJD(TDB) conversion.

FREQ = 1/.15956200154633035
EPOCH = 2459499.9473282346
TARGET = 333114400007906
CONTROL = 333114400172086


def matrix(t, frequencies):
    angle = 2*np.pi*(np.asarray(t)-EPOCH)
    return np.column_stack([np.ones(len(angle))] + [a for f in frequencies
                         for a in [np.cos(angle*f),np.sin(angle*f)]])


def weighted(X,y,w):
    sw=np.sqrt(w)
    return np.linalg.lstsq(X*sw[:,None],np.asarray(y)*sw,rcond=None)[0]


def photometry():
    a=pd.read_csv(ROOT/"data/ztf-field.csv")
    a=a[(a.catflags==0)&(a.magerr<.18)&(a.mag>12)]
    left=a[a.oid==TARGET]
    right=a[a.oid==CONTROL]
    m=left.merge(right,on="expid",suffixes=("_t","_c")).sort_values("hjd_t").copy()
    m=m[m.hjd_t>2451544.5+(2017.6-2000)*365.25+400]
    m["delta"]=m.mag_t-m.mag_c
    m["delta_error"]=np.hypot(m.magerr_t,m.magerr_c)
    t=m.hjd_t.to_numpy()
    y=m.delta.to_numpy()
    low=np.searchsorted(t,t-25)
    high=np.searchsorted(t,t+25)
    m["res"]=y-np.array([np.median(y[a:b]) for a,b in zip(low,high)])
    cutoff=5*max(.04,np.quantile(np.abs(m.res),.75))+.1
    m=m[np.abs(m.res)<cutoff].copy()
    t=bjd(m.hjd_t.to_numpy(),"sun")
    X=matrix(t,[FREQ])
    y=m.res.to_numpy()
    w=1/(m.delta_error.to_numpy()**2+.06**2)
    coef=weighted(X,y,w)
    r=np.corrcoef(y,X@coef)[0,1]
    phase=float(np.arctan2(coef[2],coef[1])/(2*np.pi))
    # The phase origin in this analysis is *brightness* maximum, where
    # negative ZTF magnitudes have expected phase 0.5 cycles.
    early=m.hjd_t.to_numpy()<2460300
    beta=weighted(X[early],y[early],w[early])
    actual=y[~early]; model=X[~early]@beta
    cor=float(np.corrcoef(actual,model)[0,1])
    weights=w[~early]
    actual=actual-np.average(actual,weights=weights)
    model=model-np.average(model,weights=weights)
    inner=np.sum(weights*actual*model)
    month=np.floor(m.hjd_t.to_numpy()[~early]/30)
    scores=np.array([np.sum(weights[month==q]*actual[month==q]*model[month==q])
                     for q in np.unique(month)])
    rng=np.random.default_rng(1234)
    random_sign=rng.choice([-1,1],size=(100000,len(scores)))
    exceed=int(np.sum((random_sign@scores)>=inner))
    # Astrometric FIRST-MOMENT blend diagnostic, with comparison star positions
    # removed *on the same exposure*, and independent 150-day offsets.
    skydec=np.deg2rad(-14.315539)
    east=(m.ra_t.to_numpy()-m.ra_c.to_numpy())*np.cos(skydec)*3600
    north=(m.dec_t.to_numpy()-m.dec_c.to_numpy())*3600
    seasons=np.floor((m.hjd_t.to_numpy()-2459500)/150).astype(int)
    M=np.column_stack([seasons==q for q in np.unique(seasons)]+
                      [(X@coef-coef[0])/np.hypot(*coef[1:])])
    Y=np.column_stack([east,north,y])
    def position(ix):
        a=np.linalg.lstsq(M[ix],Y[ix],rcond=None)[0][-1]
        return -a[:2]/(.4*np.log(10)*a[2])
    pos=position(np.arange(len(m)))
    groups=[np.flatnonzero(np.floor(m.hjd_t.to_numpy()/30)==i)
            for i in np.unique(np.floor(m.hjd_t.to_numpy()/30))]
    rng=np.random.default_rng(1312)
    boots=[]
    for _ in range(1500):
        idx=np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))])
        boots.append(position(idx))
    out=dict(n=int(len(m)),control_oid=str(CONTROL),matched_exposures=int(len(m)),
             period_days=1/FREQ,semiamplitude_mag=float(np.hypot(*coef[1:])),
             mag_phase_cycles=phase,correlation=float(r),
             train_before_hjd=2460300,train_n=int(early.sum()),test_n=int((~early).sum()),
             train_amp_mag=float(np.hypot(*beta[1:])),
             heldout_correlation=cor,heldout_month_blocks=int(len(scores)),
             heldout_one_sided_signflip_exceedances=exceed,heldout_signflip_trials=100000,
             heldout_p=(exceed+1)/100001,
             inferred_variable_offset_east_north_arcsec=pos.tolist(),
             monthly_bootstrap_offset_95_east_north_arcsec=np.quantile(boots,[.025,.975],axis=0).tolist(),
             centroid_model="PSF-fit coordinates approximated as first moments; systematics not included")
    years={}
    for year in sorted(set(np.floor((m.hjd_t-2451544.5)/365.25+2000).astype(int))):
        q=np.floor((m.hjd_t-2451544.5)/365.25+2000).astype(int)==year
        if q.sum()<10:continue
        beta_y=weighted(X[q],y[q],w[q])
        years[str(year)]=dict(n=int(q.sum()),semiamplitude_mag=float(np.hypot(*beta_y[1:])),
              mag_phase_cycles=float(np.arctan2(beta_y[2],beta_y[1])/(2*np.pi)))
    out["seasons"]=years
    return out


def tess_archived():
    o={}
    for s in (80,92):
        d=pd.read_csv(ROOT/f"data/tess-{s}.csv.gz")
        t=d.btjd.to_numpy()+2457000
        y=d.corrected_e_per_s.to_numpy()
        X=matrix(t,[FREQ/2,FREQ,2*FREQ,3*FREQ])
        a=np.linalg.lstsq(X,y,rcond=None)[0]
        row=dict(n=len(t),cadence_minutes=float(np.median(np.diff(t))*1440))
        for i,key in enumerate(("f/2","f","2f","3f")):
            b=a[2*i+1:2*i+3]
            row[key]=dict(semiamplitude_e_per_s=float(np.linalg.norm(b)),
                phase_cycles=float((np.arctan2(b[1],b[0])/(2*np.pi)+.5)%1-.5))
        o[str(s)]=row
    return o


def asassn(directory):
    import hashlib
    ids=("661428620105","661431894529","661431894524")
    paths=[directory/f"asas-{k}.csv" for k in ids]
    if any(not p.exists() for p in paths):
        return {"status":"unavailable","missing":[p.name for p in paths if not p.exists()]}
    series=[pd.read_csv(p) for p in paths]
    shared=[]
    for i in range(3):
        for j in range(i+1,3):
            x=series[i].merge(series[j],on="image_id",suffixes=("_a","_b"))
            g=(x.quality_a=="G")&(x.quality_b=="G")&(x.phot_filter_a=="g")&(x.phot_filter_b=="g")
            shared.append(dict(pair=[ids[i],ids[j]],n=int(g.sum()),
                magnitude_correlation=float(np.corrcoef(x.mag_a[g],x.mag_b[g])[0,1])))
    v=series[0]
    v=v[(v.quality=="G")&(v.mag_err>0)&(v.mag_err<.2)&(v.phot_filter=="g")&
        (v.jd>2458400)&(v.jd<2460809)].sort_values("jd").copy()
    t=v.jd.to_numpy();y=v.mag.to_numpy()
    lo=np.searchsorted(t,t-20);hi=np.searchsorted(t,t+20)
    v["res"]=y-np.array([np.median(y[i:j]) for i,j in zip(lo,hi)])
    v=v[np.abs(v.res)<.7]
    t=v.jd.to_numpy()+.0008   # Approximate JD(UTC) -> TDB; no site ephemeris in ASAS download.
    y=v.res.to_numpy()
    w=1/(v.mag_err.to_numpy()**2+.04**2)
    X=matrix(t,[FREQ])
    b=weighted(X,y,w)
    # For an external source, absolute phase comparisons require observatory
    # barycentric correction that is not provided. No phase confirmation is made.
    return dict(status="tested_without_site_barycentric_correction",photometry_sha256=[
        __import__("hashlib").sha256(p.read_bytes()).hexdigest() for p in paths],
        n=int(len(v)),semiamplitude_mag=float(np.hypot(*b[1:])),
        phase_approx_jd_cycles=float(np.arctan2(b[2],b[1])/(2*np.pi)),
        residual_rms_mag=float(np.std(y-X@b)),
        shared_aperture_pairs=shared)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path)
    p.add_argument("--asassn-dir",type=Path)
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    gaia=pd.read_csv(ROOT/"data/gaia-field.csv",dtype={"source_id":str})
    target=gaia[gaia.source_id=="4104113350446549888"].iloc[0]
    out={"ztf_exposure_paired":photometry(),"tess_original_pca_harmonic":tess_archived(),
         "gaia_dr3_target_has_epoch_photometry":str(target.has_epoch_photometry).lower()=="true"}
    if args.asassn_dir:out["asas_sn"]=asassn(args.asassn_dir)
    if args.check:
        z=out["ztf_exposure_paired"];tt=out["tess_original_pca_harmonic"]
        assert z["n"]==246 and z["test_n"]==94 and z["train_n"]==152
        assert .042<z["semiamplitude_mag"]<.052
        assert z["heldout_correlation"]>.50 and z["heldout_p"]<.003
        assert all(.4<tt[str(s)]["2f"]["semiamplitude_e_per_s"]<.8 for s in (80,92))
        assert not out["gaia_dr3_target_has_epoch_photometry"]
        print("Passed exposure-paired and TESS archive regression checks.",file=sys.stderr)
    report=json.dumps(out,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(report)
    else:print(report)


if __name__=="__main__":
    main()
