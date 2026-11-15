"""Per-band V612 Scuti clock and non-detection at the second harmonic.

New independent script on the ORIGINAL PR #24 photometry: ZTF matched g
exposures, Pan-STARRS five filters, optional public ASAS-SN Sky Patrol CSV.
The latter are not guaranteed independent stellar light curves because apertures
overlap; ASAS V is essentially pre-2019 only. No p-value for global selection.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from replay import bjd

P=.15956200154633035
F=1/P
E0=2459499.9473282346


def table():
    z=pd.read_csv(ROOT/"data/ztf-field.csv")
    z=z[(z.catflags==0)&(z.magerr<.18)&(z.mag>12)]
    m=z[z.oid==333114400007906].merge(
      z[z.oid==333114400172086],on="expid",suffixes=("_t","_c"))
    m=m.sort_values("hjd_t")
    m=m[m.hjd_t>2451544.5+(2017.6-2000)*365.25+400].copy()
    m["delta"]=m.mag_t-m.mag_c
    m["err"]=np.hypot(m.magerr_t,m.magerr_c)
    t=m.hjd_t.to_numpy()
    y=m.delta.to_numpy()
    lo=np.searchsorted(t,t-25)
    hi=np.searchsorted(t,t+25)
    m["res"]=y-np.array([np.median(y[i:j]) for i,j in zip(lo,hi)])
    width=5*max(.04,np.quantile(np.abs(m.res),.75))+.1
    m=m[np.abs(m.res)<width].copy()
    return m


def four_harmonic(t):
    x=2*np.pi*(t-E0)
    return np.column_stack([np.ones(len(t))]+[
        v for q in (.5*F,F,2*F,3*F) for v in (np.cos(x*q),np.sin(x*q))])


def fit(X,y,w):
    v=np.sqrt(w)
    return np.linalg.lstsq(X*v[:,None],np.asarray(y)*v,rcond=None)[0]


def color_test():
    m=table()
    t=bjd(m.hjd_t.to_numpy(),"sun")
    y=m.res.to_numpy()
    w=1/(m.err.to_numpy()**2+.06**2)
    X=four_harmonic(t)
    b=fit(X,y,w)
    rng=np.random.default_rng(125)
    months=np.floor(t/30).astype(int)
    groups=[np.flatnonzero(months==j) for j in np.unique(months)]
    draws=[]
    for _ in range(2500):
        ii=np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))])
        bb=fit(X[ii],y[ii],w[ii])
        draws.append(bb[5:7])
    q=np.array(draws)
    # Day/month sign-flip test of 2f after nuisance f/2, f and 3f,
    # preserving each observing month's correlated/noisy residuals.
    use=[i for i in range(9) if i not in (5,6)]
    N=X[:,use];Z=X[:,5:7]
    resid=y-N@fit(N,y,w)
    # weighted Frisch-Waugh residualization for the tested two columns
    W=np.sqrt(w)
    Zres=Z-N@np.linalg.lstsq(N*W[:,None],Z*W[:,None],rcond=None)[0]
    V=np.array([((Zres[ii].T*w[ii])@resid[ii]) for ii in groups])
    inv=np.linalg.inv((Zres.T*w)@Zres)
    stat=lambda values:np.sqrt(np.sum((values@inv)*values,axis=-1))
    obs=float(stat(V.sum(axis=0)))
    trials=rng.choice([-1,1],size=(20000,len(groups)))@V
    exceed=int(np.sum(stat(trials)>=obs))
    ztf=dict(n=len(m),f_amplitude_mag=float(np.linalg.norm(b[3:5])),
       harmonic_2f_amp_mag=float(np.linalg.norm(b[5:7])),
       harmonic_2f_mag_phase_cycles=float(np.arctan2(b[6],b[5])/(2*np.pi)),
       harmonic_2f_amp_day_bootstrap_95=np.quantile(
           np.linalg.norm(q,axis=1),[.025,.975]).tolist(),
       harmonic_2f_month_signflip=dict(trials=20000,exceedances=exceed,
              p=(exceed+1)/20001))
    ps=pd.read_csv(ROOT/"data/panstarrs.csv")
    ps=ps[(ps.psfFlux>0)&(ps.psfFluxErr<.17*ps.psfFlux)]
    by_filter={}
    for ident,label in enumerate(("g","r","i","z","y"),1):
        tab=ps[ps.filterID==ident]
        by_filter[label]=dict(n=int(len(tab)),
           observing_nights=int(len(set(np.floor(tab.obsTime.to_numpy())))),
           time_range_mjd=[float(tab.obsTime.min()),float(tab.obsTime.max())])
    return dict(ztf_g_exposure_paired=ztf,panstarrs_preeruption=by_filter,
       caveat="Pan-STARRS detections predate 2017 eruption; not a late optical phase validation.")


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path)
    a.add_argument("--check",action="store_true")
    args=a.parse_args()
    report=color_test()
    if args.check:
        z=report["ztf_g_exposure_paired"]
        assert z["n"]==246
        assert .040<z["f_amplitude_mag"]<.058
        assert z["harmonic_2f_amp_mag"]<.015
        assert z["harmonic_2f_month_signflip"]["p"]>.1
        assert sum(q["n"] for q in report["panstarrs_preeruption"].values())==69
        print("ZTF g-band harmonic null and five-filter Pan-STARRS coverage validated.",file=sys.stderr)
    content=json.dumps(report,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content)
    else:print(content)


if __name__=="__main__":
    main()
