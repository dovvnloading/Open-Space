"""Remeasure the short-cadence trigger and its other visits at the discovery period."""
import argparse
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from pixel_evidence import measure_pixels

def remeasure(rec,root,ep,full_pixels):
    p,e,width=ep["period_days"],ep["epoch_btjd"],ep["duration_min"]/1440
    with fits.open(root/rec["filename"]) as h:
        data=h[1].data
        t=np.asarray(data["TIME"],float)
        yy,xx=np.indices(h[2].data.shape)
        x,y=rec["target_xy"]
        ap=(xx-x)**2+(yy-y)**2<=2.25
        flux=np.asarray(data["FLUX"][:,ap],float).sum(axis=1)
        ok=(data["QUALITY"]==0)&np.isfinite(t+flux)
        pix=measure_pixels(t,data["FLUX"],data["QUALITY"],[x,y],e,width*1440,p) if full_pixels else None
        t,f=t[ok],flux[ok];f=f/np.median(f)-1
        phase=((t-e+.5*p)%p)/p-.5
        # Remove the mean of each orbit from data *and* each regressor to fit orbit intercepts.
        _,g=np.unique(np.rint((t-e)/p).astype(int),return_inverse=True)
        counts=np.bincount(g)
        center=lambda a:a-(np.bincount(g,weights=a)/counts)[g]
        cols=[np.sin(2*np.pi*phase),np.cos(2*np.pi*phase),-(abs(phase)<=width/(2*p)).astype(float)]
        X=np.stack([center(a) for a in cols],axis=1)
        yy0=center(f)
        def fit(x):
            beta=np.linalg.lstsq(x,yy0,rcond=None)[0]
            resid=yy0-x@beta
            # Cluster by orbital cycle to allow correlated noise within each orbit.
            meat=np.array([np.bincount(g,weights=v*resid) for v in x.T])
            bread=np.linalg.pinv(x.T@x)
            cov=bread@(meat@meat.T)@bread
            return dict(coefficients_ppm=(1e6*beta).tolist(),cluster_se_ppm=(np.sqrt(np.diag(cov))*1e6).tolist())
        bins=np.minimum(39,((phase+.5)*40).astype(int))
        return dict(filename=rec["filename"],tic=rec["tic"],sector=rec["sector"],cadence_s=rec["cadence_s"],period_days=p,epoch_btjd=e,
                    duration_minutes=width*1440,points=len(t),cycles=len(counts),box_only=fit(X[:,2:]),harmonic_plus_box=fit(X),
                    flux_by_phase_ppm=(np.bincount(bins,weights=yy0,minlength=40)/np.bincount(bins,minlength=40)*1e6).tolist(),pixels=pix)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("results",type=Path);p.add_argument("cache",type=Path);p.add_argument("out",type=Path);a=p.parse_args()
    rows=[json.loads(f.read_text()) for f in sorted(a.results.glob("*.json"))]
    out=[]
    for r in rows:
        for trig in r.get("search",{}).get("periodic",[]):
            out.append(dict(trigger=r["filename"],crosschecks=[remeasure(q,a.cache,trig,r==q) for q in rows if q["tic"]==r["tic"] and q["status"]=="searched"]))
    a.out.write_text(json.dumps(out,indent=1,allow_nan=False)+"\n")
