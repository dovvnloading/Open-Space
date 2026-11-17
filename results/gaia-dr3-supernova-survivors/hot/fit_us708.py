"""Independent centroid fit to US 708 in its archived SDSS optical spectrum."""
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares

from screen import CACHE, HERE
sys.path.insert(0,str(HERE.parent/"spectra"))
from fit_control import factor, vacuum

URL = "https://dr17.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/0870/spec-0870-52325-0575.fits"
AIR = {"HeI_4471":4471.50,"HeII_4542":4541.59,"HeII_4686":4685.68,"HeII_5412":5411.52,
       "HeII_4100":4100.04,"HeII_4339":4338.67,"HeII_4859":4859.32,"HeII_6560":6560.10}


def fit(wave,flux,weight,lines):
    segments=[]
    init=[]
    lo=[]
    hi=[]
    for name in lines:
        center=float(vacuum(np.array([AIR[name]]))[0])
        keep=(abs(wave-center)<32)&(weight>0)&np.isfinite(flux)
        x,y,w=wave[keep],flux[keep],weight[keep]
        continuum=np.median(y)
        segments.append((center,x,y,w))
        init.extend([continuum,0,max(1.,0.3*continuum),2.5])
        lo.extend([-np.inf,-np.inf,0.,0.7])
        hi.extend([np.inf,np.inf,np.inf,10.])
    def residual(p):
        delta=factor(p[0])
        out=[]
        for i,(center,x,y,w) in enumerate(segments):
            cont,slope,depth,width=p[4*i+1:4*i+5]
            model=cont+slope*(x-center)-depth*np.exp(-0.5*((x-center*delta)/width)**2)
            out.append((y-model)*w)
        return np.concatenate(out)
    best=min([least_squares(residual,[v,*init],bounds=([-1700,*lo],[1700,*hi]),x_scale='jac')
              for v in np.arange(-1500,1501,300)],key=lambda result:result.cost)
    dof=len(best.fun)-len(best.x)
    cov=np.linalg.pinv(best.jac.T@best.jac) * max(1,2*best.cost/dof)
    return dict(lines=list(lines),rv_kms=float(best.x[0]),formal_err_kms=float(np.sqrt(cov[0,0])),
                chi2=float(2*best.cost),dof=dof,amplitudes=[float(x) for x in best.x[3::4]],
                widths_A=[float(x) for x in best.x[4::4]])


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    path=CACHE/"us708_sdss_dr17.fits"
    if not path.exists():
        with urllib.request.urlopen(URL,timeout=90) as stream:
            path.write_bytes(stream.read())
    with fits.open(path) as hdus:
        coadd=hdus["COADD"].data
        w=10.**coadd["loglam"].astype(float)
        y=coadd["flux"].astype(float)
        weight=np.sqrt(np.maximum(0,coadd["ivar"].astype(float)))
        weight[coadd["and_mask"]!=0]=0
        results=[fit(w,y,weight,lines) for lines in (list(AIR), list(AIR)[:4],list(AIR)[4:],*[ [a] for a in AIR ])]
        meta=dict(url=URL,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  plate=870,mjd=52325,fiber=575,run2d="26",air_rest_A=AIR,sdss_heliorv_header_kms=hdus[0].header.get("HELIO_RV"),
                  results=results)
    (HERE/"us708_fit.json").write_text(json.dumps(meta,indent=2)+"\n")
    print([(a["lines"],round(a["rv_kms"],2),round(a["formal_err_kms"],2)) for a in results])

if __name__=="__main__": main()
