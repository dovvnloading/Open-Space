"""Fit the radial profile, then measure the frozen outer-field and energy predictions."""
import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
from scipy.special import ndtr

from search import event_path, geometry, photons, rows, write_csv

DISTANCE = 157.3
BANDS = np.array([0.7, 0.9, 1.1, 1.3])


def fit_profile(array, out):
    """Use fixed radial coverage and a smoothly varying background in squared radius."""
    t, dx, dy, e = array
    r0 = np.sqrt(826.6*t/DISTANCE)
    u = (dx*dx + dy*dy) / r0**2
    use = (e>=70) & (e<=300) & (r0*np.sqrt(0.5)>50) & (r0*np.sqrt(1.65)<550)
    edges = np.arange(.5, 1.65001, .025)
    mid = (edges[1:]+edges[:-1])/2
    def fit(mask):
        y = np.histogram(u[mask], edges)[0]
        def model(p):
            base = np.exp(p[0] + p[1]*(mid-1) + p[2]*(mid-1)**2)
            line = np.diff(ndtr((edges - p[3])/p[4])) * p[5]
            return base, line
        def nll(p):
            a,b=model(p)
            return np.sum(a+b-y*np.log(a+b))
        bounds = [(0,9),(-5,5),(-6,6),(.79,1.19),(.022,.22),(0,2500)]
        fit = minimize(nll, [np.log(max(1,np.median(y))),-.6,0,1,.07,150],
                       method="Nelder-Mead", bounds=bounds, options={"maxiter":3500,"xatol":1e-7})
        # Start bounded quasi-Newton at the simplex optimum to check convergence.
        final = minimize(nll, fit.x, method="L-BFGS-B", bounds=bounds)
        p=final.x
        return y, model(p), p
    y,(b,l),p=fit(use)
    results = []
    for label, mask in [("all",use),("before_4h",use & (t<4*3600)),
                        ("after_4h",use & (t>=4*3600))]:
        _,_,par = fit(mask)
        results.append(dict(subset=label, distance_pc=round(DISTANCE/par[3],2),
                            sigma_u=round(par[4],4), fitted_excess=round(par[5],1)))
    write_csv(out/"profile-fits.csv", results)
    write_csv(out/"profile.csv", [dict(u=round(c,5), photons=int(n),background=round(bi,4),
                                      line=round(li,4)) for c,n,bi,li in zip(mid,y,b,l)])


def snapshot_counts(hd, row, file_id):
    events = hd["EVENTS"]
    x0,y0,sx,sy,t0 = geometry(row, events)
    d = events.data
    time = d["TIME"]-t0
    xx, yy = (d["X"]-x0)*sx, (d["Y"]-y0)*sy
    rr = np.hypot(xx,yy)
    bad = np.zeros((600,600), dtype=bool)
    for pix in hd["BADPIX"].data:
        x, y, length = (int(pix[k]) for k in ("RAWX","RAWY","YEXTENT"))
        bad[y:y+max(1,length),x] = True
    yy_grid,xx_grid = np.mgrid[2:598:2,2:598:2]
    valid = ~bad[yy_grid,xx_grid]
    pixels = np.stack([xx_grid[valid],yy_grid[valid],np.ones(valid.sum())],axis=1)
    results=[]
    for k,gti in enumerate(hd["GTI"].data):
        start,stop = gti["START"]-t0, gti["STOP"]-t0
        active = (time>=start)&(time<=stop)
        if min(abs(start),abs(stop))>432000 or stop<200 or active.sum()<40:
            continue
        observed = np.stack([d["RAWX"][active],d["RAWY"][active],np.ones(active.sum())],axis=1)
        positions = np.stack([xx[active],yy[active]],axis=1)
        # The residual includes dither. Fit each GTI separately using all events.
        coeff = np.linalg.lstsq(observed,positions,rcond=None)[0]
        residual = np.sqrt(np.mean(np.sum((observed@coeff-positions)**2,axis=1)))
        grid_r=np.hypot(*(pixels@coeff).T)
        for label,emin,emax,rmin,rmax in (("outer_soft",70,300,550,850),
                     ("outer_hard",301,700,550,850),("inner_0.5_0.7",50,69,50,550),
                     ("inner_soft",70,300,50,550),("inner_hard",301,700,50,550)):
            if np.sqrt(826.6*max(stop,0)/DISTANCE*BANDS[-1])<rmin:
                continue
            area=np.zeros(3)
            grid=grid_r[(grid_r>=rmin)&(grid_r<rmax)]**2
            lo,hi=max(start,200),min(stop,432000)
            if hi<=lo:continue
            for frac in (np.arange(5)+.5)/5:
                area+=np.histogram(grid*DISTANCE/(826.6*(lo+frac*(hi-lo))), BANDS)[0]
            if area[1]==0 or area[0]+area[2]==0:
                continue
            keep=(active & (time>=lo)&(time<=hi)&(d["GRADE"]<=12)&(d["PI"]>=emin)&
                  (d["PI"]<=emax)&(rr>=rmin)&(rr<rmax))
            n=np.histogram(rr[keep]**2*DISTANCE/(826.6*time[keep]), BANDS)[0]
            alpha=area[1]/(area[0]+area[2])
            results.append(dict(obsid=file_id, gti=k, band=label, hours=round((lo+hi)/7200,4),
                          rms_arcsec=round(residual,2),inner=int(n[0]),on=int(n[1]),outer=int(n[2]),
                          alpha=round(alpha,5),expected=round(alpha*(n[0]+n[2]),4)))
    return results


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",type=Path,default=Path("xrt-cache"))
    ap.add_argument("--out",type=Path,default=Path("."))
    args=ap.parse_args()
    args.cache.mkdir(exist_ok=True,parents=True)
    args.out.mkdir(exist_ok=True,parents=True)
    row=next(r for r in rows("sample.csv") if r["name"]=="180728A")
    arrays,counts=[],[]
    for rec in rows("event-files.csv"):
        if rec["name"]!=row["name"]:continue
        path=event_path(args.cache,rec)
        arrays.append(photons(path,row))
        with fits.open(path) as hd:
            counts.extend(snapshot_counts(hd,row,rec["obsid"]))
    fit_profile(np.concatenate(arrays,axis=1),args.out)
    write_csv(args.out/"annuli.csv",counts)
    for key in sorted({c["band"] for c in counts}):
        sub=[c for c in counts if c["band"]==key]
        on=sum(c["on"] for c in sub)
        bg=sum(c["expected"] for c in sub)
        variance=on+sum(c["alpha"]**2*(c["inner"]+c["outer"]) for c in sub)
        print(key,len(sub),on,round(bg,1),round((on-bg)/np.sqrt(variance),2))


if __name__=="__main__":
    main()
