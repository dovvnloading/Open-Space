"""Reproduce Rubin source vetting, uncertainty stress tests and astrometric plots."""
import csv
from datetime import date, timedelta
from itertools import product

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cutouts import ROOT, DATA, choose_controls, rows, read_cutout

def number(r, key):
    return float(r[key])

def output_csv(name, records):
    with (ROOT/"output"/name).open("w", newline="") as stream:
        w = csv.DictWriter(stream,fieldnames=list(records[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(records)

def template_blend(mid):
    """Measure local template signal at the target and weighted PSF flux per component."""
    (h, a), (_, v), (_, _) = read_cutout(mid,"Template")
    yy,xx=np.indices(a.shape)
    rho=np.hypot(xx-(float(h["CRPIX1"])-1), yy-(float(h["CRPIX2"])-1))
    core=(rho<=4)&np.isfinite(a)&np.isfinite(v)&(v>0)
    sky=np.nanmedian(a[(rho>=8)&(rho<12)])
    sig=np.sqrt(np.sum(v[core]))
    return float(np.sum(a[core]-sky)),float(sig)
    
def main():
    fitfile=ROOT/"output"/"fits.csv"
    with fitfile.open() as stream:
        allfits=list(csv.DictReader(stream))
    targets=rows("target.csv")
    tids={r["measurement_id"] for r in targets}
    bykey={(r["id"],r["kind"],r["variant"]): r for r in allfits}
    control=choose_controls(targets,rows("peer_pool.csv"))
    if [(r["target_id"],r["control_id"]) for r in rows("controls.csv")] != [
        (r["target_id"],r["control_id"]) for r in control]:
        raise ValueError("Control selection does not reproduce")
    result=[]
    for t in targets:
        mid=t["measurement_id"]
        fit=bykey[mid,"Science","nominal"]
        peers=[bykey[r["control_id"],"Science","nominal"] for r in control if r["target_id"]==mid]
        template_flux,template_error=template_blend(mid)
        result.append(dict(id=mid,utc_date=str(date(1858,11,17)+timedelta(days=int(float(t["mjd"])))),
                       mjd=t["mjd"],band={"3":"i","4":"z"}[t["band"]],visit=t["visit"],
                       detector=t["detector"],snr=t["snr"],toporange_AU=t["toporange"],
                       east=fit["east_arcsec"],north=fit["north_arcsec"],pa=fit["pa_deg"],
                       separation=fit["separation_arcsec"],flux_fraction=fit["faint_flux_fraction"],
                       delta_chi2=fit["delta_chi2"],peers_above=sum(number(p,"delta_chi2")>=number(fit,"delta_chi2") for p in peers),
                       n_peers=len(peers),template_aperture_nJy=template_flux,template_aperture_sigma_nJy=template_error))
    output_csv("astrometry.csv",result)
    # Fit fixed separation vector on all independent exposures, minimize over 180-degree
    # component ambiguity, with an added independent error per coordinate per exposure.
    sci=[bykey[t["measurement_id"],"Science","nominal"] for t in targets]
    xy=np.array([[number(r,"east_arcsec"),number(r,"north_arcsec")] for r in sci])
    covariance=np.array([[[number(r,"sigma_east")**2,number(r,"cov_east_north")],
                           [number(r,"cov_east_north"),number(r,"sigma_north")**2]] for r in sci])
    floor_results=[]
    for floor in (0.,.03,.06,.10):
        weights=np.linalg.inv(covariance + np.eye(2)*floor**2)
        combined=np.linalg.inv(weights.sum(axis=0))
        best=np.inf
        for signs in product((1,-1),repeat=len(sci)-1):
            vec=xy*np.array([1,*signs])[:,None]
            center=combined@np.einsum("ijk,ik->j",weights,vec)
            delta=vec-center
            chi=np.einsum("ij,ijk,ik->",delta,weights,delta)
            best=min(chi,best)
        floor_results.append(dict(independent_floor_arcsec=floor,constant_vector_chi2=best,dof=2*len(sci)-2))
    output_csv("constant_vector.csv",floor_results)
    sensitivity=[]
    for kind,variant in [("Science","nominal"),("Science","inner"),("Science","broadened"),("Difference","nominal")]:
        fits=[bykey[t["measurement_id"],kind,variant] for t in targets]
        sensitivity.append(dict(kind=kind,variant=variant,
            endpoint_pa_change_deg=number(fits[-1],"pa_deg")-number(fits[0],"pa_deg"),
            min_delta_chi2=min(number(f,"delta_chi2") for f in fits),
            max_delta_chi2=max(number(f,"delta_chi2") for f in fits)))
    output_csv("sensitivity.csv",sensitivity)
    fig,axes=plt.subplots(1,2,figsize=(9.0,3.8),constrained_layout=True)
    for label,color,cmarker in [("Science","#183e6c","o"),("Difference","#b14a16","s")]:
        fits=[bykey[t["measurement_id"],label,"nominal"] for t in targets]
        x=[number(t,"mjd")-61218 for t in targets]
        for ax,y,e in [(axes[0],"pa_deg","sigma_pa_deg"),(axes[1],"separation_arcsec","sigma_sep")]:
            ax.errorbar(x,[number(f,y) for f in fits],yerr=[number(f,e) for f in fits],
                        color=color,fmt=cmarker,fillstyle=("full" if label=="Science" else "none"),
                        capsize=2,label=label,alpha=.88)
    for ax in axes:
        ax.set_xlabel("Days after MJD 61218"); ax.grid(alpha=.2)
    axes[0].set_ylabel("Position angle east of north (deg)")
    axes[1].set_ylabel("Projected separation (arcsec)")
    axes[0].legend(frameon=False)
    fig.savefig(ROOT/"output"/"rubin_astrometry.png",dpi=160)
    plt.close(fig)
    return result,floor_results
if __name__=="__main__":
    main()
