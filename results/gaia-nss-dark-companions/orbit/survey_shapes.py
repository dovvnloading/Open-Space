"""External survey RV shape checks on all strict, repeat-observed screened stars."""
import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from fit_controls import HERE
sys.path.insert(0,str(HERE.parent/"lamost"))
from check_orbits import predictions


def read(p):
    with p.open() as f:return list(csv.DictReader(f))


def main(catalog, extension):
    directories=[HERE.parent]+([HERE.parent/"extension"] if extension else [])
    mass={}
    visit_groups=defaultdict(list)
    seen=set()
    hashes={}
    for directory in directories:
        mpath=directory/"mass_screen.csv" if directory.name=="extension" else directory/"gaia/mass_screen.csv"
        for r in read(mpath):
            mass[(r["source_id"],r["nss_solution_type"])]=float(r["m2_screen_msun"])
        for survey in ["sdss","desi"]:
            p=directory/"epochs"/f"{survey}_visits.csv"
            if not p.exists():
                continue
            hashes[str(p.relative_to(HERE.parent))]=hashlib.sha256(p.read_bytes()).hexdigest()
            for r in read(p):
                if not int(r["quality"]):continue
                sid,ins=r["source_id"],r["instrument"]
                unique=(sid,ins,r.get("spectrum_pk",r.get("expid")))
                if unique in seen:continue
                seen.add(unique)
                err=float(r["err_kms"])
                floor=0 if survey=="desi" else (1. if "APOGEE" in ins else 5.)
                t=float(r["time_jd"])-2400000.5
                visit_groups[sid,ins,int(t)].append([t,float(r["rv_kms"]),err,floor])
    grouped=defaultdict(lambda:defaultdict(list))
    for (sid,ins,_),visits in visit_groups.items():
        t,y,s,f=np.array(visits).T
        w=s**-2
        grouped[sid][ins].append([np.average(t,weights=w),np.average(y,weights=w),
                                  np.sqrt(f[0]**2+1/w.sum())])
    results=[]
    skipped=[]
    with gzip.open(catalog,"rt") as f:
        reader=csv.DictReader(l for l in f if not l.startswith("#"))
        for r in reader:
            sid,kind=r["source_id"],r["nss_solution_type"]
            if sid not in grouped or (sid,kind) not in mass:continue
            for inst,nightly in grouped[sid].items():
                if len(nightly)<2:continue
                t,y,s=np.array(nightly).T
                try:
                    models=predictions(r,t,np.random.default_rng(int(sid)))
                except (ValueError,np.linalg.LinAlgError) as e:
                    skipped.append([sid,kind,inst,str(e)]);continue
                w=s**-2
                def stat(paths):
                    z=y-paths
                    z-=np.average(z,axis=1,weights=w)[:,None]
                    return ((z/s)**2).sum(axis=1)
                xs=stat(models)
                if kind.startswith("Orbital"):
                    xs=np.minimum(xs,stat(-models))
                dof=len(y)-1
                q=chi2.sf(xs,dof)
                null=stat(np.zeros((1,len(y))))[0]
                results.append(dict(source_id=sid,kind=kind,instrument=inst,n_nights=len(t),
                    period=float(r["period"]),m2_screen=mass[(sid,kind)],baseline_d=round(np.ptp(t),3),
                    nominal_range_kms=round(np.ptp(models[0]),3),measured_range_kms=round(np.ptp(y),3),
                    chi2_const=round(null,3),chi2_nominal=round(xs[0],3),dof=dof,
                    frac_gaia_draws_p_gt_p01=np.mean(q[1:]>.01),q95=np.quantile(q[1:],.95)))
    out=HERE/("survey_union_shapes.csv" if extension else "survey_initial_shapes.csv")
    with out.open("w") as stream:
        w=csv.DictWriter(stream,fieldnames=results[0].keys());w.writeheader()
        w.writerows(sorted(results,key=lambda r:(r["source_id"],r["kind"],r["instrument"])))
    out.with_suffix(".json").write_text(json.dumps(dict(inputs=hashes,n=len(results),skipped=skipped),indent=2)+"\n")
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("nss",type=Path)
    p.add_argument("--extension",action="store_true")
    a=p.parse_args()
    main(a.nss,a.extension)
