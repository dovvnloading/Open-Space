"""Phase-averaged detectability of the four measured control velocity curves."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import chi2, ncx2

from fit_controls import HERE, AU_KM
from sdss_sb1 import curve

BASE=HERE.parent


def read(p):
    with p.open() as stream:
        return list(csv.DictReader(stream))


def epochs(dirs):
    groups=defaultdict(list)
    hashes={}
    seen=set()
    for directory in dirs:
        for survey,name in [("LAMOST","lamost/exposures.csv"),("SDSS","epochs/sdss_visits.csv"),
                            ("DESI","epochs/desi_visits.csv")]:
            path=directory/name
            if not path.exists():
                raise FileNotFoundError(path)
            hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
            for r in read(path):
                if survey=="LAMOST":
                    t=float(r["mjd_utc_start"])+2400000.5
                    inst="LAMOST"
                    date=r["obsdate"]
                    unique=(r["obsid"],r["lmjm"])
                    err=float(r["err_kms"])
                    floor=1.
                else:
                    if int(r["quality"])!=1:continue
                    t=float(r["time_jd"])
                    inst=r["instrument"]
                    date=str(int(t-.5))
                    if survey=="SDSS":
                        unique=(r["spectrum_pk"],)
                        floor=1. if "APOGEE" in inst else 5.
                        err=float(r["err_kms"])
                    else:
                        unique=(r["targetid"],r["expid"])
                        floor=float(r["floor_kms"])
                        err=float(r["formal_err_kms"])
                key=(r["source_id"],inst,*unique)
                if key in seen:continue
                seen.add(key)
                groups[(r["source_id"],inst,date)].append((t,err,floor))
    bysource=defaultdict(lambda:defaultdict(list))
    for (sid,inst,date),vals in groups.items():
        t,s,floor=np.array(vals).T
        w=s**-2
        bysource[sid][inst].append((np.average(t,weights=w),np.sqrt(1/w.sum()+max(floor)**2)))
    return bysource,hashes


def probabilities(orbit,data):
    p,e,tp,a,ci,w=orbit[:6]
    k=2*np.pi*AU_KM*a*np.sqrt(1-ci**2)/(p*86400*np.sqrt(1-e**2))
    phase=(np.arange(1024)+.5)/1024
    start=2457389.+p*phase
    pars=np.column_stack([np.broadcast_to(v,phase.shape) for v in [p,0.,k,e,np.rad2deg(w),start]])
    delta=np.zeros(len(phase))
    dof=0
    for dates in data.values():
        if len(dates)<2:continue
        t,s=np.array(dates).T
        rv=curve(pars,t)
        rv-=np.average(rv,axis=1,weights=s**-2)[:,None]
        delta+=((rv/s)**2).sum(axis=1)
        dof+=len(t)-1
    if dof==0:return 0.,0
    threshold=chi2.isf(.001,dof)
    return float(ncx2.sf(threshold,dof,delta).mean()),dof


def main(extension):
    dirs=[BASE]+([BASE/"extension"] if extension else [])
    data,hashes=epochs(dirs)
    controls=json.loads((HERE/"control_fits.json").read_text())
    output=[]
    for sid,observations in sorted(data.items()):
        row=dict(source_id=sid,epochs=sum(map(len,observations.values())),
                 instruments="|".join(sorted(observations)))
        for name,fit in controls.items():
            row[name],row["dof"]=probabilities(fit["best_fit"],observations)
        output.append(row)
    label="union" if extension else "initial"
    path=HERE/f"cadence_{label}.csv"
    with path.open("w") as f:
        writer=csv.DictWriter(f,fieldnames=output[0].keys())
        writer.writeheader()
        writer.writerows(output)
    out=dict(input_sha256=hashes,n_targets_with_any_clean_epoch=len(output),
             n_with_within_instrument_repeat=sum(r["dof"]>0 for r in output),controls={})
    for name in controls:
        arr=np.array([r[name] for r in output])
        out["controls"][name]=dict(at_least_50pct=int((arr>=.5).sum()),at_least_90pct=int((arr>=.9).sum()),
                                   sum_of_detection_probabilities=float(arr.sum()))
    (HERE/f"cadence_{label}.json").write_text(json.dumps(out,indent=2)+"\n")
    print(out["n_targets_with_any_clean_epoch"],out["n_with_within_instrument_repeat"],out["controls"])
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--extension",action="store_true")
    main(parser.parse_args().extension)
