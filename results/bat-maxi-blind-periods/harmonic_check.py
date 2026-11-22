"""Test whether a longer published orbit predicts unequal candidate-length cycles."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from field_checks import load_archives
from search import bounded_uncertainties


def model(t,period,orders):
    theta=2*np.pi*(t-55000.)/period
    return np.column_stack([np.ones(len(t))]+[v for h in orders for v in (np.sin(theta*h),np.cos(theta*h))])


def predict_gain(t,y,e,period,order,training):
    # Both models contain the strongest three multiples of the target frequency.
    # The second permits other harmonics of the longer orbit. Fit only the training half.
    small=[order,2*order,3*order]
    large=sorted(set(small)|set(range(1,9)))
    test=~training
    w=1/e**2
    def loss(orders):
        x=model(t,period,orders)
        b=np.linalg.lstsq(x[training]*np.sqrt(w[training,None]),y[training]*np.sqrt(w[training]),rcond=None)[0]
        residual=y[test]-x[test]@b
        # Permit a new mean rate in the held-out interval equally for both models.
        residual-=np.average(residual,weights=w[test])
        return float(np.sum(w[test]*residual**2))
    l0,l1=loss(small),loss(large)
    return l0,l1,l0-l1


def main(folder,out):
    decisions=json.loads((folder/"decisions.json").read_text())
    selected=[r for r in decisions if r.get("harmonic_order") and r.get("diagnostic_type")=="orbital_harmonic"]
    keys={r["source_id"]+"."+b for r in selected for b in ("bat","maxi")}
    curves=load_archives(sorted((folder/"evidence").glob("curves-*.tar.gz")),keys)
    rows=[]
    for r in selected:
        h=int(r["harmonic_order"])
        reference=next(ref for ref in r["references"] if ref.get("role")=="orbital_fundamental")
        for name,p in (("published orbit",reference["period_days"]),("candidate period times harmonic order",r["period_days"]*h)):
            for band in ("bat","maxi"):
                t,y,err=curves[r["source_id"]+"."+band]
                e,_=bounded_uncertainties(t,err)
                half=t<=(t[0]+t[-1])/2
                for label,train in (("predict later",half),("predict earlier",~half)):
                    loss0,loss1,gain=predict_gain(t,y,e,p,h,train)
                    rows.append(dict(source=r["source_id"],instrument=band,heldout=label,period_choice=name,
                                     parent_period_days=p,target_period_days=r["period_days"],harmonic_order=h,
                                     target_family_loss=loss0,parent_family_loss=loss1,
                                     heldout_loss_reduction=gain,relative_reduction=gain/loss0,
                                     reference=reference["url"]))
    out.mkdir(parents=True,exist_ok=True)
    if rows:
        with (out/"orbital-harmonics.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    main(a.folder,a.out)
