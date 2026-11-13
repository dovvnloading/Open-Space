"""Post-screen comparison of fixed periodic photometric models on held-out cadences."""

import argparse
import json
from pathlib import Path

from .tess import flatten, load_product
from .tess_archive import atomic_bytes


def evaluate(train, test, ephem):
    """Compare predictive sums of squared standardized residuals to a held-out mean.

    The period, phase, stellar harmonic coefficients and box shape are fit or
    selected in the discovery sector only. A single DC offset is allowed per
    test model because the sector flux scales were normalized separately.
    """
    import numpy as np

    period, epoch, duration = (ephem[k] for k in ("period_days", "epoch_btjd", "duration_days"))
    def basis(t, model):
        angle=2*np.pi*(t-epoch)/period
        if model=="box":
            phase=(t-epoch+period/2)%period-period/2
            return np.column_stack((np.ones(len(t)), (abs(phase)<=duration/2).astype(float)))
        n=int(model.split("_")[-1])
        columns=[np.ones(len(t))]
        for k in range(1,n+1):
            columns.extend([np.sin(angle*k),np.cos(angle*k)])
        return np.column_stack(columns)
    t1, f1, err1 = train
    t2, f2, err2 = test
    w2=err2**-2
    def heldout_loss(predicted):
        delta=f2-predicted
        return float(np.sum(w2*(delta-np.average(delta, weights=w2))**2))
    flat_loss=heldout_loss(0)
    outcomes=[]
    for name in ("box","fourier_1","fourier_4"):
        design=basis(t1,name)
        coeff=np.linalg.lstsq(design/err1[:,None], f1/err1,rcond=None)[0]
        predicted=basis(t2,name) @ coeff
        loss=heldout_loss(predicted)
        harmonic=None
        if name.startswith("fourier"):
            harmonic=int(np.argmax(np.hypot(coeff[1::2],coeff[2::2])))+1
        outcomes.append({"shape":name,"fitted_coefficients":len(coeff),
                         "dominant_harmonic":harmonic,
                         "improvement_over_flat":flat_loss-loss})
    return {"baseline_test_loss":flat_loss,"test_samples":len(t2), "models":outcomes,
            "best_test_shape":max([{"shape":"flat","improvement_over_flat":0}]+outcomes,
                                  key=lambda r:r["improvement_over_flat"])["shape"]}


def run(plan, evidence, cache):
    targets={t["tic_id"]:t for t in plan["targets"]}
    results=[]
    for row in evidence["targets"]:
        if row["disposition"] not in {"two_sector_periodic_dip", "periodic_dip_with_veto"}:
            continue
        tic=row["tic_id"]
        target=targets[tic]
        curves=[]
        for p, saved in zip(target["products"],row["sources"],strict=True):
            raw,meta=load_product(p,cache,tic)
            if meta["sha256"]!=saved["sha256"]:
                raise ValueError(f"Input digest changed for {tic}")
            curves.append(flatten(raw))
        results.append({"tic_id":tic,"role":row["role"],"disposition":row["disposition"],
                        "period_days":row["ephemeris"]["period_days"],
                        **evaluate(*curves,row["ephemeris"])})
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair",nargs=2,type=Path,action="append",required=True,metavar=("PLAN","EVIDENCE"))
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from hashlib import sha256
    cases, inputs=[],[]
    for plan,evidence in args.pair:
        packet=json.loads(evidence.read_bytes())
        plan_hash=sha256(plan.read_bytes()).hexdigest()
        if packet["plan_sha256"] != plan_hash:
            raise ValueError("Evidence differs from plan")
        inputs.append({"plan_sha256":plan_hash,"evidence_sha256":sha256(evidence.read_bytes()).hexdigest()})
        cases.extend(run(json.loads(plan.read_bytes()),packet,args.cache))
    atomic_bytes(args.output,(json.dumps({"experiment":"post_screen_shape_check_v1","inputs":inputs,
                                          "cases":cases},indent=2,allow_nan=False)+"\n").encode())


if __name__=="__main__":
    main()
