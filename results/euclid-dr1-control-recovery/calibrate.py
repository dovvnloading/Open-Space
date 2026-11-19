#!/usr/bin/env python3
"""Per-control calibrated empirical receipt (conditional to catalogue matches).
A Wilson interval is not survey/injection completeness. Resubstitution warning.
"""
import json,math
from pathlib import Path
import pyarrow as pa
import selection
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"

def wilson(x,n,z=1.95996398454):
    p=x/n
    d=1+z*z/n
    c=(p+z*z/(2*n))/d
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [round(c-h,5),round(c+h,5)]

def main():
    audit=json.loads((OUT/"control_audit.json").read_text())
    freeze=json.loads((ROOT/"freeze.json").read_text())
    if freeze["source_index_sha256"]!="f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a":
        raise AssertionError("Changed source index")
    result=[]
    for r in audit:
        x=dict(name=r["name"],class_label=r["class"],z=r["z"],matched="source" in r,
               passed=False,old_strict_recovered=r.get("original_strict_pass",False),
               cut_failures=r.get("failed_original_cuts",[]))
        if "source" in r:
            # Independent read: one-source Parquet replay through the identical
            # array mask used in streaming scan.
            tab=pa.Table.from_pylist([r["source"]])
            arr=selection.arrays(tab)
            x["passed"]=bool(arr["eligible"][0])
            x["new_score"]=round(float(arr["score"][0]),5)
            x["YJ_used"]=float(arr["YJ_USED"][0]) if math.isfinite(arr["YJ_USED"][0]) else None
            x["VIS_SNR"]=r["derived"]["VIS_SNR_3fwhm"]
            x["VIS_DET"]=r["source"]["vis_det"]
            x["flag_j"]=r["source"]["flag_j"]
            x["flag_h"]=r["source"]["flag_h"]
        result.append(x)
    q=[r for r in result if r["class_label"]=="quasar" and r["matched"]]
    d=[r for r in result if r["class_label"]=="dwarf" and r["matched"]]
    kq=sum(r["passed"] for r in q);kd=sum(r["passed"] for r in d)
    if len(q)!=24 or len(d)!=14 or kq!=21 or kd!=5:
        raise AssertionError("Frozen control operating point diverged: "+str((len(q),len(d),kq,kd)))
    receipt={
      "status":"MATCHED_CONTROL_APPEARANCE_REPLAY_PASS",
      "selection_id":freeze["selection_id"],
      "matched_quasar_recovery":dict(numerator=kq,denominator=len(q),fraction=kq/len(q),
                                    wilson95=wilson(kq,len(q))),
      "matched_dwarf_retention":dict(numerator=kd,denominator=len(d),fraction=kd/len(d),
                                     wilson95=wilson(kd,len(d))),
      "new_z_gt7_matched_quasars":dict(numerator=sum(r["passed"] for r in q if r["z"]>7),
                                      denominator=sum(r["z"]>7 for r in q)),
      "baseline":dict(quasar_recovery="3/24",borrowed_NISP_shape="1/24"),
      "tuned_using_same_labels":True, "independent_holdout":False,
      "footprint_and_flux_complete":False,
      "completeness_calibrated_for_field":False,
      "dwarf_population_purity_calibrated":False,
      "spectroscopic_discoveries":0,
      "interpretation":"Apparent recovery on the matched control cohort, with Wilson binomial uncertainty; not an imaging-area selection function."
    }
    OUT.mkdir(exist_ok=True)
    (OUT/"selection_control_recovery.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    (OUT/"calibration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("FROZEN_CONTROL_REPLAY",json.dumps(receipt),flush=True)
if __name__=="__main__":main()
