"""Make a source-by-source evidence register and the counts used in the README."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def summarize(out):
    cat = pd.read_csv(out/"catalog_sources.csv.gz",low_memory=False)
    leads = pd.read_csv(out/"input.csv",low_memory=False)
    vet = pd.read_csv(out/"vetting/vetting.csv",low_memory=False)
    pairs = pd.read_csv(out/"comparisons.csv.gz",low_memory=False)
    limits = pd.read_csv(out/"limits.csv",low_memory=False)
    register = leads.merge(vet.drop(columns=["ra","dec"]),on="id",validate="one_to_one",how="outer",indicator=True)
    if not register._merge.eq("both").all():
        raise ValueError("Each photometric lead must have a source-level vetting decision")
    register.drop(columns=["_merge"],inplace=True)
    # Use any qualifying pair from this detection, including a runner-up pair
    # whose result survives an independent spectral stress test.
    ok = pairs.loc[pairs.photometric_lead].copy()
    robust = ok.loc[ok.model_floor_ratio.ge(10)].sort_values("model_floor_ratio",ascending=False)
    robust = robust.drop_duplicates("detection_uid").set_index("detection_uid")
    for k in ("high","low","low_kind","ratio","guarded_ratio","model_floor_ratio"):
        register["robust_"+k] = register.detection_uid.map(robust[k])
    register["disposition"] = register.vet_status
    mask=register.vet_status.eq("candidate") & ~register.detection_uid.isin(robust.index)
    register.loc[mask,"disposition"]="inconclusive"
    register.loc[mask,"rationale"]="Tenfold change depends on the selected band conversion or lacks a spectral stress test"
    register.to_csv(out/"review_register.csv.gz",index=False)
    candidates=register.loc[register.disposition.eq("candidate")].sort_values("robust_model_floor_ratio",ascending=False)
    candidates.to_csv(out/"candidates.csv",index=False)
    controls=pd.read_csv(out/"vetting/host_shifted_controls.csv") if (out/"vetting/host_shifted_controls.csv").exists() else pd.DataFrame()
    shifts={}
    if len(controls):
        for k,g in controls.groupby("catalog"):
            obs=g.real_with_counterpart.sum()
            ng=g.shifted_checked.sum()
            nr=g.real_checked.sum()
            ns=g.shifted_with_counterpart.sum()
            shifts[k]={"real_matches":int(obs),"real_checked":int(nr),"shifted_matches":int(ns),
                       "shifted_checked":int(ng),"chance_per_position":float(ns/ng) if ng else None,
                       "approx_false_fraction":float((ns/ng)/(obs/nr)) if obs and nr and ng else None}
    services={}
    for encoded in register.service_outcomes:
        for key,value in json.loads(encoded).items():
            bucket=services.setdefault(key,{"ok":0,"incomplete":0})
            bucket["ok" if value=="ok" else "incomplete"]+=1
    def frequencies(values):
        return {str(k):int(v) for k,v in values.value_counts().items()}
    def json_load(name):
        return json.loads((out/name).read_text())
    summ={
        "catalogues":json_load("catalog_statistics.json"),
        "model":json_load("spectral_model.json"),"map_status":json_load("limits_summary.json")["counts"],
        "positive_comparisons":len(pairs),"quality_comparisons":int(pairs.quality_ok.sum()),
        "compared_slew_rows":int(pairs.detection_uid.nunique()),
        "quality_compared_slew_rows":int(pairs.loc[pairs.quality_ok,"detection_uid"].nunique()),
        "nominal_tenfold":int(pairs.loc[pairs.nominal_tenfold,"detection_uid"].nunique()),
        "quality_nominal_tenfold":int(pairs.loc[pairs.nominal_tenfold & pairs.quality_ok,"detection_uid"].nunique()),
        "qualified_pairs":len(ok),"photometric_lead_detections":len(register),
        "photometric_lead_source_designations":int(register.source_group.nunique()),
        "model_stress_test_leads":len(robust),
        "photometric_pairs":frequencies(ok.high+" over "+ok.low+" "+ok.low_kind),
        "dispositions":frequencies(register.disposition),"service_completeness":services,
        "candidate_detections":len(candidates),"candidate_source_designations":int(candidates.source_group.nunique()),
        "xray_shifted":{k:{c:v for c,v in d.items() if c!="offsets"} for k,d in json_load("shifted_controls.json").items()},
        "host_shifted":shifts,
    }
    (out/"summary.json").write_text(json.dumps(summ,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summ, indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    summarize(p.parse_args().out)
