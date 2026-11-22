"""Report map coverage and compare numerical values to the first reduction."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main(out, earlier):
    t=pd.read_csv(out/"limits.csv",low_memory=False)
    ok=t.loc[t.status.eq("ok")].copy()
    if ok.limit_flux_cgs.isna().any() or (ok[["limit_flux_cgs","exp_s"]]<=0).any().any():
        raise ValueError("Invalid flux limit or exposure in accepted map product")
    percentiles=[0,.05,.5,.95,1]
    checks=[]
    for (m,method),g in ok.groupby(["mission","method"]):
        checks.append({"mission":m,"method":method,"n":len(g),"percentiles":percentiles,
                       "flux_limit":g.limit_flux_cgs.quantile(percentiles).tolist(),
                       "exposure_s":g.exp_s.quantile(percentiles).tolist(),
                       "counts":g.counts.quantile(percentiles).tolist()})
    compare={}
    if earlier and (earlier/"limits.csv").exists():
        old=pd.read_csv(earlier/"limits.csv",usecols=["detection_uid","mission","status","limit_flux_cgs"])
        joined=ok.merge(old,on=["detection_uid","mission"],how="left",suffixes=("","_previous"),
                        validate="one_to_one")
        for m,g in joined.groupby("mission"):
            repeated=g.loc[g.status_previous.eq("ok")]
            change=(repeated.limit_flux_cgs/repeated.limit_flux_cgs_previous-1).abs()
            compare[m]={"valid_both_runs":len(repeated),"newly_valid":len(g)-len(repeated),
                        "relative_difference_over_0p01":int((change>.01).sum()),
                        "max_relative_difference":float(change.max()) if len(change) else None}
    note={"map_distributions":checks,"comparison_to_first_reduction":compare}
    (out/"limits_qa.json").write_text(json.dumps(note,indent=2,allow_nan=False)+"\n")
    print(json.dumps(note,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("--earlier",type=Path)
    a=p.parse_args()
    main(a.out,a.earlier)
