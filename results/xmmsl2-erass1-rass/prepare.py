"""Compare measured epochs and one-sided limits, retaining failed cuts."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

MISSIONS = ("xmm", "ero", "rass")
PAIRS = (("xmm", "ero"), ("xmm", "rass"), ("ero", "rass"))
Z = 1.96


def numeric(table, key, default=np.nan):
    value = table[key] if key in table else pd.Series(default, index=table.index)
    return pd.to_numeric(value, errors="coerce")


def truth(values):
    return values.fillna(False).astype(str).str.lower().isin(["true", "1", "1.0", "yes", "t"])


def present(values):
    return values.notna() & values.astype(str).str.strip().ne("")


def prepare(out):
    cat = pd.read_csv(out/"catalog_sources.csv.gz", low_memory=False)
    if cat.detection_uid.isna().any() or cat.detection_uid.duplicated().any():
        raise ValueError("Catalogue detection keys must be present and unique")
    limits = pd.read_csv(out/"limits.csv", low_memory=False)
    if len(limits) and (limits.duplicated(["detection_uid", "mission"]).any()
                        or not limits.detection_uid.isin(cat.detection_uid).all()):
        raise ValueError("Invalid map-to-source join")
    for mission in ("ero", "rass"):
        cols = ["detection_uid", "limit_flux_cgs", "status"]
        if "model_max_relative" in limits:
            cols.append("model_max_relative")
        part = limits.loc[limits.mission.eq(mission), cols].copy()
        part = part.rename(columns={k: f"{mission}_ul_{k}" for k in cols[1:]})
        cat = cat.merge(part, how="left", on="detection_uid", validate="one_to_one")
    supplemental=truth(cat["xmm_full_only_corroborated"]) if "xmm_full_only_corroborated" in cat else pd.Series(False,index=cat.index)
    cat["selection_stratum"]=np.where(truth(cat.clean),"clean",np.where(supplemental,"corroborated full-only","other full-only"))
    admissible=truth(cat.clean) | supplemental
    quality={m:truth(cat[m+"_quality"]) for m in MISSIONS}
    quality["xmm"] |= supplemental
    records = []
    for m1,m2 in PAIRS:
        for high,low in [(m1,m2),(m2,m1)]:
            hf, he = numeric(cat,high+"_flux"), numeric(cat,high+"_err")
            lf, le = numeric(cat,low+"_flux"), numeric(cat,low+"_err")
            low_is_detection = pd.Series(True,index=cat.index) if low=="xmm" else present(cat[low+"_name"])
            ul = numeric(cat,low+"_ul_limit_flux_cgs")
            dl = low_is_detection & lf.gt(0)
            use_limit = ~low_is_detection & ul.gt(0)
            possible = hf.gt(0) & (dl | use_limit)
            low_value = lf.where(dl,ul)
            ratio = hf / low_value
            high_lo = (hf-Z*he).clip(lower=0).where(he.gt(0))
            low_hi = (lf+Z*le).where(le.gt(0)).where(dl,ul)
            guarded = high_lo / low_hi
            # Independent worst-case factors allow the spectrum to change at each epoch.
            # These factors describe only the finite set of alternative models in catalogs.py.
            hfactor = numeric(cat,high+"_model_min")
            lfactor = numeric(cat,low+"_model_max").where(dl,numeric(cat,low+"_ul_model_max_relative"))
            model = guarded * hfactor / lfactor
            good = admissible & quality[high] & (
                (dl & quality[low]) | use_limit)
            for m in (high,low):
                if m!="xmm":
                    good &= ~truth(cat[m+"_ambiguous"])
            if {high,low}=={"ero","rass"}:
                good &= ~dl | truth(cat.ero_rass_consistent)
            idx = possible & ratio.ge(1)
            rows = pd.DataFrame({
                "detection_uid": cat.detection_uid, "name": cat.name, "high": high, "low": low,
                "selection_stratum":cat.selection_stratum,
                "low_kind": np.where(use_limit,"upper_limit","detection"),
                "high_flux":hf,"high_error":he,"low_flux_or_limit":low_value,"low_error":le.where(dl),
                "ratio":ratio,"guarded_ratio":guarded,"model_floor_ratio":model,
                "quality_ok":good,"nominal_tenfold":ratio.ge(10),
                "photometric_lead":good & guarded.ge(10),
            }).loc[idx]
            records.append(rows)
    pairs = pd.concat(records,ignore_index=True)
    pairs.to_csv(out/"comparisons.csv.gz",index=False)
    leads = pairs.loc[pairs.photometric_lead].sort_values("guarded_ratio",ascending=False)
    peaks = leads.drop_duplicates("detection_uid")
    cat=cat.rename(columns={"ratio":"catalog_only_ratio","ratio_lower":"catalog_only_ratio_lower","pair":"catalog_only_pair"})
    selected = peaks.merge(cat, on=["detection_uid","name","selection_stratum"], validate="one_to_one")
    selected["id"] = selected.detection_uid
    selected["xmm_ra"], selected["xmm_dec"] = selected.ra, selected.dec
    # A good eROSITA match may provide a tighter location for counterpart checks.
    err = numeric(selected,"xmm_poserr").where(lambda s:s.gt(0),8.)
    ee = numeric(selected,"ero_poserr")
    er, ed = numeric(selected,"ero_ra"), numeric(selected,"ero_dec")
    refine = (truth(selected.ero_quality) & ~truth(selected.ero_ambiguous)
              & ee.gt(0) & ee.lt(err) & er.between(0,360) & ed.between(-90,90))
    selected["ra"] = selected.ra.where(~refine,er)
    selected["dec"] = selected.dec.where(~refine,ed)
    selected["poserr_arcsec"] = err.where(~refine,ee)
    selected["position_mission"] = np.where(refine,"ero","xmm")
    selected.to_csv(out/"input.csv",index=False)
    print(f"{len(cat)} slew rows, {len(pairs)} comparable pair directions, {len(leads)} qualifying comparisons in {len(peaks)} slew rows",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    prepare(p.parse_args().out)
