"""Publish the final per-detection register while retaining automated screening decisions."""
import argparse
import json
from pathlib import Path

import pandas as pd


def annotate(out):
    study=out.parent
    reg=pd.read_csv(out/"review_register.csv.gz",low_memory=False).fillna("")
    pre=reg.set_index("detection_uid").disposition
    review=pd.read_csv(study/"published_event_review.csv",keep_default_na=False)
    review["detection_uid"]=review.reviewed_detection_uids.str.split(r"\s*\|\s*")
    individual=review.explode("detection_uid")
    if individual.detection_uid.duplicated().any() or not individual.detection_uid.isin(reg.detection_uid).all():
        raise ValueError("A reviewed prior event has an ambiguous or missing detection key")
    individual=individual.set_index("detection_uid")
    for column in ("prior_identifier","already_published_phenomenon","same_slew_high_state",
                   "source_specific_evidence","primary_reference","additional_reference"):
        reg[column]=reg.detection_uid.map(individual[column]).fillna("")
    reg["screening_disposition"]=reg.disposition
    same=reg.same_slew_high_state.eq("yes")
    reg.loc[same,"disposition"]="previously_reported_flare"
    reg.loc[same,"rationale"]=reg.loc[same,"source_specific_evidence"]
    reg["named_in_li_slew_variability_table"]=[
        f"J/MNRAS/512/3858/table {name}" in refs for name,refs in zip(reg.name,reg.event_review)]
    original=reg.named_in_li_slew_variability_table & reg.high.eq("xmm") & reg.disposition.eq("inconclusive")
    reg.loc[original,"disposition"]="previously_catalogued_slew_variable"
    reg.loc[original,"rationale"]=("Named in the complete Li et al. XMMSL2 variability table. "
          "The new eROSITA flux/limit is a later epoch. A separate new outburst is unestablished.")
    reg.loc[original,"primary_reference"]="https://doi.org/10.1093/mnras/stac714"
    # Published source identity alone cannot decide whether a later flux change
    # has previously appeared in the literature.
    optional=study/"source_decisions.csv"
    if optional.exists():
        extra=pd.read_csv(optional,keep_default_na=False)
        required={"detection_uid","source","decision","evidence","reference"}
        if not required.issubset(extra) or extra.detection_uid.duplicated().any():
            raise ValueError("Invalid source-level decision table")
        for row in extra.itertuples():
            flag=reg.detection_uid.eq(row.detection_uid)
            if flag.sum()!=1 or reg.loc[flag,"name"].iloc[0]!=row.source or not row.evidence or not row.reference:
                raise ValueError(f"Untraceable individual review: {row.detection_uid}")
            if row.decision not in {"candidate","inconclusive","previously_reported_flare","known_or_stellar"}:
                raise ValueError("Unknown source decision")
            reg.loc[flag,["disposition","rationale","decision_reference"]]=[row.decision,row.evidence,row.reference]
    host_file=out/"host_adjudication.csv"
    field_cols=[]
    if host_file.exists():
        fields=pd.read_csv(host_file,keep_default_na=False).set_index("detection_uid")
        if fields.index.duplicated().any() or set(fields.index)!=set(reg.detection_uid):
            raise ValueError("Host adjudications do not cover the photometric selection")
        for key in ("host","galaxy_support","host_high_sep_arcsec","NWAY_p_any","NWAY_p_i",
                    "Gaia_astrometric_or_variable_in_high95","other_eROSITA_main_in_high95"):
            if key in fields:
                reg["high_field_"+key]=reg.detection_uid.map(fields[key])
                field_cols.append("high_field_"+key)
    broker=json.loads((out/"broker_cones.json").read_text())
    cones={r["detection_uid"]:r for r in broker["fields"] if r.get("detection_uid")}
    colours=pd.read_csv(out/"rass_hardness/source_diagnostics.csv",low_memory=False)
    for i,r in reg.loc[reg.disposition.eq("candidate")].iterrows():
        cone=cones.get(r.detection_uid,{})
        loci=cone.get("loci",[])
        gate=float(r.match_radius_arcsec)
        in_gate=[a for a in loci if float(a["sep_arcsec"])<=gate]
        core=json.loads(r.service_outcomes)
        source_checks=("Gaia","Legacy","NED","SIMBAD","SIMBAD_bibliography","VSX","WISE","papers")
        measured=colours.loc[colours.xmmsl2_name.eq(r["name"]) & colours.low_mission.eq(r.low)]
        floor=pd.to_numeric(measured.min_rass_only_guarded_ratio,errors="coerce").max()
        valid=(all(core.get(c)=="ok" for c in source_checks) and not r.stellar_evidence
               and not r.named_in_li_slew_variability_table and float(r.robust_model_floor_ratio)>=10
               and str(r.high_field_galaxy_support).startswith("published primary +")
               and float(r.high_field_Gaia_astrometric_or_variable_in_high95)==0
               and float(r.high_field_other_eROSITA_main_in_high95)==0 and floor>=10
               and cone.get("status")=="ok" and in_gate and broker.get("controls_recovered")
               and not any(l.get("association_marked") or l.get("tns_records") for l in loci))
        if not valid:
            raise ValueError(f"Independent candidate evidence has changed for {r.detection_uid}; review required")
        reg.loc[i,"rass_colour_guarded_floor_fixed_faint_ecf"]=floor
        reg.loc[i,"tns_broker_loci_in_35arcsec"]=len(loci)
        reg.loc[i,"tns_broker_nearest_locus_arcsec"]=min(l["sep_arcsec"] for l in in_gate)
        reg.loc[i,"tns_broker_associated_records"]=0
        reg.loc[i,"novelty_scope"]="No report of this specific decline identified in reviewed sources. Complete official TNS census unavailable."
    field_cols.extend(k for k in ("rass_colour_guarded_floor_fixed_faint_ecf","tns_broker_loci_in_35arcsec",
                     "tns_broker_nearest_locus_arcsec","tns_broker_associated_records","novelty_scope") if k in reg)
    compact=pd.read_csv(out/"review_photometry.csv",keep_default_na=False,dtype=str)
    result=reg.set_index("detection_uid")
    for k in ["named_in_li_slew_variability_table","screening_disposition","disposition","rationale","prior_identifier","source_specific_evidence","primary_reference","decision_reference"]+field_cols:
        if k in result:
            compact[k]=compact.detection_uid.map(result[k])
    if len(compact)!=len(result) or set(compact.detection_uid)!=set(result.index):
        raise ValueError("Photometry and adjudications do not cover the same observations")
    compact.to_csv(out/"final_source_register.csv",index=False)
    positive=reg.loc[reg.disposition.eq("candidate")]
    positive.to_csv(out/"candidates.csv",index=False)
    counts=reg.disposition.value_counts().to_dict()
    summary={"photometric_lead_detections":len(reg),"distinct_slew_designations":reg.source_group.nunique(),
             "dispositions":counts,"distinct_designations_per_disposition":reg.groupby("disposition").source_group.nunique().to_dict(),
             "candidates":len(positive),"candidate_designations":int(positive.source_group.nunique()),
             "changes_from_automated_screen":int((reg.disposition!=reg.detection_uid.map(pre)).sum()),
             "scope":"A published flare designation can repeat across slew visits. Per-category source-name counts may overlap. Competing hosts, unsuccessful searches, and prior host identifications retain their evidence in the discovery register."}
    required={"Gaia","Legacy","NED","SIMBAD","SIMBAD_bibliography","TNS","VSX","WISE","papers"}
    complete=reg.service_outcomes.map(lambda text: all(json.loads(text).get(k)=="ok" for k in required))
    summary["named_in_li_variability_table"]=int(reg.named_in_li_slew_variability_table.sum())
    summary["complete_core_vetting_detections"]=int(complete.sum())
    summary["spectrally_robust_and_complete_core_vetting"]=int((complete & pd.to_numeric(reg.robust_model_floor_ratio,errors="coerce").ge(10)).sum())
    summary["only_incompleteness_in_automated_rationale"]=int(reg.rationale.str.fullmatch(r"incomplete [^|]+",na=False).sum())
    (out/"final_summary.json").write_text(json.dumps(summary,indent=2,default=int)+"\n")
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    print(json.dumps(annotate(p.parse_args().out),default=int))
