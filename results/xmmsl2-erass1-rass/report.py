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
    legacy_file = out/"rass_legacy.csv"
    legacy_confusion = set()
    if legacy_file.exists():
        legacy=pd.read_csv(legacy_file,low_memory=False)
        if len(legacy):
            plausible=legacy.plausible.astype(str).str.lower().isin(["true","1"])
            unresolved=legacy.needs_review.astype(str).str.lower().isin(["true","1"])
            legacy_confusion=set(legacy.loc[plausible | unresolved,"detection_uid"])
    legacy_pair=ok.low.eq("rass") & ok.low_kind.eq("upper_limit") & ok.detection_uid.isin(legacy_confusion)
    robust = ok.loc[ok.model_floor_ratio.ge(10) & ~legacy_pair].sort_values("model_floor_ratio",ascending=False)
    robust = robust.drop_duplicates("detection_uid").set_index("detection_uid")
    for k in ("high","low","low_kind","ratio","guarded_ratio","model_floor_ratio"):
        register["robust_"+k] = register.detection_uid.map(robust[k])
    register["legacy_rosat_neighbor"] = register.detection_uid.isin(legacy_confusion)
    register["disposition"] = register.vet_status
    mask=register.vet_status.eq("candidate") & ~register.detection_uid.isin(robust.index)
    register.loc[mask,"disposition"]="inconclusive"
    register.loc[mask,"rationale"]="Tenfold change lacks a surviving spectral stress test or relies on a disputed ROSAT upper limit"
    decision_file=out/"event_history"/"published_tde_matches.csv"
    if decision_file.exists():
        known=pd.read_csv(decision_file)
        if len(known):
            source_radius=register.set_index("detection_uid").match_radius_arcsec
            known=known.loc[known.separation_arcsec.le(known.detection_uid.map(source_radius))]
            notes=known.groupby("detection_uid").published_name.agg(lambda g:" | ".join(sorted(set(g))))
            register["published_tde"] = register.detection_uid.map(notes)
            mask=register.disposition.eq("candidate") & register.published_tde.notna()
            register.loc[mask,"disposition"]="inconclusive"
            register.loc[mask,"rationale"]="Published TDE in the positional error circle; separate outburst unestablished"
    controls_file=out/"positive_controls.csv"
    if controls_file.exists():
        known=pd.read_csv(controls_file).dropna(subset=["detection_uid"]).drop_duplicates("detection_uid").set_index("detection_uid")
        register["known_slew_flare_reference"]=register.detection_uid.map(known.bibcode)
        flag=register.known_slew_flare_reference.notna()
        register.loc[flag,"disposition"]="previously_reported_flare"
        register.loc[flag,"rationale"]="This position belongs to a published slew-flare source in positive_controls.csv"
    field_file=out/"published_hosts"/"local_galaxies.csv"
    if field_file.exists():
        galaxies=pd.read_csv(field_file,dtype={"ls_id":str})
        first=galaxies.loc[galaxies["rank"].eq(1)].set_index("detection_uid")
        for k in ("ls_id","type","host_ra","host_dec","sep_arcsec","pcc","shifted_at_least_as_good","shifted_total","empirical_fraction"):
            register["local_galaxy_"+k]=register.detection_uid.map(first[k])
    salvato=out/"published_hosts"/"LS10.csv.gz"
    if salvato.exists() and salvato.stat().st_size>50:
        prob=pd.read_csv(salvato,dtype=str,low_memory=False)
        if {"pub_NWAY_p_i","pub_NWAY_p_any"}<=set(prob):
            prob["sort_probability"]=pd.to_numeric(prob.pub_NWAY_p_i,errors="coerce")
            best=prob.sort_values("sort_probability",ascending=False).drop_duplicates("ero_name").set_index("ero_name")
            for k in ("pub_LS10_FULLID","pub_LS10_RA","pub_LS10_DEC","pub_LS10_TYPE","pub_NWAY_p_i","pub_NWAY_p_any","pub_bestCTP"):
                if k in best:
                    register["salvato_"+k[4:]]=register.ero_name.map(best[k])
    cases_file=out/"casework"/"case_index.csv"
    if cases_file.exists() and cases_file.stat().st_size>4:
        cases=pd.read_csv(cases_file).set_index("detection_uid")
        for k in ("interpretation","prior_event_leads","publications"):
            register["case_"+k]=register.detection_uid.map(cases[k])
        manual=register.disposition.eq("candidate") & (register.case_prior_event_leads.fillna(0).gt(0) | register.case_interpretation.fillna("").str.contains("bright-star alternative"))
        register.loc[manual,"disposition"]="inconclusive"
        register.loc[manual,"rationale"]="Source-specific publication or optical-loading alternative in casebook requires review"
    # Flat comparison sheet for every qualifying observation. Keep the
    # separate catalogue/limit status so a blank can never pass for zero flux.
    display = ["detection_uid","source_group","name","ra","dec","position_mission",
               "xmm_time","xmm_flux","xmm_err","ero_name","ero_flux","ero_err",
               "ero_ul_limit_flux_cgs","ero_ul_status","rass_name","rass_flux","rass_err",
               "rass_ul_limit_flux_cgs","rass_ul_status","high","low","low_kind",
               "ratio","guarded_ratio","model_floor_ratio","robust_model_floor_ratio",
               "host_catalog","host_id","host_type","host_offset_arcsec","match_radius_arcsec",
               "salvato_NWAY_p_i","local_galaxy_pcc","local_galaxy_empirical_fraction",
               "known_slew_flare_reference","published_tde","case_publications","prior_bibcodes",
               "event_review","stellar_evidence","service_outcomes","disposition","rationale"]
    register.sort_values("guarded_ratio",ascending=False)[
        [k for k in display if k in register]
    ].to_csv(out/"review_photometry.csv",index=False)
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
    selected_cols=["detection_uid","name","ra","dec","high","low","low_kind","ratio","guarded_ratio","model_floor_ratio","robust_model_floor_ratio","host_id","host_type","xmm_counterpart","xmm_xtyp","known_associations","rationale","disposition","case_interpretation","case_publications","salvato_NWAY_p_any","salvato_NWAY_p_i","local_galaxy_pcc","local_galaxy_empirical_fraction"]
    examples=register.sort_values("guarded_ratio",ascending=False).drop_duplicates("source_group")
    examples[[c for c in selected_cols if c in examples]].to_csv(out/"review_examples.csv",index=False)
    missing_patterns={}
    for encoded in register.service_outcomes:
        outcomes=json.loads(encoded)
        missing=",".join(sorted(k for k,v in outcomes.items() if v!="ok")) or "complete"
        missing_patterns[missing]=missing_patterns.get(missing,0)+1
    snapshots={}
    ex=examples.fillna("")
    for label,part in [("candidates",ex.loc[ex.disposition.eq("candidate")]),("optical_hosts",ex.loc[ex.host_type.astype(str).str.contains(r"Galaxy|Sy[12]|AGN|QSO|^G$|XRB|HXB|LXB",case=False,regex=True)]),("brightest_unresolved",ex.loc[ex.disposition.eq("inconclusive")]),("known_or_stellar",ex.loc[ex.disposition.eq("known_or_stellar")])]:
        display_cols=selected_cols+["service_outcomes","event_review","prior_bibcodes","xmm_flux","ero_flux","rass_flux","ero_ul_limit_flux_cgs","rass_ul_limit_flux_cgs"]
        snapshots[label]=part[[c for c in display_cols if c in part]].head(15).to_dict("records")
    (out/"review_shortlist.json").write_text(json.dumps(snapshots,indent=2,default=str)+"\n")
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
        "dispositions":frequencies(register.disposition),"selection_strata":frequencies(register.selection_stratum),
        "legacy_1rxs_positional_warnings":int(register.legacy_rosat_neighbor.sum()),"service_completeness":services,
        "candidate_detections":len(candidates),"candidate_source_designations":int(candidates.source_group.nunique()),
        "xray_shifted":{k:{c:v for c,v in d.items() if c!="offsets"} for k,d in json_load("shifted_controls.json").items()},
        "host_shifted":shifts,"missing_check_combinations":missing_patterns,
    }
    (out/"summary.json").write_text(json.dumps(summ,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summ, indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    summarize(p.parse_args().out)
