#!/usr/bin/env python3
"""Merge primary HEASARC QPE receipts with conservative source-product retries."""
import datetime as dt
import json
from collections import Counter
from heasarc_batch import OUT, manifest, save

FULL="position_verified_net_pn_soft_screened"
PARTIAL="position_verified_pn_single_exposure_net_soft_screened_PARTIAL"

def main():
    cohort,fr=manifest()
    expected={(oid,r["detid"]) for oid,rows in cohort.items() for r in rows}
    initial={}; full={}; partial={}; ids=set(); batch_info=[]; errors=[]
    for file in sorted(OUT.glob("heasarc_batch_cohort_*.json")):
        r=json.loads(file.read_text())
        if r.get("cohort_manifest_sha256")!=fr["final_csv_sha256"]:
            errors.append("cohort sha mismatch "+file.name);continue
        batch_info.append({"file":file.name,"stage":"primary","attempted":len(r.get("observations",[])),
                           "status":r.get("status")})
        for o in r.get("observations",[]):
            if o["obs_id"] in ids:errors.append("duplicate obsid "+o["obs_id"])
            ids.add(o["obs_id"])
            for s in o.get("source_results",[]):
                k=(o["obs_id"],s.get("detid"))
                if k not in expected:errors.append("foreign primary row "+repr(k))
                elif s.get("status")==FULL:initial[k]=full[k]=s
    for file in sorted(OUT.glob("heasarc_rescue_*.json")):
        r=json.loads(file.read_text())
        if r.get("manifest_sha256")!=fr["final_csv_sha256"]:
            errors.append("rescue sha mismatch "+file.name);continue
        batch_info.append({"file":file.name,"stage":"retry","attempted":len(r.get("observations",[])),
                           "status":r.get("status")})
        for o in r.get("observations",[]):
            for s in o.get("source_results",[])+o.get("partial_source_results",[]):
                k=(o["obs_id"],s.get("detid"))
                if k not in expected:errors.append("foreign retry row "+repr(k));continue
                if s.get("status")==FULL:full.setdefault(k,s)
                if s.get("status")==PARTIAL and k not in full:partial.setdefault(k,s)
    for k in full:partial.pop(k,None)
    flags={}; bg=Counter(); surrogates=0; recoveries=0; bgflag=0
    for k,s in list(full.items())+list(partial.items()):
        p=s.get("independent_field_particle_background") or {}
        bg[p.get("status","not_available")]+=1
        b=s.get("offsource_scaled_background_control") or {}
        bgflag+=int(b.get("recurrence_unvetted") is True)
        inj=s.get("poisson_count_domain_injection") or {}
        surrogates+=int(inj.get("n_trials",0));recoveries+=int(inj.get("n_recovered",0))
        if s.get("unvetted_repeated_peak_trigger"):
            f={"obs_id":k[0],"detid":k[1],"srcid":s.get("srcid"),"ra":s.get("ra"),
               "dec":s.get("dec"),"source_product":s.get("pn_product"),
               "product_sha256":s.get("merged_pn_fits_sha256",s.get("pn_fits_sha256")),
               "coverage":"full_merged" if k in full else "partial_per_exposure",
               "finder":s.get("finder"),"field_particle":p,"scaled_background":b,
               "coincident_particle_flare":any(
                 x.get("source_masked_particle_flare_coincident") is True for x in p.get("peak_checks",[])),
               "validated_new_qpe":False,"classification":"UNVETTED_FLAG; no complete independent soft-event-region, pileup and novelty gates"}
            if k==("0144570101","101445701010002"):
                f["classification"]="KNOWN_SEYFERT_2_AGN + PARTICLE FLARE VETO; NOT NEW QPE"
                f["name"]="XMMU J134736.6+173403; SDSS J134736.39+173404.6"
                f["literature"]=["https://arxiv.org/abs/0801.4246","https://arxiv.org/abs/1803.10976"]
            flags[k]=f
    covered=set(full)|set(partial)
    summary={
      "as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
      "cohort_sha256":fr["final_csv_sha256"],
      "frozen_obsids":670,"frozen_detections":967,
      "obsids_attempted_by_primary_batches":len(ids),
      "all_frozen_obsids_attempted":len(ids)==670,
      "original_full_pn_source_curves":len(initial),
      "new_full_curves_recovered":len(full)-len(initial),
      "verified_full_merged_pn_source_curves":len(full),
      "partial_single_exposure_pn_source_curves":len(partial),
      "unique_frozen_detection_rows_with_pn_data":len(covered),
      "obsids_with_full_merged_pn_data":len({k[0] for k in full}),
      "obsids_with_any_pn_data":len({k[0] for k in covered}),
      "remaining_unscreened_detection_rows":967-len(covered),
      "remaining_obsids_without_any_pn_data":670-len({k[0] for k in covered}),
      "field_particle_control_statuses":dict(bg),
      "scaled_offsource_aperture_repeat_flags":bgflag,
      "count_domain_poisson_surrogate_trials":surrogates,
      "count_domain_poisson_surrogate_recoveries":recoveries,
      "preliminary_recurrent_triggers":len(flags),"triggers":list(flags.values()),
      "validated_new_qpes":0,"remarkable_allowed":False,
      "SAS_pileup_complete":False,"independent_soft_event_region_complete":False,
      "event_list_injection_complete":False,"trial_adjusted_false_positive_rate":None,
      "complete_SIMBAD_4LAC_legacy_TNS_ATel_multiepoch":False,
      "ready_for_review":False,"merge_allowed":False,
      "missing_frozen_rows":[{"obs_id":o,"detid":d} for o,d in sorted(expected-covered)],
      "source_receipts":batch_info,"integrity_issues":errors,
      "referee":"KEEP DRAFT; source-masked field particle and source BACK1-3 are diagnostics, not fully validated independent soft-background or event-list controls; incomplete host and novelty vetting; no confirmed QPE."
    }
    save("heasarc_postrescue_referee.json",summary)
    txt=["# HEASARC 4XMM-DR14 pn full and partial source-coverage referee","",
      f"**Frozen sources:** 967 detections in 670 observations; **primary attempted:** {len(ids)} obsIDs.",
      f"**Full merged pn light curves:** {len(full)} ({len(full)-len(initial)} recovered).",
      f"**Additional partial per-exposure pn light curves:** {len(partial)}; these cannot stand for complete observation coverage.",
      f"**ObsIDs with some pn:** {len({k[0] for k in covered})}; **still missing:** {967-len(covered)} detection rows.",
      f"**Triggers:** {len(flags)} preliminary, **0 confirmed novel QPEs**.",
      f"**Background scaled aperture recurring flags:** {bgflag}. **Poisson count-domain surrogate pair recoveries:** {recoveries}/{surrogates}. Not event-list sensitivity or false-alarm rate.","",
      "## Individual triggers",""]
    for f in flags.values():
        txt.extend([f"### {f['obs_id']} / {f['detid']}",
           f"Coverage {f['coverage']}; recurrence {f['finder'].get('intervals_h')}; particle coincident {f['coincident_particle_flare']}.",
           f"**Verdict: {f['classification']}**.",""])
        if "literature" in f:txt.extend(["References: "+", ".join(f["literature"]),""])
    txt+=["## Referee","",
      "**KEEP PR #205 DRAFT; DO NOT MERGE; NO [REMARKABLE].** Failed or unavailable source products are not negative QPE observations. Unmatched per-exposure PN time series are not complete-obs screens. An independent matched-energy soft background event-region check, SAS epatplot, held-out controls, event-list Poisson injections, trial-adjusted FPR, full 4LAC/SIMBAD+Legacy host vet and eRASS1/Swift/Chandra/TNS/ATel vetting of any flagged candidate remain open.","",
      "All request/download SHA-256 and retries are in the primary and recovery JSON receipts. Machine result: heasarc_postrescue_referee.json.",""]
    (OUT/"HEASARC_POSTRESCUE_REFEREE.md").write_text("\n".join(txt))
    print("HEASARC_POSTRESCUE",json.dumps({k:summary[k] for k in (
       "verified_full_merged_pn_source_curves","partial_single_exposure_pn_source_curves",
       "preliminary_recurrent_triggers","validated_new_qpes")}),flush=True)

if __name__=="__main__":
    main()
