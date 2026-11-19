#!/usr/bin/env python3
"""Independent integrity, coverage, and referee aggregation of frozen HEASARC pn batches."""
import csv, hashlib, json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs"
def dump(path,obj):
    (OUT/path).write_text(json.dumps(obj,indent=2,sort_keys=True,default=str,allow_nan=False)+"\n")

def main():
    manifest=OUT/"archive_provisional_nuclei_sdss_gaia.csv"
    sha=hashlib.sha256(manifest.read_bytes()).hexdigest()
    freeze=json.loads((OUT/"archive_final_host_cohort.json").read_text())
    assert sha==freeze["final_csv_sha256"],"Frozen host sample changed"
    rows=list(csv.DictReader(manifest.open(newline="")))
    expected={(r["obs_id"],r["detid"]) for r in rows}
    obsids={r["obs_id"] for r in rows}
    assert len(rows)==len(expected)==967 and len(obsids)==670
    files=sorted(OUT.glob("heasarc_batch_cohort_*.json"))
    seen=set(); screened=[]; flags=[]; errors=[]; batches=[]; missing=[]
    statuses=Counter(); background=Counter(); tested=set(); trials=0; recovered=0
    threshold=None
    for path in files:
        r=json.loads(path.read_text())
        o=r.get("observations",[])
        batches.append({"file":path.name,"start":r.get("start"),"requested":r.get("count_requested"),
                        "attempted":len(o),"status":r.get("status")})
        if r.get("cohort_manifest_sha256")!=sha:
            errors.append(path.name+": cohort hash mismatch");continue
        if threshold is None: threshold=r.get("thresholds_sha256")
        elif threshold!=r.get("thresholds_sha256"):
            errors.append(path.name+": finder threshold hash mismatch");continue
        if r.get("status")!="BATCH_FINISHED_UNVETTED" or len(o)!=r.get("count_selected"):
            errors.append(path.name+": incomplete batch")
        for obs in o:
            oid=obs["obs_id"]
            if oid not in obsids or oid in seen:
                errors.append(path.name+": duplicate or foreign obsid "+oid);continue
            seen.add(oid)
            statuses[obs.get("status","unknown")]+=1
            n=0
            for s in obs.get("source_results",[]):
                key=(oid,s["detid"])
                if key not in expected or key in tested:
                    errors.append(path.name+": bad or repeated source "+str(key));continue
                tested.add(key)
                if s.get("status")!="position_verified_net_pn_soft_screened":continue
                screened.append((oid,s));n+=1
                b=s.get("independent_field_particle_background",{})
                background[b.get("status","unavailable")]+=1
                i=s.get("poisson_count_domain_injection",{})
                trials+=i.get("n_trials",0);recovered+=i.get("n_recovered",0)
                if not s.get("unvetted_repeated_peak_trigger"):continue
                found=s.get("finder",{})
                x={"obs_id":oid,"detid":s["detid"],"srcid":s.get("srcid"),
                   "ra_deg":s.get("ra"),"dec_deg":s.get("dec"),
                   "pn_file":s.get("pn_product"),"pn_fits_sha256":s.get("merged_pn_fits_sha256"),
                   "position_error_arcsec":s.get("source_sep_arcsec"),
                   "pn_pileup_catalogue_diagnostic":s.get("pn_pileup_diagnostic"),
                   "finder":found,"field_particle":b,
                   "field_particle_peaks_coincident":sum(p.get("source_masked_particle_flare_coincident") is True for p in b.get("peak_checks",[])),
                   "source_specific_offsource_background":s.get("offsource_scaled_background_control"),
                   "count_domain_surrogate":i,
                   "pileup_epatplot":"NOT DONE","soft_offsource_event_extraction":"NOT DONE",
                   "event_list_injection":"NOT DONE","temperature_evolution":"NOT MEASURED",
                   "simbad_4LAC":"NOT FULLY CHECKED","eRASS1_TNS_ATel":"NOT FULLY CHECKED",
                   "validated_new_QPE":False,
                   "classification":"UNVETTED_FLAG_NOT_QPE_CANDIDATE"}
                if oid=="0144570101" and s["detid"]=="101445701010002":
                    x.update(classification="KNOWN_SEYFERT2_AGN_AND_PARTICLE_FLARE_NOT_NOVEL_QPE",
                             published_id="XMMU J134736.6+173403 / SDSS J134736.39+173404.6",
                             prior_literature=[
                               "https://arxiv.org/abs/0801.4246",
                               "https://arxiv.org/abs/1803.10976",
                               "https://doi.org/10.1051/0004-6361/202142444"],
                             swift_chandra_xmm_history="Swift 29 epochs 2008 (published 23.82 and 71.44 h QPO), Chandra 2008 source location, XMM 2003 and 2010-2012; not a novel QPE",
                             literature_status="reviewed_AGN_not_QPE")
                flags.append(x)
            if not n:
                missing.append({"obs_id":oid,"reason":obs.get("status"),
                                "merged_products":obs.get("pn_merged_srctsr8_available"),
                                "uninspected_products":obs.get("possible_uninspected_source_products")})
    gaia=json.loads((OUT/"archive_gaia_crossmatch_00000_00989.json").read_text())
    for f in flags:
        f["gaia_dr3_counterparts"]=[s for s in gaia["source_rows"] if s["detid"]==f["detid"]]
        f["gaia_astrometric_star_flag"]=any(s["detid"]==f["detid"] for s in gaia["stellar_flags"])
    screen_obs={obs for obs,_ in screened}
    audit={"timestamp_utc":datetime.now(timezone.utc).isoformat(),
           "cohort_sha256":sha,"finder_threshold_sha256":threshold,
           "frozen_cohort_obsids":670,"frozen_detections":967,
           "obsids_attempted":len(seen),"attempt_coverage_complete":len(seen)==670 and not errors,
           "obsids_with_valid_pn_curve":len(screen_obs),"pn_lightcurves_screened":len(screened),
           "detection_rows_not_screened":967-len(screened),
           "obsids_without_valid_pn_curve":670-len(screen_obs),
           "obsid_statuses":dict(statuses),"field_particle_statuses":dict(background),
           "offsource_scaled_background_flags":sum(s.get("offsource_scaled_background_control",{}).get("recurrence_unvetted") is True for _,s in screened),
           "poisson_count_surrogate_trials":trials,"poisson_count_surrogate_train_recoveries":recovered,
           "unvetted_recurrent_triggers":len(flags),"trigger_details":flags,
           "validated_new_qpes":0,"empirical_archive_FPR":None,"event_list_injection_complete":False,
           "soft_band_offsource_event_regions_complete":False,"SAS_epatplot_complete":False,
           "all_SIMBAD_4LAC_multi_epoch_TNS_ATel_vetting_complete":False,
           "ready_for_review":False,"merge_allowed":False,
           "batches":batches,"missing_valid_curve_obsids":missing,
           "unattempted_obsids":sorted(obsids-seen),"errors":errors,
           "referee":"REQUEST_CHANGES: NASA PPS retrieval works but survey still has missing/nonverified curves and no event-level background/pileup/FPR and complete novelty gates. KEEP PR DRAFT, no [REMARKABLE]."}
    dump("heasarc_aggregate_referee.json",audit)
    md=["# NASA HEASARC QPE frozen-archive referee audit","",
        f"**Attempted:** {len(seen)}/670 observations. **Verified:** {len(screen_obs)} observations, **{len(screened)}** source EPIC-pn 0.2–2 keV light curves of 967 eligible detection rows.",
        f"**Preliminary flags:** {len(flags)}; **validated novel QPEs: 0**; no [REMARKABLE].",
        f"**Independent source-masked field particle-curve comparisons:** {background.get('compared_source_masked_particle_flare_field',0)}. **Scaled off-source aperture repeat flags:** {audit['offsource_scaled_background_flags']}.",
        f"**Poisson count-domain proxy recoveries:** {recovered}/{trials}. Not a spatial/event-list injection, independent population FPR, or validated sensitivity.",
        f"**Missing verified pn curves:** {audit['obsids_without_valid_pn_curve']} observation IDs / {audit['detection_rows_not_screened']} detection rows. No missing product is a negative.",
        "", "## Batch receipts", "", "| Batch receipt | Attempted | Status |","|---|---:|---|"]
    for b in batches:
        md.append(f"| [{b['file']}]({b['file']}) | {b['attempted']}/{b['requested']} | {b['status']} |")
    md+=["","## Trigger-by-trigger assessment",""]
    for f in flags:
        md += [f"### {f['obs_id']} / {f['detid']}",
               f"FITS-position-verified at ({f['ra_deg']}, {f['dec_deg']}); recurrence {f['finder'].get('intervals_h')} hours, field particle-flare coincidences {f['field_particle_peaks_coincident']}.",
               f"**Disposition:** {f['classification']}. Gaia DR3 counterpart rows: {len(f['gaia_dr3_counterparts'])}; significant astrometric star flag: {f['gaia_astrometric_star_flag']}.",
               "SAS pile-up/soft off-source event-list re-extraction, full SIMBAD/4LAC and eRASS1/TNS/ATel are outstanding unless separately documented. QPE vs AGN/TDE remains unconfirmed for any unknown flags.",""]
        for u in f.get("prior_literature",[]): md.append("- Published study: "+u)
    md+=["","## Referee decision","",
         "**KEEP DRAFT / DO NOT MERGE.** Verified HEASARC pn source curves are observational evidence; unverified source rows are not QPE nulls. A source-masked particle-flare curve and scaled source BACK1–3 are useful controls, but not a full independent soft-band event-region background check. Count-domain injections do not provide trial-adjusted false-alarm probability or event-list completeness. No validated novel QPE and no remarkable entry.",""]
    (OUT/"HEASARC_REFEREE.md").write_text("\n".join(md))
    print("HEASARC_AGGREGATE",json.dumps({k:audit[k] for k in ["obsids_attempted","obsids_with_valid_pn_curve","pn_lightcurves_screened","unvetted_recurrent_triggers","validated_new_qpes"]}),flush=True)

if __name__=="__main__":
    main()
