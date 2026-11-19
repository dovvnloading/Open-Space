#!/usr/bin/env python3
"""Freeze the Gaia-screened SDSS-only subset, not an archive QPE survey."""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
OUT=BASE/"outputs"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for data in iter(lambda:f.read(1024**2),b""):
            h.update(data)
    return h.hexdigest()

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg.get("stage")=="archive_finalize"
    frozen=json.loads((OUT/"archive_pn_point_freeze.json").read_text())
    sdss=json.loads((OUT/"archive_sdss_stellar_blazar_audit.json").read_text())
    gaia=json.loads((OUT/"archive_gaia_crossmatch_00000_00989.json").read_text())
    assert frozen["status"]=="PROVISIONAL_SAMPLE_FROZEN"
    assert sdss["status"]=="SDSS_SUBSET_EXCLUSIONS_COMPLETE"
    assert gaia["status"]=="COMPLETE" and not gaia["missing_chunks"]
    assert gaia["start"]==0 and gaia["stop"]==989
    parent=OUT/"archive_pn_point_freeze.csv"
    sdss_input=OUT/"archive_sdss_after_exclusions.csv"
    assert sha(parent)==frozen["manifest_sha256"]
    assert gaia["input_manifest_sha256"]==sha(sdss_input)
    with parent.open(newline="") as f: all_rows=list(csv.DictReader(f))
    with sdss_input.open(newline="") as f: sdss_rows=list(csv.DictReader(f))
    assert len(all_rows)==1108 and len(sdss_rows)==989
    assert len(set(x["detid"] for x in sdss_rows))==989
    gaia_flags={x["detid"] for x in gaia["stellar_flags"]}
    assert len(gaia_flags)==gaia["n_stellar_flagged_detections"]==22
    eligible=[r for r in sdss_rows if r["detid"] not in gaia_flags]
    assert len(eligible)==967
    outfile=OUT/"archive_provisional_nuclei_sdss_gaia.csv"
    with outfile.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(sdss_rows[0]))
        w.writeheader();w.writerows(eligible)
    ids=set(x["detid"] for x in eligible)
    assert not ids.intersection(gaia_flags)
    t={"as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "status":"PROVISIONAL_HOST_COHORT_FROZEN__NOT_QPE_SCREENED",
       "parent_4xmm_xray_soft_long":33488,
       "sdss_galaxy_centroid_matches":2510,
       "sdss_pointlike_longpn_detections":1108,
       "sdss_spectroscopic_qso_star_exclusions":119,
       "post_sdss_detections":989,
       "gaia_dr3_targets_audited":989,
       "gaia_dr3_counterparts":gaia["n_gaia_matches"],
       "gaia_dr3_astrometric_star_flags":22,
       "post_gaia_provisional_detections":len(eligible),
       "post_gaia_provisional_unique_xray_srcids":len(set(r["srcid"] for r in eligible)),
       "post_gaia_provisional_unique_obsids":len(set(r["obs_id"] for r in eligible)),
       "input_sdss_sha256":sha(sdss_input),
       "input_gaia_receipt_sha256":sha(OUT/"archive_gaia_crossmatch_00000_00989.json"),
       "final_csv_sha256":sha(outfile),
       "frozen_thresholds_sha256":sha(BASE/"control_thresholds.json"),
       "pn_source_lightcurves_screened":0,
       "pn_obsids_screened":0,
       "pn_obsids_attempted_before_gaia":3,
       "validated_new_qpes":0,
       "unvetted_recurrent_candidate_count":"NOT MEASURED (no pn curves)",
       "empirical_false_positive_rate":None,
       "archive_screen_permitted":False,
       "known_blazars_completely_excluded":False,
       "legacy_surveys_footprint_completely_screened":False,
       "independent_background_regions_tested":False,
       "candidate_literature_TNS_ATel_checked":False,
       "limitations":[
          "Only clean SDSS DR18 primary morphological galaxy hosts; Legacy DR10 complement unqueried",
          "Astrometric Gaia matches flagged at parallax >5sigma or component motion >5sigma are conservative contamination exclusions; spurious Gaia astrometry for extended optical objects requires adjudication",
          "SDSS QSO/star veto is not a full 4LAC/SIMBAD known-blazar or Gaia epoch-aware stellar-vetting census",
          "PPS new-source pn products unavailable via ESA XSA; three priority observations all disconnected, including repeat-control probe",
          "Source-specific independent background GTIs/FBKTSR, piled-event epatplot, Poisson event injection, negative AGN/stellar population and trial-corrected FPR unavailable",
          "Zero pn-screened obsIDs is not a null QPE detection and cannot support a rate or [REMARKABLE] finding"]}
    out=OUT/"archive_final_host_cohort.json"
    out.write_text(json.dumps(t,indent=2,allow_nan=False)+"\n")
    audit=OUT/"ARCHIVE_SCREEN_AUDIT.md"
    existing=audit.read_text()
    if "## Gaia DR3 post-SDSS astrometric exclusion" not in existing:
        extra=("\n## Gaia DR3 post-SDSS astrometric exclusion\n\n"+
          f"One complete 25-query uploaded Gaia DR3 positional audit checked **989/989** SDSS-surviving provisional detections within 1 arcsec; found **{t['gaia_dr3_counterparts']}** Gaia associations and **22** distinct detections with Gaia parallax/proper-motion significance >=5σ. After flag-exclusion, **{len(eligible)}** provisional host detections, **{t['post_gaia_provisional_unique_xray_srcids']}** unique X-ray SRCIDs and **{t['post_gaia_provisional_unique_obsids']}** observations remain. This is **not a spectroscopically verified nucleus/blazar-clean all-sky survey**, nor does Gaia prove all other candidates lack stellar contaminants. Audit: [Gaia raw checks](archive_gaia_crossmatch_00000_00989.json), [final provisional CSV](archive_provisional_nuclei_sdss_gaia.csv), [machine gate](archive_final_host_cohort.json).\n\n"+
          "**pn-source screen remains 0**, because XSA disconnected for all attempted new-source observations; no candidate counts can be inferred. **Referee: KEEP DRAFT / NO MERGE**.\n")
        audit.write_text(existing.rstrip()+"\n"+extra)
    print(json.dumps({k:t[k] for k in (
        "status","post_gaia_provisional_detections",
        "post_gaia_provisional_unique_xray_srcids","post_gaia_provisional_unique_obsids",
        "pn_source_lightcurves_screened","validated_new_qpes")},indent=2),flush=True)

if __name__=="__main__":main()
