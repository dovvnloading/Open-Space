#!/usr/bin/env python3
"""Network-free, independent frozen-receipt referee. No discovery by assertion."""
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
EXPECTED_ROWS=888152165
EXPECTED_INDEX="f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def num(r,k):
    try:
        f=float(r[k])
        return f if math.isfinite(f) else None
    except (ValueError,TypeError,KeyError):return None

def validate_row(r):
    j,je,h,he,y,ye,v,ve=[num(r,k) for k in (
       "flux_j_3fwhm_aper","fluxerr_j_3fwhm_aper",
       "flux_h_3fwhm_aper","fluxerr_h_3fwhm_aper",
       "flux_y_3fwhm_aper","fluxerr_y_3fwhm_aper",
       "flux_vis_3fwhm_aper","fluxerr_vis_3fwhm_aper")]
    if any(q is None for q in (j,je,he,ye,v,ve)):return False
    if min(j,je,he,ye,ve)<=0 or j/je<6-1e-7 or abs(v/ve)>=2:return False
    if num(r,"vis_det")!=0 or num(r,"flag_j")!=0 or num(r,"flag_h")!=0:return False
    if num(r,"ellipticity") is None or num(r,"ellipticity")>=.25:return False
    sp=num(r,"spurious_prob")
    if sp is not None and sp>=.5:return False
    yy=(-2.5*math.log10(y/j) if y is not None and y/ye>=2
        else -2.5*math.log10(2*ye/j))
    hh=(-2.5*math.log10(j/h) if h is not None and h/he>=2
        else -2.5*math.log10(j/(2*he)))
    return yy>=.35-1e-6 and hh<=.15+1e-6 and abs(yy-r["Y_minus_J_for_selection"])<1e-4 and abs(hh-r["J_minus_H_for_selection"])<1e-4

def main():
    summary=read(OUT/"summary.json")
    if summary["status"]!="SCAN_COMPLETE":raise RuntimeError("Not a full-release scan")
    controls=read(OUT/"controls_recovery.json")
    ctrl_counters={}
    for kind in ("cold_dwarf","warm_T_control","quasar"):
        c=[r for r in controls if r["class"]==kind]
        m=[r for r in c if r["catalog_matched_3arcsec"]]
        ctrl_counters[kind]=(len(c),len(m),sum(r["mask_recovered"] for r in m))
        declared=summary["control_recovery"][kind]
        if (declared["frozen_controls"],declared["matched_catalog_within_3arcsec"],
            declared["mask_recovered"])!=ctrl_counters[kind]:
            raise AssertionError("Control summary differs from receipts")
    if sha(ROOT/"selection.py")!=read(OUT/"batches/batch_000_manifest.json")["selection_sha256"]:
        raise RuntimeError("Science selection fingerprint differs")
    identifiers=set()
    partition_keys=set()
    measured=Counter()
    for n in range(32):
        m=read(OUT/"batches"/f"batch_{n:03d}_manifest.json")
        p=OUT/"batches"/f"batch_{n:03d}_candidates.jsonl.gz"
        if m["status"]!="COMPLETE" or m["batch"]!=n or m["index_sha256"]!=EXPECTED_INDEX:
            raise AssertionError("Bad partition manifest")
        if sha(p)!=m["candidates_sha256"]:
            raise AssertionError("Bad packet digest")
        local=0
        with gzip.open(p,"rt") as f:
            for line in f:
                r=json.loads(line)
                if not validate_row(r):raise AssertionError("Selected row fails independent mask")
                ident=str(r["object_id"])
                if ident in identifiers:raise AssertionError("Duplicate MER identity")
                identifiers.add(ident)
                local+=1
        if local!=m["candidate_count"]:raise AssertionError("Packet count mismatch")
        measured["strict"]+=local
        for part in m["partitions"]:
            k=(part["Norder"],part["Npix"])
            if k in partition_keys:raise AssertionError("Duplicate publisher HEALPix leaf")
            partition_keys.add(k)
            measured["rows"]+=part["rows"]
            measured["broad"]+=part["broad"]
            measured["partitions"]+=1
    if measured["rows"]!=EXPECTED_ROWS or measured["partitions"]!=768:
        raise AssertionError("No complete release coverage")
    if measured["strict"]!=summary["colour_mask_rows"] or measured["rows"]!=summary["scanned_rows"]:
        raise AssertionError("Summary not backed by frozen packets")
    shortlist=read(OUT/"unvetted_shortlist.json")
    for r in shortlist:
        if r["id"] not in identifiers:raise AssertionError("Untracked shortlist ID")
        if r["not_a_confirmed_cold_dwarf"] is not True:
            raise AssertionError("Unjustified discovery assertion")
    vet=read(OUT/"vetting_status.json")
    if vet["discovery_level"] or vet["new_objects_validated_colder_than_WISE_0855"]:
        raise AssertionError("Unsupported extraordinary claim")
    report={
        "status":"INDEPENDENT_RECEIPT_REPLAY_PASS_NO_DISCOVERY_CLAIM",
        "source_rows":measured["rows"],"partitions":measured["partitions"],
        "independently_rechecked_selected_rows":measured["strict"],
        "controls":ctrl_counters,"ranked_rows_checked":len(shortlist),
        "sha256_of_summary":sha(OUT/"summary.json"),
        "sha256_of_external_receipts":sha(OUT/"external_receipts.json"),
        "NISP_PSF_assessment":"NOT_CHECKED",
        "multi_epoch_image_artifact_rejection":"NOT_CHECKED",
        "unmasked_survey_area":"NOT_MEASURED",
        "spectroscopically_confirmed_novel_Y_dwarfs":0,
        "literature_completeness":"DATED_REFERENCE_SETS_ONLY",
        "discovery_level":False}
    (OUT/"referee.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("REFEREE",json.dumps(report),flush=True)

if __name__=="__main__":main()
