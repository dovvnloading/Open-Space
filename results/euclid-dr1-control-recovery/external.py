#!/usr/bin/env python3
"""Bounded CatWISE/unWISE, Gaia, public literature/spectrum metadata receipts.
A metadata response or missing match is NOT a validated spectrum / flux limit.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,sys
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
sys.path.insert(0,str(ROOT.parent/"euclid-dr1-wide-highz-vetting"))
import vet

def save(name,x): (OUT/name).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+"\n")

def main():
    candidates=json.loads((OUT/"top_candidates.json").read_text())
    receipts=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        fut={pool.submit(vet.vet_one,r,i):i for i,r in enumerate(candidates[:16])}
        for f in as_completed(fut):
            i=fut[f]
            try:got=f.result()
            except Exception as e:
                got={"id":candidates[i]["id"],"status":"QUERY_EXCEPTION",
                     "error":type(e).__name__+": "+str(e)[:350],
                     "spectroscopic_redshift_confirmed_by_this_study":False}
            receipts.append((i,got))
            print("CATWISE_WISE_SPECTRUM",i+1,"/",len(fut),flush=True)
    receipts=[x for i,x in sorted(receipts)]
    pixel_file=OUT/"image_pixel_receipts.json"
    pixels=json.loads(pixel_file.read_text()) if pixel_file.exists() else []
    pix={r["id"]:r for r in pixels}
    by={r["id"]:r for r in receipts}
    enriched=[]
    for r in candidates:
        rec=by.get(r["id"])
        x=dict(r)
        x["external_vetting"]=rec
        x["WISE_color_vega"]=((rec or {}).get("diagnostics",{}).get("W1_minus_W2_vega_mag"))
        x["CatWISE_motion_diagnostic"]=((rec or {}).get("diagnostics",{}).get("foreground_motion_diagnostic","NOT_ESTABLISHED"))
        p=pix.get(r["id"])
        x["pixel_vetting_status"]=p.get("status") if p else "NOT_QUERIED"
        x["pixel_science_bands"]=p.get("science_bands_with_pixels",0) if p else 0
        x["image_artefacts_cleared"]=False
        x["new_spectrum_confirmed"]=False
        x["field_quasar_probability"]=None
        x["followup_priority"]=(("REJECT_FOREGROUND_MOTION" if x["CatWISE_motion_diagnostic"].startswith("REJECT")
                                 else "PHOTOMETRIC_PRIORITY_NEEDS_IMAGE_SPECTRUM") if rec else
                                "NOT_EXTERNALLY_VETTED")
        enriched.append(x)
    save("external_receipts.json",receipts)
    save("top_candidates_with_checks.json",enriched)
    stats={"targets_queried":len(receipts),
           "targets_with_any_service_error":sum(any(y.get("status")=="QUERY_ERROR"
                for y in r.get("queries",{}).values()) or r.get("status")=="QUERY_EXCEPTION"
                for r in receipts),
           "CatWISE_with_position_matched_W1_W2":sum(
                r.get("diagnostics",{}).get("W1_minus_W2_vega_mag") is not None
                for r in receipts),
           "significant_catalog_motion_flags":sum(
                r.get("diagnostics",{}).get("foreground_motion_diagnostic","").startswith("REJECT")
                for r in receipts),
           "actual_public_spectra_examined":0,
           "new_spectroscopic_confirmations":0,
           "missing_catalog_row_not_flux_limit":True,
           "simulated_photometric_quasar_probabilities":False,
           "status":"INCOMPLETE_EXTERNAL_VETTING_NO_NEW_SPECTROSCOPY"}
    save("external_summary.json",stats)
    print("EXTERNAL_SUMMARY",json.dumps(stats),flush=True)
if __name__=="__main__":main()
