#!/usr/bin/env python3
"""Recheck and persist transparent 300s CSV source/offsource curves for every flag.

These products are derivatives of verified HEASARC compressed FITS inputs, not
replacements for the archival FITS. The science finder itself is unchanged.
"""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
import control_validation as cv
from heasarc_batch import OUT, get_file, net_merged, save

def main():
    path=OUT/"heasarc_postrescue_referee.json"
    if not path.exists():path=OUT/"heasarc_aggregate_referee.json"
    summary=json.loads(path.read_text())
    flags=summary.get("triggers",summary.get("trigger_details",[]))
    receipts=[]
    for x in flags:
        oid=x["obs_id"];detid=x["detid"]
        product=x.get("source_product",x.get("pn_file",x.get("product")))
        expected=x.get("product_sha256",x.get("pn_fits_sha256"))
        rec={"obs_id":oid,"detid":detid,"source_product":product,"expected_fits_sha256":expected,
             "downloads":[]}
        if not product:
            rec["status"]="no_qualified_FITS_filename";receipts.append(rec);continue
        raw=get_file(oid,product,rec["downloads"],attempts=3)
        if raw is None:
            rec["status"]="download_failed";receipts.append(rec);continue
        rec["actual_fits_sha256"]=hashlib.sha256(raw).hexdigest()
        if not expected:
            rec["status"]="missing_original_sha256_not_compared";receipts.append(rec);continue
        if expected!=rec["actual_fits_sha256"]:
            rec["status"]="FITS_HASH_MISMATCH_NO_CSV";receipts.append(rec);continue
        try:
            pos,curve,ref,hdr=net_merged(raw)
            rr=np.asarray(curve["rate"],float)
            eb=np.asarray(curve["error"],float)
            br=np.asarray(curve["background"],float)
            ex=np.asarray(curve["exposure"],float)
            t=np.asarray(curve["t_s"],float)
            hard=curve.get("hard")
            filename=f"heasarc_trigger_curve_{oid}_{detid}.csv"
            file=OUT/filename
            with file.open("w",newline="") as fp:
                writer=csv.writer(fp)
                writer.writerow(["mjd_tt","time_s_since_mjdref","net_0p2_2keV_cps",
                   "net_error_cps","scaled_offsource_background_cps","hard_2_12keV_cps",
                   "effective_exposure_s"])
                for i in np.flatnonzero(curve["good"]):
                    if not all(np.isfinite(z) for z in (t[i],rr[i],eb[i],br[i],ex[i])):
                        continue
                    writer.writerow([round(float(ref+t[i]/86400),8),
                       round(float(t[i]),3),round(float(rr[i]),8),
                       round(float(eb[i]),8),round(float(br[i]),8),
                       round(float(hard[i]),8) if hard is not None and np.isfinite(hard[i]) else "",
                       round(float(ex[i]),3)])
            rec["status"]="FITS_HASH_VERIFIED_300S_CURVE_SAVED"
            rec["curve_csv"]=filename
            rec["csv_sha256"]=hashlib.sha256(file.read_bytes()).hexdigest()
            rec["good_300s_bins"]=int(np.count_nonzero(curve["good"]))
            rec["source_position_in_fits_deg"]=list(pos)
        except Exception as exc:
            rec["status"]="FITS_CURVE_ERROR";rec["error"]=repr(exc)[:300]
        receipts.append(rec)
    save("heasarc_trigger_curve_receipts.json",{
      "trigger_count":len(flags),"derived_curves_verified":sum(
        x["status"]=="FITS_HASH_VERIFIED_300S_CURVE_SAVED" for x in receipts),
      "note":"300s NET RATE1-3 + scaled PPS BACK1-3; not raw event data or independent soft event extraction",
      "sources":receipts})
    print("HEASARC_TRIGGER_CURVES",sum(x["status"]=="FITS_HASH_VERIFIED_300S_CURVE_SAVED" for x in receipts),"/",len(flags),flush=True)

if __name__=="__main__":main()
