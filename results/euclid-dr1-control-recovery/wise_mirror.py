#!/usr/bin/env python3
"""Alternative public CatWISE2020 TAP catalog route and unWISE retry.
Never convert empty/error responses into flux limits or stationary objects.
Original IRSA receipts remain unchanged and separately citable.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import json,math,sys
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
sys.path.insert(0,str(ROOT.parent/"euclid-dr1-wide-highz-vetting"))
import vet

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")

def nl_catwise(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    h=.00084
    half=h/max(0.1,math.cos(math.radians(dec)))
    # All 16 top objects are far from RA=0 boundary; fail closed if not.
    if ra-half<0 or ra+half>360:
        return {"status":"RA_WRAP_REQUIRES_SEPARATE_QUERY","catalog_match":None}
    sql=("SELECT TOP 15 ra,dec,w1mpro,w2mpro,pmra,pmdec,sigpmra,sigpmdec "
         "FROM catwise2020.main "
         f"WHERE ra BETWEEN {ra-half:.8f} AND {ra+half:.8f} "
         f"AND dec BETWEEN {dec-h:.8f} AND {dec+h:.8f}")
    return vet.tap_request("https://datalab.noirlab.edu/tap/sync",sql,
                           "CatWISE2020_NOIRLab_mirror",("ra","dec"),r,timeout=18)

def unwise_gator(r):
    ra,dec=float(r["ra"]),float(r["dec"])
    url="https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    out={"status":"QUERY_ERROR","service":"unWISE2019_Gator_fallback",
         "endpoint":url,"catalog_match":None}
    try:
        import requests,hashlib
        v=requests.get(url,params={"catalog":"unwise_2019","spatial":"Cone",
                                   "objstr":f"{ra:.8f} {dec:.8f}",
                                   "radius":"3","radunits":"arcsec","outfmt":"1"},
                       timeout=(10,25))
        out["status_code"]=v.status_code
        out["response_sha256"]=hashlib.sha256(v.content).hexdigest()
        if v.status_code!=200 or "<html" in v.text.lower():
            raise RuntimeError("GATOR_HTTP_"+str(v.status_code))
        # Gator's IPAC table has | delimited schema and raw data; do not
        # guess positional columns or fluxes without a parsed table.
        from astropy.io import ascii
        t=ascii.read(v.text,format="ipac")
        out["returned_rows"]=len(t)
        out["status"]="CATALOG_ROWS_NOT_ASSOCIATED" if len(t) else "NO_ROWS_NOT_FLUX_LIMIT"
        matches=[]
        for row in t:
            d={k:str(row[k]) for k in t.colnames}
            mra=vet.good(d.get("ra"));mdec=vet.good(d.get("dec"))
            if mra is None or mdec is None:continue
            sep=vet.sphere_distance(ra,dec,mra,mdec)
            if sep<=3:matches.append((sep,d))
        if matches:
            sep,match=min(matches,key=lambda x:x[0])
            out["catalog_match"]=match
            out["separation_arcsec"]=round(sep,4)
            out["status"]="POSITION_MATCHED_CATALOG_FLUX_NOT_FORCED"
    except Exception as e:
        out["error"]=type(e).__name__+": "+str(e)[:260]
    return out

def one(row,i):
    cw=nl_catwise(row)
    uw=unwise_gator(row)
    q={"catwise_2020":cw,"unwise_2019":uw}
    diagnostics=vet.motion_and_wise(q)
    return {"rank":i+1,"id":row["id"],"CatWISE_NOIRLab":cw,
            "unWISE_Gator":uw,"diagnostics":diagnostics}

def main():
    top=json.loads((OUT/"top_candidates_with_checks.json").read_text())
    receipts=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        fs={pool.submit(one,row,i):i for i,row in enumerate(top[:16])}
        for f in as_completed(fs):
            i=fs[f]
            try:val=f.result()
            except Exception as e:
                val={"rank":i+1,"id":top[i]["id"],"status":"MIRROR_QUERY_EXCEPTION",
                     "error":type(e).__name__+": "+str(e)[:300]}
            receipts.append((i,val))
            print("WISE_MIRROR",i+1,val.get("diagnostics",{}).get("wise_color_status"),flush=True)
    receipts=[x for _,x in sorted(receipts)]
    save("wise_mirror_receipts.json",receipts)
    for i,rec in enumerate(receipts):
        top[i]["wise_mirror_receipt"]=rec
        d=rec.get("diagnostics") or {}
        if d.get("W1_minus_W2_vega_mag") is not None:
            top[i]["WISE_color_vega"]=d["W1_minus_W2_vega_mag"]
            top[i]["WISE_color_status"]="CATWISE2020_MIRROR_MATCH_NOT_FORCED_FLUX"
        elif not top[i].get("WISE_color_status"):
            top[i]["WISE_color_status"]="CATWISE_UNAVAILABLE_NO_FORCED_LIMIT"
        if d.get("foreground_motion_diagnostic","").startswith("REJECT"):
            top[i]["CatWISE_motion_diagnostic"]=d["foreground_motion_diagnostic"]
            top[i]["followup_priority"]="REJECT_FOREGROUND_MOTION"
    save("top_candidates_with_checks.json",top)
    status={"mirrors_attempted":len(receipts),
       "CatWISE2020_mirror_position_matched_W1_W2":sum(
          rec.get("diagnostics",{}).get("W1_minus_W2_vega_mag") is not None
          for rec in receipts),
       "CatWISE_significant_component_motion":sum(
          rec.get("diagnostics",{}).get("foreground_motion_diagnostic","").startswith("REJECT")
          for rec in receipts),
       "unWISE_fallback_position_matches":sum(
          rec.get("unWISE_Gator",{}).get("catalog_match") is not None for rec in receipts),
       "unWISE_missing_is_not_non_detection":True,
       "calibrated_midIR_forced_flux_limits":False}
    save("wise_mirror_summary.json",status)
    print("WISE_MIRROR_SUMMARY",json.dumps(status),flush=True)
if __name__=="__main__":main()
