#!/usr/bin/env python3
"""Trace slit-level provenance from DJA coadds. Never treat as native detector DQ."""
from pathlib import Path
import argparse,json,sys
import numpy as np
from astropy.io import fits
sys.path.insert(0,str(Path(__file__).resolve().parent))
import download_targets as down
HERE=Path(__file__).resolve().parent
GROUPS=json.loads((HERE/"literature_referee_2026-11-19.json").read_text())["entries"]
WORK=Path("/tmp/dja-pr197-nine-fits")
FIELDS=["filename","source_id","source_name","source_ra","source_dec","source_xpos","source_ypos",
       "shutter_state","shutter_id","slitlet_id","quadrant","detector","xstart","ystart","xsize",
       "ysize","msa_configuration_id","msa_metadata_id","nints","nod_type","x_offset","y_offset"]
def clean(x):
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,bytes):return x.decode("utf-8","replace").strip()
    if isinstance(x,float) and not np.isfinite(x):return None
    return x if isinstance(x,(str,bool,float,int)) or x is None else str(x)
def run(out):
    report={"scope":"Only 10 DJA coadded .spec.fits SLITS tables; no native integration DQ or complete MSA neighbouring shutters",
            "status":"METADATA_AUDIT_NOT_NATIVE_DQ_OR_CONTAMINATION_CLEARANCE","sources":[]}
    for g in GROUPS:
        if g["referee_class"]!="UNRESOLVED_NOT_NEW":continue
        for s in g["individual_spectra"]:
            path=WORK/s["root"]/s["file"]
            item={"group_id":g["group_id"],"file":s["file"]}
            if not path.exists():
                root,name,url=down.safe(s)
                down.download_one({"root":root,"file":name,"url":url},WORK,cap_mb=400)
            try:
                with fits.open(path,memmap=True) as h:
                    item["no_native_DQ"] = not ("DQ" in h and h["DQ"].data is not None)
                    slits=h["SLITS"].data
                    cols=h["SLITS"].columns.names
                    item["records"]=[{k:clean(row[k]) for k in FIELDS if k in cols}
                                     for row in slits]
                    item["unique_ids"]=sorted(set(str(r.get("source_id")) for r in item["records"]))
                    item["detectors"]=sorted(set(str(r.get("detector")) for r in item["records"]))
                    item["shutter_states"]=sorted(set(str(r.get("shutter_state")) for r in item["records"]))
                    item["no_conflicting_source_id_in_this_coadd"]=len(item["unique_ids"])<=1
                    item["potential_undetected_neighbour_slit"]=True
            except Exception as e:
                item["error"]=type(e).__name__+": "+str(e)[:250]
            report["sources"].append(item)
            out.mkdir(parents=True,exist_ok=True)
            (out/"followup_slits.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("SLIT_METADATA "+str(len(report["sources"]))+" spectra; NOT native integration DQ")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args();run(a.out)
