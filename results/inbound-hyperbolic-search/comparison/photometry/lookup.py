"""Match the frozen pre-discovery Palomar positions to public ZTF quadrant metadata."""
import base64
import csv
from datetime import date
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import Request,urlopen

root=Path(__file__).resolve().parents[2]
out=Path(sys.argv[1])
out.mkdir(parents=True,exist_ok=True)
source=root/"orbits/data/2025-V1.obs.gz"
lines=gzip.decompress(source.read_bytes()).decode().splitlines()
metadata={}
receipts=[]
errors=[]
for row in [r for r in lines if r[77:80]=="I41"]:
    y,m,dd=row[15:32].split()
    jd=date(int(y),int(m),1).toordinal()+1721424.5+float(dd)-1
    h,minute,sec=map(float,row[32:44].split())
    ra=15*(h+minute/60+sec/3600)
    d,minute,sec=map(float,row[45:56].split())
    dec=(-1 if row[44]=="-" else 1)*(d+minute/60+sec/3600)
    url="https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?"+urlencode(
        {"POS":f"{ra:.9f},{dec:.9f}","WHERE":f"obsjd BETWEEN {jd-0.002:.9f} AND {jd+0.002:.9f}",
         "COLUMNS":"filefracday,field,ccdid,qid,filtercode,obsjd,infobits,exptime","ct":"csv"})
    receipt={"url":url,"mpc_mid_jd":jd,"ra_deg":ra,"dec_deg":dec}
    receipts.append(receipt)
    try:
        with urlopen(Request(url,headers={"User-Agent":"Open-Space-archive/1"}),timeout=45) as response:
            raw=response.read()
            receipt.update(sha256=hashlib.sha256(raw).hexdigest(),date=response.headers.get("Date"),
                           http_status=response.status,bytes=len(raw))
        import astropy.io.ascii as asc
        decoded=raw.decode("utf-8")
        table=csv.DictReader(io.StringIO(decoded)) if "filefracday," in decoded[:1200] else asc.read(decoded,format="ipac")
        for r in table:
            fields="filefracday field ccdid qid filtercode obsjd infobits exptime".split()
            item={k:str(r[k]).strip() for k in fields}
            mid=float(item["obsjd"])+float(item["exptime"])/172800
            if abs(jd-mid)*86400>8:
                continue
            item.update(pred_ra_deg=ra,pred_dec_deg=dec,matched_mpc_mid_jd=jd,reported_mag=row[65:70].strip())
            key=(item["filefracday"],item["ccdid"],item["qid"])
            metadata[key]=item
    except Exception as exc:
        receipt["error"]=str(exc)
        errors.append(str(exc))
res={"source":"MPC I41 positions already in the fit, used only as image controls",
     "positions":len(receipts),"query_failures":errors,"receipts":receipts,
     "ztf":{"ztf_nominal_footprints":list(metadata.values())}}
f=out/"archive-results.json"
f.write_text(json.dumps(res,indent=2)+"\n")
print("MATCH",len(metadata),"quadrants;",len(errors),"failures",flush=True)
print("RESULT","mpc_ztf_matches.json",base64.b64encode(f.read_bytes()).decode())
