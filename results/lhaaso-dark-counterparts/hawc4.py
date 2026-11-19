"""Freeze and parse the 2026 4HWC 85-source table from its public arXiv rendering.

The raw HTML is compressed and SHA-256 receipted. Name, sky position and
angular extension are parsed from the authors' *Table 1*; 4HWC's proposed
physical associations are precedence, not new identifications.
"""
from pathlib import Path
import gzip
import hashlib
import json
import re
import io
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from astropy.coordinates import SkyCoord
import astropy.units as u

OUT=Path(__file__).resolve().parent/"output"
CAT=OUT/"catalogs"
FROZEN=OUT/"frozen"
for p in (OUT,CAT,FROZEN):p.mkdir(parents=True,exist_ok=True)
URL="https://ar5iv.labs.arxiv.org/html/2602.00263"
KEY="hawc4_arxiv_2602_00263"

def first_number(s):
    v=re.search(r"[-+]?\d+(?:\.\d+)?",str(s).replace("\u2212","-"))
    return float(v.group(0)) if v else float("nan")

def main():
    path=OUT/"receipts.json"
    receipts=json.loads(path.read_text())
    # Supersede the old HAWC public-catalogue placeholder with 4HWC.
    receipts.pop("hawc_latest",None)
    status={"status":"pending","source":"Alfaro et al. 2026 arXiv:2602.00263 Table 1","url":URL}
    try:
        r=requests.get(URL,timeout=(12,60),
                       headers={"User-Agent":"Open-Space public catalogue research (arxiv cite)"})
        r.raise_for_status()
        raw=r.content
        if not 100_000<len(raw)<12_000_000:raise ValueError(f"unexpected HTML size {len(raw)}")
        html=BeautifulSoup(raw,"html.parser")
        parsed=[]
        for tr in html.find_all("tr"):
            cells=tr.find_all(["td","th"],recursive=False)
            if len(cells)<8:continue
            texts=[x.get_text(" ",strip=True) for x in cells]
            m=re.search(r"4HWC\s+J(\d{4}[+-]\d{4})",texts[0])
            if m is None:continue
            ra=first_number(texts[1])
            dec=first_number(texts[2])
            ext=first_number(texts[3])
            if not(np.isfinite(ra) and np.isfinite(dec) and 0<=ra<360 and abs(dec)<=90):
                raise ValueError("malformed position "+repr(texts[:4]))
            c=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,frame="icrs").galactic
            parsed.append({"name":"4HWC J"+m.group(1),"ra":ra,"dec":dec,
                           "l":float(c.l.deg),"b":float(c.b.deg),
                           "radius_deg":float(ext) if np.isfinite(ext) else 0.0,
                           "extension_model":"Gaussian_sigma" if np.isfinite(ext) else "point_source",
                           "ts":first_number(texts[-2]),"source":"Alfaro2026_Table1"})
        cat=pd.DataFrame(parsed).drop_duplicates(subset="name")
        if len(cat)!=85:
            raise ValueError(f"expected 85 catalogue sources from Table 1; got {len(cat)}")
        assert cat.ra.between(0,360).all() and cat.dec.between(-90,90).all()
        cat.to_csv(CAT/"hawc4_sources.csv",index=False)
        (FROZEN/"hawc4_2026_ar5iv.html.gz").write_bytes(gzip.compress(raw,compresslevel=8))
        status.update({"status":"loaded","retrieved_utc":datetime.now(timezone.utc).isoformat(),
                       "raw_sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),
                       "frozen_path":"frozen/hawc4_2026_ar5iv.html.gz",
                       "entries":len(cat),"arxiv_id":"2602.00263",
                       "survey_note":"TeV class veto, not proof of which accelerator powers a UHE source."})
    except Exception as exc:
        status["error"]=repr(exc)[:900]
    receipts[KEY]=status
    path.write_text(json.dumps(receipts,indent=2,sort_keys=True)+"\n")
    # Keep run_status consistent with changed receipts.
    p=OUT/"run_status.json"
    if p.exists():
        s=json.loads(p.read_text())
        s["loaded"]=sorted(k for k,v in receipts.items() if v.get("status")=="loaded")
        s["pending"]=sorted(k for k,v in receipts.items() if v.get("status")!="loaded")
        p.write_text(json.dumps(s,indent=2,sort_keys=True)+"\n")
    print("HAWC 4HWC 2026:",status["status"],status.get("entries",status.get("error","")))

if __name__=="__main__":main()
